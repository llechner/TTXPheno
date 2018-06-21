#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
ROOT.gROOT.SetBatch(True)
from math                                import sqrt, cos, sin, pi, isnan, sinh, cosh
import copy
import imp
import numpy as np

# RootTools
from RootTools.core.standard             import *

# TTXPheno
from TTXPheno.Tools.user                 import plot_directory
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR 
from TTXPheno.Tools.WeightInfo           import WeightInfo
from TTXPheno.Tools.cutInterpreter       import cutInterpreter

# Import samples
from TTXPheno.samples.benchmarks         import *

# Import helpers
from plot_helpers                        import *

# Import process variables
import process_variables

# Import additional
from array                               import array

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--plot_directory',     action='store',      default = 'gen')
argParser.add_argument('--sample',             action='store',      default = 'fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--process',            action='store',      default = 'ttZ')
argParser.add_argument('--order',              action='store',      default = 2)
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--selectPlots',        action='store',      default = [], type=int, nargs='+', help='Run only on a small subset of the data?')
argParser.add_argument('--selection',          action='store',      default = 'lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--variables',          action='store',      default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--parameters',         action='store',      default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150)

args = argParser.parse_args()

if len(args.parameters) == 0: args.parameters = None
if len(args.selectPlots) == 0: args.selectPlots = None

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample = getattr( samples, args.sample )

# Scale the plots with number of events used (implemented in ref_lumiweight1fb)
event_factor = 1.
fisher_directory = 'fisher_information'
if args.small:
    sample.reduceFiles( to = 10 )
    event_factor = sample.nEvents / float(sample.chain.GetEntries())
    fisher_directory += '_small'

# Polynomial parametrization
w = WeightInfo(sample.reweight_pkl)
w.set_order(int(args.order))
if len(args.variables) == 0: args.variables = w.variables
else: args.variables = [ item for item in w.variables if item in args.variables ]

# Format input parameters to dict
WC_string = 'SM'
WC = {}
if args.parameters is not None:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals   = list( map( float, str_vals ) )
    WC = { coeff:vals[i] for i, coeff in enumerate(coeffs) }
    WC_string = '_'.join( args.parameters )

def get_reweight_function():
    ''' return a weight function
    '''

    def reweight( event, sample ):
        return event.ref_lumiweight1fb * args.luminosity * event_factor

    return reweight

weightFunction = get_reweight_function()
weightString = 'ref_lumiweight1fb*%s*%s' %( str(args.luminosity), str(event_factor) )
selectionString = cutInterpreter.cutString( args.selection )

# Additional plot variables from process_variables_ROC.py (defined for ttZ, ttW and ttgamma)
plotVariables = getattr( process_variables, args.process )['2D']

# reduce plotting variables
if args.selectPlots is not None:
    plotVariables = [ var for var in plotVariables if var['index'] in args.selectPlots ]

if len(plotVariables) == 0:
    print( 'Variable indices not found. No plot variable selected! Exiting...' )
    exit()

expo = 1. / len(args.variables)

def appendPlotInformation( VariableList ):
    ''' append normalized detI and x value to dict in VariableList
    '''

    for plotVariable in VariableList:

        # check if selection string already contains cuts on this parameter
        selectionList = [ item for item in selectionString.split('&&') if plotVariable['var'].replace('abs(','').replace(')','') in item ]
        plotVariable['pre-selection'] = ','.join(selectionList).replace( plotVariable['var'], plotVariable['plotstring'] ) if len(selectionList) != 0 else ''

        #remove initial selection string
        sample.setSelectionString('1')

        coeffList = get2DCoeffPlotFromDraw( sample, args.order, plotVariable['var'], plotVariable['binning'], selectionString, weightString=weightString )

        detIList  = [ np.linalg.det( w.get_fisherInformation_matrix( coeffs, args.variables, **WC )[1] ) for coeffs in coeffList ]
#        detI0     = sum( detIList )
        detIList  = [ abs(detI)**expo for detI in detIList ]

        def binningToXList( binningList ):
            bins, start, stop = tuple( binningList )
            xList = list( np.linspace( start=start, stop=stop, num=bins+1 ) )
            return [ 0.5*(x + xList[i+1]) for i,x in enumerate( xList[:-1] ) ]

        xRange = binningToXList( plotVariable['binning'] )
        plotVariable['x_graph'] = array( 'd', xRange )
        plotVariable['y_graph'] = array( 'd', detIList )


# Additional plot variables from process_variables.py (defined for ttZ, ttW and ttgamma)
appendPlotInformation( plotVariables )

# Plots
def drawPlot( variable, log = False ): 
    ''' Plotting function
    '''

    # Canvas
    c1 = ROOT.TCanvas()

    mg = ROOT.TGraph( len(variable['x_graph']), variable['x_graph'], variable['y_graph'] )
    mg.SetFillColor( 30 )
    ROOT.gStyle.SetBarWidth( 1.2 )
    mg.Draw('AB')

    # Layout
    mg.GetXaxis().SetLabelSize(0.035)
    mg.GetXaxis().SetTitleSize(0.035)
    mg.GetXaxis().SetTitle( variable['plotstring'] )
    mg.GetXaxis().SetTitleOffset(1.3)
    mg.GetYaxis().SetTitleOffset(1.8)

    mg.GetYaxis().SetTitle( 'det(I_{ij})^{(1/%i)}'%len(args.variables) if len(args.variables) > 1 else 'det(I_{ij})' )
    mg.GetYaxis().SetLabelSize(0.035)
    mg.GetYaxis().SetTitleSize(0.035)
    mg.GetYaxis().SetRangeUser(0, max(variable['y_graph'])*1.4)

    if log:
        c1.SetLogy()
        mg.GetYaxis().SetRangeUser(1e-5, max(variable['y_graph'])*10)

    # Info
    t1 = ROOT.TPaveText(0.55, .80, .95, .92, "blNDC")
    t1.SetTextAlign(33)
    t1.SetTextSize(0.035)
    t1.SetTextFont(1)
    t1.SetBorderSize(0)
    t1.SetFillColor(ROOT.kWhite)
    if variable['pre-selection'] != '': t1.AddText('pre-selection on var: ' + variable['pre-selection'])
    t1.AddText('Restricted to: ' + ', '.join(args.variables) if len([i for i, j in zip(args.variables, w.variables) if i != j]) > 0 else 'No Restriction to Variables' )
    t1.AddText('WC: ' + ' '.join(args.parameters) if args.parameters is not None else 'Standard Model' )
    t1.Draw()

    c1.Update()

    plot_directory_ = os.path.join(\
        plot_directory,
        args.plot_directory, 
        sample.name, 
        fisher_directory, 
        'kinematics', 
        args.selection,
        WC_string,
        '_'.join(args.variables) if len([i for i, j in zip(args.variables, w.variables) if i != j]) > 0 else 'all',
        'log' if log else 'lin')

    if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)
    c1.Print(os.path.join(plot_directory_, '%s.png'%variable['var']))
    
    del c1

# Plot lin and log plots for each variable with full and fps scaling (cut plots)
for plotVariable in plotVariables:
    for log in [ True, False ]:
        drawPlot( plotVariable, log )


