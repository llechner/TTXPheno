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
import process_variables_ROC

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
    sample.reduceFiles( to = 1 )
    event_factor = sample.nEvents / float(sample.chain.GetEntries())
    fisher_directory += '_small'

# Polynomial parametrization
w = WeightInfo(sample.reweight_pkl)
w.set_order(int(args.order))
if len(args.variables) == 0: args.variables = w.variables

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


# Make sure that weightString contains the same as weightFunction!!!
weightString = 'ref_lumiweight1fb*%s*%s' %( str(args.luminosity), str(event_factor) )
weightFunction = get_reweight_function()

selectionString = cutInterpreter.cutString( args.selection )

# Additional plot variables from process_variables_ROC.py (defined for ttZ, ttW and ttgamma)
plotVariablesBin = getattr( process_variables_ROC, args.process )['bin']

# reduce plotting variables
if args.selectPlots is not None:
    plotVariablesBin = [ item for item in plotVariablesBin if item['index'] in args.selectPlots ]

if len(plotVariablesBin) == 0:
    print( 'Variable indices not found. No plot variable selected! Exiting...' )
    exit()

# fps info for normalization to full-pre-selection
coeff_sel = getCoeffListFromEvents( sample, selectionString = selectionString, weightFunction = weightFunction )
detI0_sel = np.linalg.det( w.get_total_fisherInformation_matrix( coeff_sel, args.variables, **WC )[1] )

# full info for normalization to full-sample
coeff_full = getCoeffListFromEvents( sample, selectionString = None, weightFunction = weightFunction )
detI0_full = np.linalg.det( w.get_total_fisherInformation_matrix( coeff_full, args.variables, **WC )[1] )

def appendPlotInformation( VariableList ):
    ''' append normalized detI and x value to dict in VariableList
    '''

    for plotVariable in VariableList:

        plotVariable['detI'] = []
        plotVariable['x_graph'] = []
        plotVariable['legendText'] = []
    
        # check if selection string already contains cuts on this parameter
        selectionList = [ item for item in selectionString.split('&&') if plotVariable['var'].replace('abs(','').replace(')','') in item ]
        plotVariable['pre-selection'] = ','.join(selectionList).replace( plotVariable['var'], plotVariable['plotstring'] ) if len(selectionList) != 0 else ''

        # x label
        plotVariable['plotstring'] = 'Bins in %s' %plotVariable['plotstring']

        # range of x
        x_vals = list( np.linspace( start=plotVariable['plotrange'][1], stop=plotVariable['plotrange'][2], num=plotVariable['plotrange'][0] ) )

        # list of cuts on variable (each cut a TGraph)
        cut_vals = plotVariable['cutPerGraph'] if len(plotVariable['cutPerGraph'])>0 else [None]
    
        for k, val in enumerate(cut_vals):

            detList = []

            for i, x in enumerate(x_vals):
        
                #remove initial selection string
                sample.setSelectionString('1')

                # add additional graphs for various cut values to the bin plot, do the usual for the cut plot
                selectionString_ = selectionString if val is None else selectionString + '&&%s>=%s' %( plotVariable['var'], str(val) )

                # add additional graphs for various bin values to the cut plot, do the usual for the bin plot
                binning_ = [int(x)] + plotVariable['binning'][1:]
        
                coeff = get2DCoeffPlotFromDraw( sample, args.order, plotVariable['var'], binning_, selectionString = selectionString_, weightString = weightString ) 
        
                # check if all entries are non-zero (otherwise FI matrix is 0-dim object -> Error)
                if all( [ all( [ v==0. for v in item ] ) for item in coeff ] ): continue

                detList.append( np.linalg.det( w.get_total_fisherInformation_matrix( coeff, args.variables, **WC )[1] ) )
        
            plotVariable['detI'].append( detList )
            plotVariable['x_graph'].append( x_vals )
            plotVariable['legendText'].append( '%s >= %s'%( plotVariable['var'], str(val) ) )


# Additional plot variables from process_variables_ROC.py (defined for ttZ, ttW and ttgamma)
appendPlotInformation( plotVariablesBin )

# Plots
def drawPlot( variable, log = False, scaleFull = False ): 
    ''' Plotting function
    '''

    expo = 1. / len(args.variables)
    normalization_string = 'full' if scaleFull else 'fps'
    #normalization
    detI0 = detI0_full if scaleFull else detI0_sel
    y = [ [ abs( val / detI0 )**expo for val in yList ] for yList in variable['detI'] ]


    def getTGraph( n, x, y, color=40, legendText='' ):
        ''' Create a TGraph object
        '''

        gr = ROOT.TGraph( n, x, y )
        gr.SetLineColor( color )
        gr.SetLineWidth( 2 )
        
        legend.AddEntry( gr, legendText, "l" )

        return gr


    # Canvas
    c1 = ROOT.TCanvas()
    mg = ROOT.TMultiGraph()

    # Legend position and margins
    leftMargin = c1.GetLeftMargin()
    rightMargin = c1.GetRightMargin()
    topMargin = c1.GetTopMargin()+0.2
    bottomMargin = c1.GetBottomMargin()
    legendWidth = 0.28
    legendHeight = 0.05*len(y)

    legend = ROOT.TLegend(1-rightMargin-legendWidth, bottomMargin, 1-rightMargin, bottomMargin+legendHeight) #bottom right
#    legend = ROOT.TLegend(leftMargin, 1-topMargin+0.2-legendHeight, leftMargin+legendWidth, 1-topMargin+0.2) #top left

    legend.SetBorderSize(1)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextFont(1)
    legend.SetTextSize(0.035)

    # Add TGraphs to Canvas
    for i, yData in enumerate(y):
        if len(variable['color'])>i: color = variable['color'][i]
        else: color = 30
        mg.Add( getTGraph( len(yData), array('d', variable['x_graph'][i]), array('d', yData), color, variable['legendText'][i] ) )

    mg.Draw('AL')
    legend.Draw()

    # Layout
    mg.GetXaxis().SetLabelSize(0.035)
    mg.GetXaxis().SetTitleSize(0.035)
    mg.GetXaxis().SetTitle( variable['plotstring'] )
    mg.GetXaxis().SetTitleOffset(1.3)
    mg.GetYaxis().SetTitleOffset(1.8)

    mg.GetYaxis().SetTitle( '((det(I_{ij}) / det(I_{ij}^{%s}))^{(1/%s)}'%(normalization_string, str(len(args.variables))) if len(args.variables) > 1 else 'det(I_{ij}) / det(I_{ij}^{%s})'%normalization_string )
    mg.GetYaxis().SetLabelSize(0.035)
    mg.GetYaxis().SetTitleSize(0.035)

    # Plot ranges
    ymin = min( [ min( yList ) for yList in y ] )
    if ymin < 1e-5: ymin = 1e-5
    ymax = max( [ max( yList ) for yList in y ] ) * 1.2 #spare area for legend
    xmin = variable['plotrange'][1]
    xmax = variable['plotrange'][2]

    if log:
        ymax *= 7   #spare area for legend
        c1.SetLogy()

    mg.GetYaxis().SetRangeUser(ymin, ymax)
    mg.GetXaxis().SetRangeUser(xmin, xmax)

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

    # Directory
    plot_directory_ = os.path.join(\
        plot_directory,
        args.plot_directory, 
        sample.name, 
        fisher_directory, 
        args.selection,
        WC_string,
        '_'.join(args.variables) if len([i for i, j in zip(args.variables, w.variables) if i != j]) > 0 else 'all',
        'ROC', 
        'binning',
        'full' if scaleFull else 'pre-selection', 
        'log' if log else 'lin')

    if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)
    c1.Print(os.path.join(plot_directory_, '%s.png'%variable['var']))
    
    del c1

# Plot lin and log plots for each variable with full and fps scaling
for plotVariable in plotVariablesBin:
    for full in [ True, False ]:
        for log in [ True, False ]:
            drawPlot( plotVariable, log, scaleFull = full )
