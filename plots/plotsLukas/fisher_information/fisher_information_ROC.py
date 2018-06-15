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
#argParser.add_argument('--scan_variable',      action='store',      default = 'Z_pt')
#argParser.add_argument('--scan_parameters',    action='store',      default = [0, 500, 20], help='list of [min, max, steps]')
argParser.add_argument('--order',              action='store',      default = 2)
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--selection',          action='store',      default = 'lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--variables',          action='store',      default = ['cpt', 'cpQM', 'ctZ', 'ctZI'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--parameters',         action='store',      default = ['ctW', '3', 'ctWI', '3', 'ctZ', '3', 'ctZI', '3'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150)

args = argParser.parse_args()

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

# Format input parameters to dict
WC_string = 'SM'
WC = {}
if len(args.parameters) != 0:
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

#coeff_sel = getCoeffListFromEvents( sample, selectionString = selectionString, weightFunction = weightFunction )
#detI_sel  = np.linalg.det( w.get_total_fisherInformation_matrix( coeff_sel, args.variables, **WC )[1] )

expo = 1. / len(args.variables)

# Additional plot variables from process_variables.py (defined for ttZ, ttW and ttgamma)
plotVariables = getattr( process_variables_ROC, args.process )

for plotVariable in plotVariables:
    plotVariable['norm_detI'] = []
    plotVariable['x_graph']   = []

    variableCut = list( np.linspace( start=plotVariable['plotrange'][1], stop=plotVariable['plotrange'][2], num=plotVariable['plotrange'][0]+1 ) )

    for i, cut in enumerate(variableCut[:-1]):
        sample.setSelectionString('1')
        coeff =  get2DCoeffPlotFromDraw( sample, args.order, plotVariable['var'], plotVariable['binning'], selectionString = selectionString + '&&%s>=%s' %( plotVariable['var'], cut ), weightString = weightString ) 
        detI =  np.linalg.det( w.get_total_fisherInformation_matrix( coeff, args.variables, **WC )[1] ) 
        if i==0: detI0 =  np.linalg.det( w.get_total_fisherInformation_matrix( coeff, args.variables, **WC )[1] ) 
        plotVariable['norm_detI'].append( abs( detI / detI0 )**expo if detI0 != 0 else 0 )
        plotVariable['x_graph'].append( cut )


# Plots
def drawPlot( variable, log = False ): 
    ''' Plotting function
    '''

    # Canvas
    c1 = ROOT.TCanvas()
#    c1.SetBottomMargin(0.05)

    mg = ROOT.TGraph( len(variable['x_graph']), array('d', variable['x_graph']), array('d', variable['norm_detI']) )
    mg.SetLineColor(30)
    mg.SetLineWidth(2)
    mg.Draw('AL')

    # Layout
    mg.GetXaxis().SetLabelSize(0.035)
    mg.GetXaxis().SetTitleSize(0.035)
    mg.GetXaxis().SetTitle( variable['plotstring'] )
    mg.GetXaxis().SetTitleOffset(1.3)
#    mg.GetXaxis().SetLabelSize(0)
#    mg.GetXaxis().SetTickLength(0.)
#    mg.GetXaxis().SetLabelOffset(999)
    mg.GetYaxis().SetTitleOffset(1.5)
    mg.GetYaxis().SetTitle( '((det(I_{ij}) / det(I_{ij}^{no-cut}))^{(1/%s)}'%len(args.variables) if len(args.variables)>1 else 'det(I_{ij}) / det(I_{ij}^{no-cut})' )
    mg.GetYaxis().SetLabelSize(0.035)
    mg.GetYaxis().SetTitleSize(0.035)
    mg.GetYaxis().SetRangeUser(0., 1.1)

    if log:
        mg.GetYaxis().SetRangeUser(.00001, 2.)
        c1.SetLogy()

    # Info
    t1 = ROOT.TPaveText(0.55, .85, .95, .92, "blNDC")
    t1.SetTextAlign(33)
    t1.SetTextSize(0.035)
    t1.SetTextFont(1)
    t1.SetBorderSize(0)
    t1.SetFillColor(ROOT.kWhite)
    t1.AddText('Restricted to: ' + ', '.join(args.variables))
    t1.AddText('WC: ' + ' '.join(args.parameters))
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
        '_'.join(args.variables),
        'ROC', 
        'log' if log else 'lin')
    
    if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)
    c1.Print(os.path.join(plot_directory_, '%s.png'%variable['var']))
    
    del c1


# Plot lin and log plots
for plotVariable in plotVariables:
    for log in [ True, False ]:
        drawPlot( plotVariable, log )
