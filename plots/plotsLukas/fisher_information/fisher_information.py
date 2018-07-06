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

# Import samples
from TTXPheno.samples.benchmarks         import *

# Import helpers
from TTXPheno.Tools.plot_helpers                        import *

# Import process variables
import process_variables

# Import additional
from array                               import array

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--version',            action='store',      default='v7')
argParser.add_argument('--level',              action='store',      default='gen',  nargs='?', choices=['reco', 'gen', 'genLep'], help='Which level of reconstruction? reco, gen, genLep')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--process',            action='store',      default='ttZ')
argParser.add_argument('--order',              action='store',      default=2)
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--variables',          action='store',      default = [], type=str, nargs='+', help = "argument variables")
argParser.add_argument('--parameters',         action='store',      default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150)
argParser.add_argument('--binThreshold',       action='store',      default=100)
argParser.add_argument('--fpsScaling',         action='store_true', help='Scale to full pre-selection')

args = argParser.parse_args()

plot_subdirectory = "%s_%s"%(args.level, args.version)

# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco': from TTXPheno.Tools.cutInterpreterReco   import cutInterpreter
elif args.level == 'genLep': from TTXPheno.Tools.cutInterpreterGenLep import cutInterpreter
else:                    from TTXPheno.Tools.cutInterpreterGen       import cutInterpreter

if len(args.parameters) == 0: args.parameters = None

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample = getattr( samples, args.sample )

fisher_directory = 'fisher_information'
if args.small:
    sample.reduceFiles( to = 10 )
    fisher_directory += '_small'

# Scale the plots with number of events used (implemented in ref_lumiweight1fb)
event_factor = sample.nEvents / float(sample.chain.GetEntries())

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


selection_string = cutInterpreter.cutString( args.selection )

# Make sure that weightString contains the same as weightFunction!!!
weightString = 'ref_lumiweight1fb*%s*%s' %( str(args.luminosity), str(event_factor) )
weightFunction = get_reweight_function()

# split selection string step by step
selections = []
if not args.fpsScaling:
    selectionElements = args.selection.split('-')
    for i, item in enumerate(selectionElements[:-1]):
        selections.append( {'plotstring':' + '.join(replace_selectionstrings(selectionElements[:i+1])) + ' selection', 'selection':'-'.join(selectionElements[:i+1])} )

selections.append( {'plotstring':'full pre-selection (fps)', 'selection':args.selection} )
    
# Additional plot variables from process_variables.py (defined for ttZ, ttW and ttgamma)
plotVariables2D = getattr( process_variables, args.process )['2D']
plotVariables3D = getattr( process_variables, args.process )['3D']
plotVariables4D = getattr( process_variables, args.process )['4D']

if args.level != 'gen':
    if args.level == 'reco': search_string, replacement_string = ( 'gen', 'reco' )
    # Take care if that is really everything you have to replace!!!
    elif args.level == 'genLep': search_string, replacement_string = ( 'genZ', 'genLepZ' )
    for item in plotVariables2D + plotVariables3D + plotVariables4D:
        item['var'] = item['var'].replace(search_string, replacement_string)


# Calculate coefficients for binned distribution
# Calculate determinant using the 'variables' submatrix of FI
for var in plotVariables2D:
    #remove initial selection string
    sample.setSelectionString('1')

    var['coeff']       = getCoeffPlotFromDraw( sample, args.order, var['var'], var['binning'], selection_string, weightString=weightString, nEventsThresh=args.binThreshold )
    # add bin information to plot labels
    var['plotstring'] = 'fps + ' + var['plotstring'] + ' (%s bins)' %str(var['binning'][0])
    var['color']       = 30

for var in plotVariables3D:
    #remove initial selection string
    sample.setSelectionString('1')

    var['coeff']       = get2DCoeffPlotFromDraw( sample, args.order, var['var'], var['binning'], selection_string, weightString=weightString, nEventsThresh=args.binThreshold )
    # add bin information to plot labels
    var['plotstring'] = 'fps + ' + var['plotstring'] + ' (%s:%s bins)' %(str(var['binning'][0]), str(var['binning'][3]))
    var['color']       = 41

for var in plotVariables4D:
    #remove initial selection string
    sample.setSelectionString('1')

    var['coeff']       = get3DCoeffPlotFromDraw( sample, args.order, var['var'], var['binning'], selection_string, weightString=weightString, nEventsThresh=args.binThreshold )
    # add bin information to plot labels
    var['plotstring'] = 'fps + ' + var['plotstring'] + ' (%s:%s:%s bins)' %(str(var['binning'][0]), str(var['binning'][3]), str(var['binning'][6]))
    var['color']       = 45

# Calculate coefficients unbinned (event loop)
# Calculate determinant using the 'variables' submatrix of FI
for selection in selections:
    sample.setSelectionString('1')
    selection['coeff'] = getCoeffListFromEvents( sample, selectionString = cutInterpreter.cutString(selection['selection']), weightFunction = weightFunction )
    selection['color'] = 46

# Full Fisher information
if not args.fpsScaling:
    sample.setSelectionString('1')
    full              = { 'plotstring':'full'}
    full['coeff']     = getCoeffListFromEvents( sample, selectionString = None, weightFunction = weightFunction )
    full['color']     = 15

expo = 1. / len(args.variables)
data = [full] if not args.fpsScaling else []
data += selections + plotVariables2D + plotVariables3D + plotVariables4D
n_data = len(data)

# Fill dictionaries with normalized fisher information and plotting data
for i,item in enumerate(data):
    # if number of events per bin < args.nEventThresh for all bins, the list will be empty
    if len(item['coeff']) == 0:
        norm_detI = 0
    else:
        detI            = np.linalg.det( w.get_total_fisherInformation_matrix( item['coeff'], args.variables, **WC )[1] )
        if i==0: detI0  = detI
        norm_detI       = abs( detI / detI0 )**expo if detI0 != 0 else 0
    item['x_graph'] = array( 'd', range( 1, n_data+1 ) )
    item['y_graph'] = array( 'd', [0]*i + [norm_detI] + [0]*(n_data-i) )

# Plots
def drawPlot( log = False ): 
    ''' Plotting function
    '''

    scalingLabel = 'full' if not args.fpsScaling else 'fps'

    def getTGraph( n, x, y, color=40 ):
        ''' Create a TGraph object
        '''

        gr = ROOT.TGraph( n, x, y )
        gr.SetFillColor( color )
        gr.SetLineWidth( 4 )

        return gr


    # Canvas
    c1 = ROOT.TCanvas()
    c1.SetBottomMargin(0.05)

    mg = ROOT.TMultiGraph()

    # Add a TGraph for each 'bin' (overkill but flexible in changing the color)
    for item in data:
        mg.Add( getTGraph( n_data, item['x_graph'], item['y_graph'], item['color'] ) )

    mg.Draw('AB')

    # Layout
    mg.GetXaxis().SetLabelSize(0)
    mg.GetXaxis().SetTickLength(0.)
    mg.GetXaxis().SetLabelOffset(999)
    mg.GetYaxis().SetTitle( '(det(I_{ij}) / det(I_{ij}^{%s}))^{(1/%s)}'%(scalingLabel, len(args.variables)) if len(args.variables)>1 else 'det(I_{ij}) / det(I_{ij}^{%s})'%scalingLabel )
    mg.GetYaxis().SetLabelSize(0.035)
    mg.GetYaxis().SetTitleSize(0.035)

    # Plot ranges
    ymin = 0
    ymax = max( [ max(item['y_graph']) for item in data ] ) * 1.2  #spare area for legend

    if log:
        ymin = min( [ min( [ y for y in list(item['y_graph'])+[10] if y > 1e-5 ] ) for item in data ] ) * 0.5
        ymax *= 7   #spare area for legend
        c1.SetLogy()

    mg.GetYaxis().SetRangeUser(ymin, ymax)

    c1.Update()

    # Info
    t1 = ROOT.TPaveText(0.55, .85, .95, .92, "blNDC")
    t1.SetTextAlign(33)
    t1.SetTextSize(0.035)
    t1.SetTextFont(1)
    t1.SetBorderSize(0)
    t1.SetFillColor(ROOT.kWhite)
    t1.AddText('Restricted to: ' + ', '.join(args.variables) if len([i for i, j in zip(args.variables, w.variables) if i != j]) > 0 else 'No Restriction to Variables' )
    t1.AddText('WC: ' + ' '.join(args.parameters) if args.parameters is not None else 'Standard Model' )
    t1.Draw()

    # Selection info
    t = ROOT.TLatex()
    t.SetTextAngle(90)
    t.SetTextSize(0.035)

    # Labeling
    for i, item in enumerate(data):
        t.DrawLatex( i+1.25, ymin*1.4 if log else 0.05, item['plotstring'] )

    # Directory
    plot_directory_ = os.path.join(\
        plot_directory,
        plot_subdirectory,
        sample.name,
        fisher_directory,
        'fisherinfo',
        args.selection,
        WC_string,
        '_'.join(args.variables) if len([i for i, j in zip(args.variables, w.variables) if i != j]) > 0 else 'all',
        scalingLabel,
        'log' if log else 'lin')
    
    if not os.path.isdir(plot_directory_): os.makedirs(plot_directory_)
    c1.Print(os.path.join(plot_directory_, 'fisherinfo_%sEventsPerBin.png'%str(args.binThreshold)))
    
    del c1


# Plot lin and log plots
for log in [ True, False ]:
    drawPlot( log )
