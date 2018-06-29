#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
ROOT.gROOT.SetBatch(True)
from math                                import sqrt, cos, sin, pi, isnan, sinh, cosh
import copy
import imp

# RootTools
from RootTools.core.standard             import *

# TTXPheno
from TTXPheno.Tools.user                 import plot_directory
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR, nanJet, nanLepton, getObjDict
from TTXPheno.Tools.WeightInfo           import WeightInfo

# Import samples
from TTXPheno.samples.benchmarks         import *


# Import helpers
from plot_helpers                        import *
from plot_kinematics_helpers             import *

# Import process variables
import process_variables

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',   nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--version',            action='store',      default='v7',     help='Appendix to plot directory')
argParser.add_argument('--process',            action='store',      default = 'ttZ')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',              action='store',      default=2,        help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',              action='store_true',                   help='Run only on a small subset of the data?')
argParser.add_argument('--level',              action='store',      default='gen',    nargs='?', choices=['reco', 'gen', 'genLep'], help='Which level of reconstruction? reco, gen, genLep')
argParser.add_argument('--variables',          action='store',      default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--parameters',         action='store',      default = ['ctW', '3', 'ctWI', '3', 'ctZ', '3', 'ctZI', '3'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150,      help='Luminosity for weighting the plots')
argParser.add_argument('--selectPlots',        action='store',      default = [], type=int, nargs='+', help='Run only on a small subset of the data?')

args = argParser.parse_args()

plot_subdirectory = "%s_%s"%(args.level, args.version)
fisher_directory = 'fisher_information'

colors = [ ROOT.kOrange, ROOT.kOrange+10, ROOT.kBlue, ROOT.kBlue-2, ROOT.kCyan+1, ROOT.kGreen+1, ROOT.kRed, ROOT.kRed+2, ROOT.kViolet+2, ROOT.kYellow+2, ROOT.kRed-7, ROOT.kPink-7, ROOT.kPink-3, ROOT.kGreen+4, ROOT.kGray+2 ]

# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco':
    from TTXPheno.Tools.cutInterpreterReco   import cutInterpreter
elif args.level == 'genLep':
    from TTXPheno.Tools.cutInterpreterGenLep import cutInterpreter
else:
    from TTXPheno.Tools.cutInterpreterGen    import cutInterpreter

preTag = 'reco' if args.level == 'reco' else 'gen'
tag    = 'reco' if args.level == 'reco' else 'genLep'

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None )
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None )

# Make subdirectory
fisher_directory = 'fisher_information'
if args.small: fisher_directory += '_small'

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
samples     = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample      = getattr( samples, args.sample )

# Import process specific variables
process_file = os.path.join( os.path.dirname( os.path.realpath( __file__ ) ), 'addons', '%s.py'%args.process )
process      = imp.load_source( "process", os.path.expandvars( process_file ) )

if args.small: sample.reduceFiles( to = 1 )
# Scale the plots with number of events used (implemented in ref_lumiweight1fb)
event_factor = sample.nEvents / float(sample.chain.GetEntries())

WC_string = 'SM'
WC = {}
params = []
if args.parameters is not None:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals   = list( map( float, str_vals ) )
    for i_param, (coeff, val, str_val) in enumerate(zip(coeffs, vals, str_vals)):
        params.append( { 
            'legendText': ' '.join([coeff,str_val]),
            'WC'        : { coeff:val },
            'color'     : colors[i_param], 
         })
#    WC = { coeff:vals[i] for i, coeff in enumerate(coeffs) }
    WC_string = '_'.join( args.parameters )
params.append( {'legendText':'SM', 'WC':{}, 'color':ROOT.kBlack} )

# Make stack 
stack = Stack(*[ [ sample ] for param in params ] )

def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'data' if hasData else "Simulation (%s)"%args.level),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( float(args.luminosity), dataMCScale ) ) if hasData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % float(args.luminosity))
    ]
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots):
  for plot in plots:
    for i_h, h in enumerate(plot.histos):
      h[0].style = styles.lineStyle(params[i_h]['color'])

  for log in [False, True]:
    # Directory structure
    plot_directory_ = os.path.join(\
        plot_directory,
        plot_subdirectory,
        sample.name,
        fisher_directory,
        'kinematics',
        args.selection,
        WC_string,
        '_'.join(args.variables) if len([i for i, j in zip(args.variables, w.variables) if i != j]) > 0 else 'all',
        'log' if log else 'lin')

    # plot the legend
    l_plot = copy.deepcopy(plots[0])
    for i_h, h in enumerate(l_plot.histos):
      h[0].legendText = params[i_h]['legendText']
      h[0].style = styles.lineStyle(params[i_h]['color'])
      h[0].Scale(0.)
      h[0].GetXaxis().SetTickLength(0.)
      h[0].GetYaxis().SetTickLength(0.)
      h[0].GetXaxis().SetLabelOffset(999.)
      h[0].GetYaxis().SetLabelOffset(999.)
    l_plot.name = "legend"
    l_plot.texX = ''
    l_plot.texY = ''
    plotting.draw(l_plot,
        plot_directory = plot_directory_,
        ratio = None,
        logX = False, logY = log, sorting = True,
        legend =  ( (0.17,0.9-0.05*sum(map(len, l_plot.histos))/3,1.,0.9), 3),
        drawObjects = drawObjects( ),
        copyIndexPHP = True,
    )

    # plot the plots
    for plot in plots:
      for i_h, h in enumerate(plot.histos):
        h[0].legendText = params[i_h]['legendText']
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = None,
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0., "auto"),
        scaling = {i:(len(params)-1) for i in range(len(params)-1)} if args.scaleLumi else {}, #Scale BSM shapes to SM (last in list)
	    legend = ( (0.17,0.9-0.05*sum(map(len, plot.histos))/3,1.,0.9), 3),
	    drawObjects = drawObjects( ),
        copyIndexPHP = True,
      )

read_variables = process.getVariableList( args.level )

logger.info( "Translating cut %s to %s", args.selection, cutInterpreter.cutString(args.selection) )
sample.setSelectionString( cutInterpreter.cutString(args.selection) )
sample.style = styles.lineStyle(ROOT.kBlue)

weightString = 'ref_lumiweight1fb*%s*%s' %( str(args.luminosity), str(event_factor) )
weight = [ [ weightString ] for param in params ]
selectionString = cutInterpreter.cutString( args.selection )

# Additional plot variables from process_variables_ROC.py (defined for ttZ, ttW and ttgamma)
plotVariables = getattr( process_variables, args.process )['2D']

# reduce plotting variables
if args.selectPlots is not None:
    plotVariables = [ var for var in plotVariables if var['index'] in args.selectPlots ]

if len(plotVariables) == 0:
    print( 'Variable indices not found. No plot variable selected! Exiting...' )
    exit()

if args.level != 'gen':
    if args.level == 'reco': search_string, replacement_string = ( 'gen', 'reco' )
    # Take care if that is really everything you have to replace!!!
    elif args.level == 'genLep': search_string, replacement_string = ( 'genZ', 'genLepZ' )
    for item in plotVariables:
        item['var'] = item['var'].replace(search_string, replacement_string)

#sequence functions
sequence = [] #process.getSequenceList( args.level )

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, addOverFlowBin=None)

plots = [] #process.getPlotList( args.scaleLumi, args.level )

for plotVariable in plotVariables:
    histos = getFisherKinematicsHisto( sample, order=args.order, var=plotVariable['var'], plotstring=plotVariable['plotstring'], binning=plotVariable['binning'], selectionString=selectionString, weightString=weightString, variables=args.variables, parameterList=params )
    plot = Plot( name = plotVariable['var'],  texX = plotVariable['plotstring'], texY = 'det(I_{ij})', binning = plotVariable['binning'] )
    plots.append( plot.fromHisto( 'histo', histos[0] ) )

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

drawPlots(plots)

