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

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',   nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--version',            action='store',      default='v7',     help='Appendix to plot directory')
argParser.add_argument('--processFile',        action='store',      default='ttZ_3l', nargs='?', choices=['ttZ_3l', 'ttZ_4l', 'ttW_2l', 'ttgamma_1l', 'ttgamma_2l'], help='Which process? ttZ, ttW, ttgamma')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',              action='store',      default=2,        help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',          action='store',      default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',              action='store_true',                   help='Run only on a small subset of the data?')
argParser.add_argument('--backgrounds',        action='store_true',                   help='include backgrounds?')
argParser.add_argument('--level',              action='store',      default='gen',    nargs='?', choices=['reco', 'gen', 'genLep'], help='Which level of reconstruction? reco, gen, genLep')
argParser.add_argument('--scaleLumi',          action='store_true',                   help='Scale lumi only?')
argParser.add_argument('--reweightPtXToSM',    action='store_true',                   help='Reweight Pt(X) to the SM for all the signals?')
argParser.add_argument('--parameters',         action='store',      default = ['ctW', '3', 'ctWI', '3', 'ctZ', '3', 'ctZI', '3'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',         action='store',      default=150,      help='Luminosity for weighting the plots')

args = argParser.parse_args()

# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco':
    from TTXPheno.Tools.cutInterpreterReco   import cutInterpreter
    from TTXPheno.Tools.objectSelection      import isGoodRecoJet      as isGoodJet
    from TTXPheno.Tools.objectSelection      import isGoodRecoLepton   as isGoodLepton

elif args.level == 'genLep':
    from TTXPheno.Tools.cutInterpreterGenLep import cutInterpreter
    from TTXPheno.Tools.objectSelection      import isGoodGenJet       as isGoodJet
    from TTXPheno.Tools.objectSelection      import isGoodGenLepton    as isGoodLepton

#remove Gen here for getting the reas cutInterpreter (this is just to check with old plots)
else:
#    from TTXPheno.Tools.cutInterpreterGen    import cutInterpreter
    from TTXPheno.Tools.cutInterpreter       import cutInterpreter
    from TTXPheno.Tools.objectSelection      import isGoodGenJet       as isGoodJet
    from TTXPheno.Tools.objectSelection      import isGoodGenLepton    as isGoodLepton

preTag = 'reco' if args.level == 'reco' else 'gen'
tag    = 'reco' if args.level == 'reco' else 'genLep'

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None )
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None )

# Make subdirectory
subDirectory = []
if args.scaleLumi:       subDirectory.append("shape")
else:                    subDirectory.append("lumi")
if args.reweightPtXToSM: subDirectory.append("reweightPtXToSM")
if args.small:           subDirectory.append("small")
subDirectory = '_'.join( subDirectory )

# Format WC input parameters
colors = [ ROOT.kOrange, ROOT.kOrange+10, ROOT.kBlue, ROOT.kBlue-2, ROOT.kCyan+1, ROOT.kGreen+1, ROOT.kRed, ROOT.kRed+2, ROOT.kViolet+2, ROOT.kYellow+2, ROOT.kRed-7, ROOT.kPink-7, ROOT.kPink-3, ROOT.kGreen+4, ROOT.kGray+2 ]
coeffs = args.parameters[::2]
str_vals = args.parameters[1::2]
vals   = list( map( float, str_vals ) )
params = []
for i_param, (coeff, val, str_val) in enumerate(zip(coeffs, vals, str_vals)):
    params.append( { 
        'legendText': ' '.join([coeff,str_val]),
        'WC'        : { coeff:val },
        'color'     : colors[i_param],
        })

# Import process specific variables
process_file = os.path.join( os.path.dirname( os.path.realpath( __file__ ) ), 'addons', '%s.py'%args.processFile )
process      = imp.load_source( "process", os.path.expandvars( process_file ) )

# Import samples
sample_file   = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
loadedSamples = imp.load_source( "samples", os.path.expandvars( sample_file ) )

signalSample = getattr( loadedSamples, args.sample )
WZSample     = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights' )
ttSample     = getattr( loadedSamples, 'fwlite_tt_lep_LO_order2_15weights_ref' )

# Polynomial parametrization
w = WeightInfo(signalSample.reweight_pkl)
w.set_order(int(args.order))

# configure samples
for s in [ signalSample, WZSample, ttSample ]:
    s.setSelectionString( cutInterpreter.cutString(args.selection) )
    s.style = styles.lineStyle(ROOT.kBlue)
    if args.small: s.reduceFiles( to = 5 )
    # Scale the plots with number of events used (implemented in ref_lumiweight1fb)
    s.event_factor = s.nEvents / float( s.chain.GetEntries() )

signal = [ signalSample for param in params ]
bg     = []

#for param in params:
#    param['n']=len(bg)+1

#if args.backgrounds: stackList = [ bg + [s] for s in signal ]
#else: stackList = [ [s] for s in signal ]

if args.backgrounds:
    bg += [ WZSample, ttSample ]
#    params += [ {'legendText':s.name.split('_')[1], 'WC':{}, 'color':ROOT.kRed, 'n':len(bg[i:])} for i, s in enumerate(bg) ]
    params += [ {'legendText':s.name.split('_')[1], 'WC':{}, 'color':ROOT.kRed} for s in bg ]

#append SM last
signal += [ signalSample ]
params.append( {'legendText':'SM', 'WC':{}, 'color':ROOT.kBlack} )
#params.append( {'legendText':'SM', 'WC':{}, 'color':ROOT.kBlack, 'n':len(bg)+1} )

#if args.backgrounds: stackList += [ bg + [signalSample] ]
#else: stackList += [ [signalSample] ]

if args.backgrounds: stack = Stack( bg )
else: stack = Stack()
stack.extend( [ [s] for s in signal ] ) 
#stack.extend( [ [signal[-1]] ] ) 

#stack = Stack( *stackList )

# reweighting of pTZ
if args.reweightPtXToSM:

    if   'ttZ' in args.processFile:     varX = "%sZ_pt"%args.level
    elif 'ttW' in args.processFile:     varX = "%sW_pt"%args.level
    elif 'ttgamma' in args.processFile: varX = "%sPhoton_pt"%args.level

    for param in params[::-1]:
        param['ptX_histo'] = signalSample.get1DHistoFromDraw(varX, [10,0,500], selectionString = cutInterpreter.cutString(args.selection), weightString = w.get_weight_string(**param['WC']))
        if param['ptX_histo'].Integral()>0: param['ptX_histo'].Scale(1./param['ptX_histo'].Integral())
        param['ptX_reweight_histo'] = params[-1]['ptX_histo'].Clone()
        param['ptX_reweight_histo'].Divide(param['ptX_histo'])
        logger.info( 'Made reweighting histogram for ptX and param-point %r with integral %f', param, param['ptX_reweight_histo'].Integral())

    def get_reweight( param, signal=True ):

        if signal:
            histo = param['ptX_reweight_histo']
            bsm_rw = w.get_weight_func( **param['WC'] )
            def reweight(event, sample):
                i_bin = histo.FindBin(getattr( event, varX ) )
                return histo.GetBinContent(i_bin)*bsm_rw( event, sample ) * event.ref_lumiweight1fb * float(args.luminosity) * float(sample.event_factor)
            return reweight

        else:
            def reweight(event, sample):
                return event.lumiweight1fb * float(args.luminosity) * float(sample.event_factor)
            return reweight

    weight = [ [ get_reweight( param, i!=0 if args.backgrounds else True ) for s in stackComponent ] for i, stackComponent in enumerate(stack) ]
#    weight = [ [ get_reweight( param, i==param['n']-1 ) for i in range(param['n']) ] for param in params ]
#    weight = [ [ get_reweight( param ) ] for param in params ]

else:
    def get_reweight( param , signal=True):

        def reweight_signal(event, sample):
            return w.get_weight_func( **param['WC'] )( event, sample ) * event.ref_lumiweight1fb * float(args.luminosity) * float(sample.event_factor)

        def reweight_bg(event, sample):
            return event.lumiweight1fb * float(args.luminosity) * float(sample.event_factor)

        return reweight_signal if signal else reweight_bg

    weight = [ [ get_reweight( param, i!=0 if args.backgrounds else True ) for s in stackComponent ] for i, stackComponent in enumerate(stack) ]
#    weight = [ [ get_reweight( param, i==param['n']-1 ) for i in range(param['n']) ] for param in params ]

def drawObjects( hasData = False ):
    print 'drawObjects check'
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'data' if hasData else "Simulation (%s)"%args.level),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( float(args.luminosity), dataMCScale ) ) if hasData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % float(args.luminosity))
    ]
    print [tex.DrawLatex(*l) for l in lines] 
    print 'drawObjects check done'
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots):
  print 'drawPlots'
  for plot in plots:
    print 'plot in plots'
    for i_h, h in enumerate(plot.histos):
      print 'histo h'
      h[0].style = styles.lineStyle(params[i_h]['color'])

  for log in [False, True]:
    print 'log ' + str(log)
    # Directory structure
    WC_directory = '_'.join(args.parameters).rstrip('0').replace('-','m').replace('.','p') if len(args.parameters)>1 else 'SM'
    plot_directory_ = os.path.join(\
        plot_directory,
        '%s_%s'%(args.level, args.version),
        signalSample.name, 
        subDirectory, 
        args.selection if args.selection is not None else 'no_selection', 
        WC_directory,
        "log" if log else "lin")

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
      print 'plot 2 in plots'
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

#sequence functions
sequence = process.getSequenceList( args.level )

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, addOverFlowBin=None)

plots = process.getPlotList( args.scaleLumi, args.level )

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

drawPlots(plots)
