#!/usr/bin/env python
''' Analysis script for standard plots '''

# Standard imports and batch mode
import ROOT, os, itertools
#from more_itertools import peekable
ROOT.gROOT.SetBatch(True)

from math import sqrt, cos, sin, pi, isnan, sinh, cosh
import copy
import imp
import pickle

# RootTools
from RootTools.core.standard import *

# TTXPheno
from TTXPheno.Tools.user import plot_directory
from TTXPheno.Tools.helpers import deltaPhi, getCollection, deltaR, nanJet, nanLepton, getObjDict
from TTXPheno.Tools.WeightInfo import WeightInfo, histo_to_list

# Import samples
from TTXPheno.samples.benchmarks import *

# Import helpers
from plot_helpers import *

# Arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',        action='store',     default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--version',         action='store',     default='v7', help='Appendix to plot directory') 
argParser.add_argument('--processFile',     action='store',     default='ttZ_3l', nargs='?', choices=['ttZ_3l', 'ttZ_4l', 'ttW_2l', 'ttgamma_1l', 'ttgamma_2l'], help='Which process? ttZ, ttW, ttgamma') 
argParser.add_argument('--sample',          action='store',     default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',           action='store',     default=2, help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',       action='store',     default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',           action='store_true', help='Run only on a small subset of the data?') 
argParser.add_argument('--backgrounds',     action='store_true', help='include backgrounds?')
argParser.add_argument('--level',           action='store',     default='gen', nargs='?', choices=['reco', 'gen', 'genLep'], help='Which level of reconstruction? reco, gen, genLep')
argParser.add_argument('--scaleLumi',       action='store_true', help='Scale lumi only?')
argParser.add_argument('--reweightPtXToSM', action='store_true', help='Reweight Pt(X) to the SM for all the signals?')
argParser.add_argument('--parameters',      action='store',     default = ['ctW', '3', 'ctWI', '3', 'ctZ', '3', 'ctZI', '3'], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',      action='store',     default=150, help='Luminosity for weighting the plots')
argParser.add_argument('--leptonFlavor',    action='store',     default='all', nargs='?', choices=['all', 'same', 'opposite', 'e', 'mu'], help='same flavor of nonZ leptons for ttZ 4l and ttgamma 2l? No effect on other processes') 
argParser.add_argument('--variables',       action='store',     default = [], type=str, nargs='+', help = "argument variables")
argParser.add_argument('--binThreshold',    action='store',     default=100)
argParser.add_argument('--addFisherInformation', action='store_true', help='include Fisher Information Plot in a.u.?')

args = argParser.parse_args()

if len(args.parameters) < 2: args.parameters = None

#Fisher Information Wilson Coeff, {} for SM
fisherInfo_WC = {}

# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco':
    from TTXPheno.Tools.cutInterpreterReco import cutInterpreter
    from TTXPheno.Tools.objectSelection import isGoodRecoJet as isGoodJet
    from TTXPheno.Tools.objectSelection import isGoodRecoLepton as isGoodLepton

elif args.level == 'genLep':
    from TTXPheno.Tools.cutInterpreterGenLep import cutInterpreter
    from TTXPheno.Tools.objectSelection import isGoodGenJet as isGoodJet
    from TTXPheno.Tools.objectSelection import isGoodGenLepton as isGoodLepton

#remove Gen here for getting the reas cutInterpreter (this is just to check with old plots)
else:
    from TTXPheno.Tools.cutInterpreterGen import cutInterpreter
    from TTXPheno.Tools.objectSelection import isGoodGenJet as isGoodJet
    from TTXPheno.Tools.objectSelection import isGoodGenLepton as isGoodLepton

preTag = 'reco' if args.level == 'reco' else 'gen'
tag = 'reco' if args.level == 'reco' else 'genLep'

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger = logger.get_logger( args.logLevel, logFile = None )
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None )

# Make subdirectory
subDirectory = []
if args.scaleLumi: subDirectory.append("shape")
else: subDirectory.append("lumi")
if args.reweightPtXToSM: subDirectory.append("reweightPtXToSM")
if args.small: subDirectory.append("small")
subDirectory = '_'.join( subDirectory )

# Format WC input parameters
colors = [ ROOT.kOrange, ROOT.kOrange+10, ROOT.kCyan+1, ROOT.kGreen+1, ROOT.kBlue,ROOT.kBlue-2, ROOT.kRed, ROOT.kRed+2, ROOT.kViolet+2, ROOT.kYellow+2, ROOT.kRed-7, ROOT.kPink-7, ROOT.kPink-3, ROOT.kGreen+4, ROOT.kGray+2 ]

params = []
if args.parameters is not None:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals = list( map( float, str_vals ) )
    for i_param, (coeff, val, str_val) in enumerate(zip(coeffs, vals, str_vals)):
        params.append( [{
            'legendText': ' '.join([coeff,str_val]).replace('gamma','#gamma'),
            'WC' : { coeff:val },
            'color' : colors[i_param],
            }])

params.append( [{'legendText':'SM', 'WC':{}, 'color':ROOT.kBlack}] )

# Import process specific variables
process_file = os.path.join( os.path.dirname( os.path.realpath( __file__ ) ), 'addons', '%s.py'%args.processFile )
process = imp.load_source( "process", os.path.expandvars( process_file ) )

#root file variables
read_variables = process.getVariableList( args.level )
#sequence functions
sequence = process.getSequenceList( args.level, args.leptonFlavor )

# Import samples
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
loadedSamples = imp.load_source( "samples", os.path.expandvars( sample_file ) )

ttXSample = getattr( loadedSamples, args.sample )
WZSample = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights' )
ttSample = getattr( loadedSamples, 'fwlite_tt_lep_LO_order2_15weights' )
ttSemiLepSample = getattr( loadedSamples, 'fwlite_tt_semilep_LO_order2_15weights' )
tWSample = getattr( loadedSamples, 'fwlite_tW_LO_order2_15weights' )
tWZSample = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights' )
tZqSample = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights' )
ZgammaSample = getattr( loadedSamples, 'fwlite_Zgamma_LO_order2_15weights' )
ttgammaSample = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights' )

if args.small:
    ttXSample.reduceFiles( to = 1 )
    WZSample.reduceFiles( to = 1 )
    ttSample.reduceFiles( to = 1 )
    ttSemiLepSample.reduceFiles( to = 1 )
    tWSample.reduceFiles( to = 1 )
    tWZSample.reduceFiles( to = 1 )
    tZqSample.reduceFiles( to = 1 )
    ZgammaSample.reduceFiles( to = 1 )
    ttgammaSample.reduceFiles( to = 1 )

# Polynomial parametrization
# ATTENTION IF U USE MORE THAN ONE SIGNAL SAMPLE!!!
w = WeightInfo(ttXSample.reweight_pkl)
w.set_order(int(args.order))
if len(args.variables) == 0: args.variables = w.variables

def checkReferencePoint( sample ):
    ''' check if sample is simulated with a reference point
    '''
    return pickle.load(file(sample.reweight_pkl))['ref_point'] != {}

# configure samples
for s in [ ttXSample, WZSample, ttSample, ttSemiLepSample, tWSample, tWZSample, tZqSample, ZgammaSample, ttgammaSample ]:
    # Scale the plots with number of events used (implemented in ref_lumiweight1fb)
    s.event_factor = s.nEvents / float( s.chain.GetEntries() )
    s.setSelectionString( cutInterpreter.cutString(args.selection) )
    if checkReferencePoint( s ):
        s.read_variables = ["ref_lumiweight1fb/F"]
        s.read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

signal = ttXSample
if args.processFile == 'ttZ_3l': bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]
elif args.processFile == 'ttZ_4l': bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]
#elif args.processFile == 'ttgamma_1l': bg = [ ttSample, ttSemiLepSample, tWSample, tWZSample, tZqSample, ZgammaSample ]
elif args.processFile == 'ttgamma_1l': bg = [ ttSample, tWSample, tWZSample, tZqSample, ZgammaSample ]
#elif args.processFile == 'ttgamma_2l': bg = [ ttSample, ttSemiLepSample, tWSample, tWZSample, tZqSample, ZgammaSample ]
elif args.processFile == 'ttgamma_2l': bg = [ ttSample, tWSample, tWZSample, tZqSample, ZgammaSample ]
else: bg = [ WZSample, ttSample, ttSemiLepSample, tWSample, tWZSample, tZqSample, ZgammaSample, ttgammaSample ]

#bg = [ WZSample, ttSample ]
#bg = [ ttSample ]
#bg = [ WZSample ]

if args.addFisherInformation:
    params.append( [ {'legendText':'FI SM ideal [a.u.]', 'WC':fisherInfo_WC, 'color':colors[-1]} ] )
    if args.backgrounds: params.append( [ {'legendText':'FI SM real', 'WC':fisherInfo_WC, 'color':colors[-1]} ] )

stackList = [ [signal] for param in params ]
if args.backgrounds: stackList += [ bg ]
stack = Stack( *stackList )


def legendtext( sample ):
    splitname = sample.name.split('_')
    if splitname[1] != 'tt': return splitname[1].replace('gamma','#gamma')
    else: return ' '.join(splitname[1:3])

if args.backgrounds: params += [[ {'legendText':legendtext(s), 'WC':{}, 'color':colors[-1-i]} for i, s in enumerate(bg) ]]

# reweighting of pTZ
if args.reweightPtXToSM:

    if 'ttZ' in args.processFile: varX = "%sZ_pt"%args.level
    elif 'ttW' in args.processFile: varX = "%sW_pt"%args.level
    elif 'ttgamma' in args.processFile: varX = "%sPhoton_pt"%args.level

    rwIndex = -2 if args.backgrounds else -1

    for i, param in enumerate( params[::-1] ):
        if i==0 and args.backgrounds: continue # no bg scaling
        param[0]['ptX_histo'] = ttXSample.get1DHistoFromDraw(varX, [10,0,500], selectionString = cutInterpreter.cutString(args.selection), weightString = w.get_weight_string(**param[0]['WC']))
        ptX_integral = param[0]['ptX_histo'].Integral()
        if ptX_integral > 0: param[0]['ptX_histo'].Scale(1./ptX_integral)
        param[0]['ptX_reweight_histo'] = params[rwIndex][0]['ptX_histo'].Clone()
        param[0]['ptX_reweight_histo'].Divide(param[0]['ptX_histo'])
        logger.info( 'Made reweighting histogram for ptX and param-point %r with integral %f', param[0], param[0]['ptX_reweight_histo'].Integral())

    def get_reweight( param, sample_, isSignal=True ):

        if isSignal:
            histo = param['ptX_reweight_histo']
            bsm_rw = w.get_weight_func( **param['WC'] )
            def reweight(event, sample):
                i_bin = histo.FindBin(getattr( event, varX ) )
                return histo.GetBinContent(i_bin)*bsm_rw( event, sample ) * event.ref_lumiweight1fb * float(args.luminosity) * float(sample.event_factor)

            return reweight

        else:
            def reweightRef(event, sample):
                return w.get_weight_func( **param['WC'] )( event, sample ) * event.ref_lumiweight1fb * float(args.luminosity) * float(sample.event_factor)

            def reweightNoRef(event, sample):
                return event.lumiweight1fb * float(args.luminosity) * float(sample.event_factor)

            return reweightRef if checkReferencePoint( sample_ ) else reweightNoRef

    weight = [ [ get_reweight( params[i][j], sample_, i != len(stack)-1 if args.backgrounds else True ) for j, sample_ in enumerate(stackComponent) ] for i, stackComponent in enumerate(stack) ]

else:
    def get_reweight( param , sample_ ):

        def reweightRef(event, sample):
            return w.get_weight_func( **param['WC'] )( event, sample ) * event.ref_lumiweight1fb * float(args.luminosity) * float(sample.event_factor)

        def reweightNoRef(event, sample):
            return event.lumiweight1fb * float(args.luminosity) * float(sample.event_factor)

        return reweightRef if checkReferencePoint( sample_ ) else reweightNoRef

    weight = [ [ get_reweight( params[i][j], sample_ ) for j, sample_ in enumerate(stackComponent) ] for i, stackComponent in enumerate(stack) ]

ttXWeightString = 'ref_lumiweight1fb*%s*%s*%s'%(str(args.luminosity), str(ttXSample.event_factor), w.get_weight_string( **fisherInfo_WC ))

def drawObjects( hasData = False ):
    titleAddon = ', FI: %s'%' '.join(args.variables) if args.addFisherInformation else ''
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'data' if hasData else "Simulation (%s)"%args.level),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( float(args.luminosity), dataMCScale ) ) if hasData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)%s' % (float(args.luminosity), titleAddon) )
    ]
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots):

  if args.backgrounds:
    ind = -3 if args.addFisherInformation else -1
    for plot in plots:
      for s, signal_histo in enumerate(plot.histos[:ind]):
        for bg_histo in plot.histos[-1]:
          signal_histo[0].Add(bg_histo)

  for plot in plots:
    histoIndexSM = len(plot.histos)-1
    if args.backgrounds and not args.addFisherInformation: histoIndexSM -= 1
    elif not args.backgrounds and args.addFisherInformation: histoIndexSM -= 1
    elif args.backgrounds and args.addFisherInformation: histoIndexSM -= 3

    for i_h, h in enumerate(plot.histos):
      for j_hi, hi in enumerate(h):
        if i_h == len(plot.histos)-1 and args.backgrounds:
            hi.style = styles.fillStyle(params[i_h][j_hi]['color'])
        else:
            hi.style = styles.lineStyle(params[i_h][j_hi]['color'])
            if i_h == histoIndexSM or (args.addFisherInformation and i_h == histoIndexSM+1): hi.SetLineWidth(2)
            if args.addFisherInformation and args.backgrounds and i_h == histoIndexSM+2:
                hi.SetLineStyle(2) #Fisher Info Bg
                hi.SetLineWidth(2)

  for log in [False, True]:
    # Directory structure
    WC_directory = '_'.join(args.parameters).rstrip('0').replace('-','m').replace('.','p') if args.parameters is not None else 'SM'
    plot_directory_ = os.path.join(\
        plot_directory,
        '%s_%s'%(args.level, args.version),
        ttXSample.name,
        'fisher_information' if args.addFisherInformation else '',
        'kinematics' if args.addFisherInformation else '',
        'backgrounds' if args.backgrounds else 'signalOnly',
        subDirectory,
        args.selection if args.selection is not None else 'no_selection',
        '%sLeptonFlavor'%args.leptonFlavor,
        WC_directory,
        '_'.join(args.variables) if args.addFisherInformation else '',
        '%sEventsPerBin'%str(args.binThreshold) if args.addFisherInformation else '',
        "log" if log else "lin")

    # plot the legend
    l_plot = copy.deepcopy(plots[0])
    for i_h, h in enumerate(l_plot.histos):
      for j_hi, hi in enumerate(h):
          hi.legendText = params[i_h][j_hi]['legendText']
          if i_h == len(plot.histos)-1 and args.backgrounds: hi.style = styles.fillStyle(params[i_h][j_hi]['color'])
          else: hi.style = styles.lineStyle(params[i_h][j_hi]['color'])
          hi.Scale(0.)
          hi.GetXaxis().SetTickLength(0.)
          hi.GetYaxis().SetTickLength(0.)
          hi.GetXaxis().SetLabelOffset(999.)
          hi.GetYaxis().SetLabelOffset(999.)
    l_plot.name = "legend"
    l_plot.texX = ''
    l_plot.texY = ''

    plotting.draw(l_plot,
        plot_directory = plot_directory_,
        ratio = None,
        logX = False, logY = log, sorting = True,
        legend = ( (0.17,0.9-0.05*sum(map(len, l_plot.histos))/3,0.9,0.9), 3),
        drawObjects = drawObjects(),
        copyIndexPHP = True,
    )

    # plot the plots
    for p, plot in enumerate(plots):
      histoIndexSM = len(plot.histos)-1
      if args.backgrounds and not args.addFisherInformation: histoIndexSM -= 1
      if args.addFisherInformation:
        if args.backgrounds: histoIndexSM -= 3
        else: histoIndexSM -= 1
        fi_Integral = plot.histos[histoIndexSM+1][0].Integral()
        fisherInfoScale = plot.histos[histoIndexSM][0].Integral() / fi_Integral if fi_Integral != 0 else 0 #scaling factor from SM histo
        fisherInfoScale *= 100 if log else 1.5 #offset

      for i_h, h in enumerate(plot.histos):
        for j_hi, hi in enumerate(h):
          if args.addFisherInformation and (i_h == histoIndexSM+1 or (i_h == histoIndexSM+2 and args.backgrounds)):
            if fisherInfoVariables[p] is not None: hi.Scale(fisherInfoScale) #signal
            else: hi.Scale(0) #no fisherInfo histo, only place holder
          hi.legendText = params[i_h][j_hi]['legendText']

      if not max( max(li.GetMaximum() for li in l) for l in plot.histos): continue # Empty plot

      scaleIndex = len(params)-1
      if args.backgrounds and not args.addFisherInformation: scaleIndex -= 1
      elif not args.backgrounds and args.addFisherInformation: scaleIndex -= 1
      elif args.backgrounds and args.addFisherInformation: scaleIndex -= 3

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = None,
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0., "auto"),
#        scaling = {i:(len(params)-1) for i in range(len(params)-1)} if args.scaleLumi else {}, #Scale BSM shapes to SM (last in list)
        scaling = {i:(scaleIndex) for i in range(scaleIndex)} if args.scaleLumi else {}, #Scale BSM shapes to SM (last in list)
	    legend = ( (0.17,0.9-0.05*sum(map(len, plot.histos))/3,0.9,0.9), 3),
	    drawObjects = drawObjects(),
        copyIndexPHP = True,
      )


#logger.info( "Translating cut %s to %s", args.selection, cutInterpreter.cutString(args.selection) )

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, addOverFlowBin=None)

plots, fisherInfoVariables = process.getPlotList( args.scaleLumi, args.level )

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

if args.addFisherInformation:
    for i, plot in enumerate(plots):
        if fisherInfoVariables[i] is None: continue

        fiHistoListIndex = len(plot.histos)-3 if args.backgrounds else len(plot.histos)-1
        bins = plot.binning
        var = fisherInfoVariables[i] #replace with better solution than FI_plots!!

        hist = w.getFisherInformationHisto( ttXSample, var, bins, selectionString=cutInterpreter.cutString(args.selection), weightString=ttXWeightString, variables=args.variables, nEventsThresh = args.binThreshold, **fisherInfo_WC)
        plot.histos[fiHistoListIndex] = [plot.fromHisto( 'histo_%i'%i, hist ).histos]

        if args.backgrounds:
            #get list of contents of bg histo
            bgContentList = histo_to_list( plot.histos_added[-1][0] )
            #get list of contents of signal histo
            sigContentList = histo_to_list( plot.histos[fiHistoListIndex-1][0] )
            #calculate purity for each bin
            purityBinList = [ sig / ( sig + bgContentList[j] ) if (sig + bgContentList[j])!=0 else 0 for j, sig in enumerate( sigContentList ) ]
            # multiply ideal fisher info histo purity
            for l in range(bins[0]):
                plot.histos[fiHistoListIndex+1][0].SetBinContent(l+1, plot.histos[fiHistoListIndex][0].GetBinContent(l+1) * purityBinList[l])

drawPlots(plots)

