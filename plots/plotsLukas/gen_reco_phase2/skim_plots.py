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
argParser.add_argument('--version',         action='store',     default='test', help='Appendix to plot directory') 
argParser.add_argument('--processFile',     action='store',     default='ttZ_3l', nargs='?', choices=['ttZ_3l', 'ttZ_4l', 'ttW_2l', 'ttgamma_1l', 'ttgamma_2l', 'ttZ_3l_small', 'ttgamma_1l_small', 'ttgamma_2l_small', 'ttZ_3l_paper', 'ttgamma_1l_paper', 'ttgamma_2l_paper'], help='Which process? ttZ, ttW, ttgamma') 
argParser.add_argument('--sample',          action='store',     default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',           action='store',     default=2, help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',       action='store',     default='lepSel3-njet3p-nbjet1p-onZ-Zpt0-mll12', help="Specify cut.")
argParser.add_argument('--small',           action='store_true', help='Run only on a small subset of the data?') 
argParser.add_argument('--backgrounds',     action='store_true', help='include backgrounds?')
argParser.add_argument('--noninfoSignal',   action='store_true', help='include non-info signal?')
argParser.add_argument('--scale14TeV',      action='store_true', help='scale 13 TeV cross-sections to 14 Tev?')
#argParser.add_argument('--level',           action='store',     default='gen', nargs='?', choices=['reco', 'gen'], help='Which level of reconstruction? reco, gen')
argParser.add_argument('--level',           action='store',     default='reco', nargs='?', choices=['reco', 'gen'], help='Which level of reconstruction? reco, gen')
argParser.add_argument('--scaleLumi',       action='store_true', help='Scale lumi only?')
argParser.add_argument('--reweightPtXToSM', action='store_true', help='Reweight Pt(X) to the SM for all the signals?')
argParser.add_argument('--parameters',      action='store',     default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--luminosity',      action='store',     default=150, help='Luminosity for weighting the plots')
argParser.add_argument('--leptonFlavor',    action='store',     default='all', nargs='?', choices=['all', 'same', 'opposite', 'e', 'mu', 'eee', 'mumumu', 'mumue', 'muee'], help='same flavor of nonZ leptons for ttZ 4l and ttgamma 2l? No effect on other processes') 
argParser.add_argument('--variables',       action='store',     default = ['cpt'], type=str, nargs='+', help = "argument variables")
argParser.add_argument('--binThreshold',    action='store',     default=100)
argParser.add_argument('--addFisherInformation', action='store_true', help='include Fisher Information Plot in a.u.?')
argParser.add_argument('--addFisherInformationBackground', action='store_true', help='include Fisher Information bg Plot in a.u.?')
argParser.add_argument('--detector',        action='store',     default='CMS', nargs='?', choices=['CMS', 'ATLAS', 'phase2_CMS'], help='Which Delphes detector simulation?') 
argParser.add_argument('--wideLegendStyle', action='store_true', help='Legend in 4 columns') 

args = argParser.parse_args()

if len(args.parameters) < 2: args.parameters = None

#Fisher Information Wilson Coeff, {} for SM
fisherInfo_WC = {}

# Import additional functions/classes specified for the level of reconstruction
if args.level == 'reco':
    from TTXPheno.Tools.cutInterpreterReco import cutInterpreter
    from TTXPheno.Tools.objectSelection import isGoodRecoJet as isGoodJet
    from TTXPheno.Tools.objectSelection import isGoodRecoLepton as isGoodLepton
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
colors = [ ROOT.kRed+1, ROOT.kGreen+2, ROOT.kOrange+1, ROOT.kViolet+9, ROOT.kSpring-7, ROOT.kRed+2,  ROOT.kPink-9, ROOT.kBlue,  ROOT.kRed-7, ROOT.kRed-10, ROOT.kRed+3,  ROOT.kGreen-7, ROOT.kGreen-10, ROOT.kGreen+3,  ROOT.kPink-7, ROOT.kPink-10, ROOT.kPink+3, ROOT.kGray+2, ROOT.kYellow-7 ]
colorsBg = [ ROOT.kAzure-3, ROOT.kGreen-2, ROOT.kCyan-9, ROOT.kRed+2, ROOT.kGray+2, ROOT.kYellow-7, ROOT.kViolet+6, ROOT.kBlue+2 ]
colorsNonInfo = [ ROOT.kRed-7, ROOT.kRed-10, ROOT.kRed+3, ROOT.kBlack]

params = []
if args.parameters is not None:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals = list( map( float, str_vals ) )
    for i_param, (coeff, val, str_val) in enumerate(zip(coeffs, vals, str_vals)):
        params.append( [{
#            'legendText': ' '.join([coeff,str_val]).replace('gamma','#gamma').replace('c','C_{').replace(' ','} = ').replace('p','#phi').replace('M','').replace('I','}^{[Im]') + '  ',
            'legendText': ' '.join([coeff]).replace('gamma','#gamma').replace('c','C_{').replace('p','#phi').replace('M','').replace('I','}^{[Im]') + '}  ',
            'WC' : { coeff:val },
            'color' : colors[i_param],
            }])

params.append( [{'legendText':args.processFile.split('_')[0].replace('gamma','#gamma'), 'WC':{}, 'color':ROOT.kBlack}] )

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

if args.processFile == 'ttgamma_1l': ttSampleName = 'fwlite_tt_nonhad_LO_order2_15weights'
else:                                ttSampleName = 'fwlite_tt_dilep_LO_order2_15weights'

ttXSample = getattr( loadedSamples, args.sample + '_%s' %args.detector)
WZSample = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights_%s' %args.detector )
#tWSample = getattr( loadedSamples, 'fwlite_tW_LO_order2_15weights_%s' %args.detector )
tWZSample = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights_%s' %args.detector )
tZqSample = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights_%s' %args.detector )
#ZgammaSample = getattr( loadedSamples, 'fwlite_Zgamma_LO_order2_15weights_%s' %args.detector )
ttgammaSample = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights_%s' %args.detector )
ttgammaSample.name = 'fwlite_ttgamma__LO_order2_15weights_%s' %args.detector
if "ttgamma" in args.processFile:
    ttSample = getattr( loadedSamples, ttSampeName + '_' + args.detector )
    ttSample.name = 'fwlite_tt__LO_order2_15weights_%s' %args.detector

# details of the categories are written in the postprocessing script
if args.processFile.split('_')[0] == 'ttgamma':

    ttgammaIsrSample  = copy.deepcopy( ttXSample ) #select ttgamma events with isolated gamma from ISR (cat a2)
    ttgammaIsrSample.name = 'fwlite_ttgamma_ISR_LO_order2_15weights_ref'

    ttgammaLepSample  = copy.deepcopy( ttSample ) #getattr( loadedSamples, ttSampleName ) #select tt events with isolated gamma from W, l, tau (cat b)
    ttgammaLepSample.name = 'fwlite_ttgamma_(W,l,tau)_LO_order2_15weights'

    ttgammabSample  = copy.deepcopy( ttSample ) #getattr( loadedSamples, ttSampleName ) #select tt events with isolated gamma from W, l, tau (cat b)
    ttgammabSample.name = 'fwlite_ttgamma_(b,j)_LO_order2_15weights'

    ttgammaFakeSample = copy.deepcopy( ttSample ) #getattr( loadedSamples, ttSampleName ) #select tt events with no gamma -> if in the plots: fake gamma (cat d)
    ttgammaFakeSample.name = 'fwlite_ttgamma_fake_LO_order2_15weights'

    # ttSample #select tt events with non-isolated gamma or gamma from bottom (cat c1 + c2)
    # signal: ttgamma events with isolated gamma from gluon or top (cat a1)

elif args.processFile.split('_')[0] == 'ttZ':
    ttZISRSample  = copy.deepcopy( ttXSample ) #select ttgamma events with isolated gamma from ISR (cat a2)
    ttZISRSample.name = 'fwlite_ttZ_(non-info)_LO_order2_15weights_ref'

nonInfo = []
if 'ttZ' in args.processFile.split('_'):
    bg = [ WZSample, tWZSample, tZqSample, ttgammaSample ]
#    bg = [ WZSample, tWZSample, tZqSample ]
    # be careful if you set nonInfo to empty list, especially with the FI plot
    if args.noninfoSignal or args.addFisherInformationBackground: nonInfo = [ ttZISRSample ]
elif 'ttgamma' in args.processFile.split('_'):
    bg = [ ttSample, tWSample, tWZSample, tZqSample, ZgammaSample ]
    # be careful if you set nonInfo to empty list, especially with the FI plot
    if args.noninfoSignal or args.addFisherInformationBackground: nonInfo = [ ttgammaIsrSample, ttgammaLepSample, ttgammabSample, ttgammaFakeSample ]

# Polynomial parametrization
# ATTENTION IF U USE MORE THAN ONE SIGNAL SAMPLE!!!
w = WeightInfo(ttXSample.reweight_pkl)
w.set_order(int(args.order))
if len(args.variables) == 0: args.variables = w.variables

def checkReferencePoint( sample ):
    ''' check if sample is simulated with a reference point
    '''
    return pickle.load(file(sample.reweight_pkl))['ref_point'] != {}

# somehow this has to be done first, not in the next loop
if args.small:
    for s in [ttXSample] + bg + nonInfo:
        s.reduceFiles( to = 30 )


# configure samples
for s in [ttXSample] + bg + nonInfo:
    # Scale the plots with number of events used (implemented in ref_lumiweight1fb)
    s.event_factor    = s.nEvents / float( s.chain.GetEntries() )
    s.xsecScaleFactor = s.xsec14 / s.xsec if args.scale14TeV else 1.
    s.setSelectionString( cutInterpreter.cutString(args.selection) )
    if checkReferencePoint( s ):
        print s.name
        s.read_variables = ["ref_lumiweight1fb/F", VectorTreeVariable.fromString('p[C/F]', nMax=2000)]

catPhoton_variables = [ "signalPhoton/I", "isrPhoton/I", "lepPhoton/I", "nonIsoPhoton/I", "fakePhoton/I"]
catZ_variables = [ "signalZ/I" ]

if args.processFile.split('_')[0] == 'ttgamma':
    # overlap removal + signal categorization for ttgamma

    ttXSample.read_variables         += catPhoton_variables
    ttXSample.addSelectionString(         "signalPhoton==1" ) #cat a1 <- the one and only signal

    ttSample.read_variables           = catPhoton_variables
    ttSample.addSelectionString(          "nonIsoPhoton==1" ) #cat c1

    if args.noninfoSignal or args.addFisherInformationBackground:
        ttgammaIsrSample.read_variables  += catPhoton_variables
        ttgammaLepSample.read_variables   = catPhoton_variables
        ttgammabSample.read_variables     = catPhoton_variables
        ttgammaFakeSample.read_variables  = catPhoton_variables

        ttgammaIsrSample.addSelectionString(  "isrPhoton==1&&abs(genPhoton_motherPdgId[0])!=21"    ) #cat a2
        ttgammaLepSample.addSelectionString(  "lepPhoton==1&&abs(genPhoton_motherPdgId[0])!=6"    ) #cat b
        ttgammabSample.addSelectionString(    "jetPhoton==1" ) #cat c2
        ttgammaFakeSample.addSelectionString( "fakePhoton==1"   ) #cat d

elif args.processFile.split('_')[0] == 'ttZ':
    # signal categorization for ttZ

    ttXSample.read_variables += catZ_variables
    ttXSample.addSelectionString( "signalZ==1" ) #Z from gluon or top

    if args.noninfoSignal or args.addFisherInformationBackground:
        ttZISRSample.read_variables += catZ_variables
        ttZISRSample.addSelectionString( "signalZ==0" ) #Z from ISR or else



def legendtext( sample ):
    splitname = sample.name.split('_')
    if splitname[1] not in ['tt','ttgamma','ttZ']: return splitname[1].replace('gamma','#gamma')
    else: return ' '.join(splitname[1:3]).replace('gamma','#gamma').replace('tau','#tau')

#first: draw all WC + SM as line
stackList = [ [ttXSample] for param in params ]
#second: draw all non-info signal red
nonInfoParams = []
if args.noninfoSignal or args.addFisherInformationBackground:
    stackList += [ nonInfo ]
    nonInfoParams = [[ {'legendText':legendtext(s), 'WC':{}, 'color':colorsNonInfo[i]} for i, s in enumerate(nonInfo) ]]
#third: draw all bg filled
bgParams = []
if args.backgrounds or args.addFisherInformationBackground:
    stackList += [ bg ]
    bgParams = [[ {'legendText':legendtext(s), 'WC':{}, 'color':colorsBg[i]} for i, s in enumerate(bg) ]]

fisherParams = []
if args.addFisherInformation:
    fisherParams = [[ {'legendText':'FI SM ideal [a.u.]', 'WC':fisherInfo_WC, 'color':ROOT.kGray+2} ]]
    if args.addFisherInformationBackground: fisherParams.append( [{'legendText':'FI SM real', 'WC':fisherInfo_WC, 'color':ROOT.kGray+2}] )
    #forth: draw FI as line
    stackList += [ [ttXSample] for param in fisherParams ]
#    allParams += fisherParams

stack = Stack( *stackList )

# reweighting of pTZ
if args.reweightPtXToSM:

    if 'ttZ' in args.processFile: varX = "%sZ_pt"%args.level
    elif 'ttW' in args.processFile: varX = "%sW_pt"%args.level
    elif 'ttgamma' in args.processFile: varX = "%sPhoton_pt"%args.level

#    rwIndex = -2 if args.backgrounds else -1

    for i, param in enumerate( params[::-1] ):
#        if i==0 and args.backgrounds: continue # no bg scaling
        param[0]['ptX_histo'] = ttXSample.get1DHistoFromDraw(varX, [50,0,500], selectionString = cutInterpreter.cutString(args.selection), weightString = w.get_weight_string(**param[0]['WC']))
        ptX_integral = param[0]['ptX_histo'].Integral()
        if ptX_integral > 0: param[0]['ptX_histo'].Scale(1./ptX_integral)
        param[0]['ptX_reweight_histo'] = params[-1][0]['ptX_histo'].Clone()
        param[0]['ptX_reweight_histo'].Divide(param[0]['ptX_histo'])
        logger.info( 'Made reweighting histogram for ptX and param-point %r with integral %f', param[0], param[0]['ptX_reweight_histo'].Integral())


allParams = params + nonInfoParams + bgParams + fisherParams

# reweighting of pTZ
if args.reweightPtXToSM:

    def get_reweight( param, sample_, isSignal=True ):

        if isSignal:
            histo = param['ptX_reweight_histo']
            bsm_rw = w.get_weight_func( **param['WC'] )
            def reweight(event, sample):
                i_bin = histo.FindBin(getattr( event, varX ) )
                return histo.GetBinContent(i_bin)*bsm_rw( event, sample ) * event.ref_lumiweight1fb * float(args.luminosity) * float(sample.event_factor) * float(sample.xsecScaleFactor)
#                return histo.GetBinContent(i_bin)*bsm_rw( event, sample ) * sample_.xsec * 1000 / sample_.nEvents / event.p_C[0] * float(args.luminosity) * float(sample.event_factor)

            return reweight

        else:
            def reweightRef(event, sample):
                return w.get_weight_func( **param['WC'] )( event, sample ) * event.ref_lumiweight1fb * float(args.luminosity) * float(sample.event_factor) * float(sample.xsecScaleFactor)
#                return w.get_weight_func( **param['WC'] )( event, sample ) * sample_.xsec * 1000 / sample_.nEvents / event.p_C[0] * float(args.luminosity) * float(sample.event_factor)

            def reweightNoRef(event, sample):
                return event.lumiweight1fb * float(args.luminosity) * float(sample.event_factor) * float(sample.xsecScaleFactor)
#                return sample_.xsec * 1000 / sample_.nEvents * float(args.luminosity) * float(sample.event_factor)

            return reweightRef if checkReferencePoint( sample_ ) else reweightNoRef

    weight = [ [ get_reweight( allParams[i][j], sample_, i < len(params) ) for j, sample_ in enumerate(stackComponent) ] for i, stackComponent in enumerate(stack) ]

else:
    def get_reweight( param , sample_ ):

        def reweightRef(event, sample):
            return w.get_weight_func( **param['WC'] )( event, sample ) * event.ref_lumiweight1fb * float(args.luminosity) * float(sample.event_factor) * float(sample.xsecScaleFactor)
#            return w.get_weight_func( **param['WC'] )( event, sample ) * sample_.xsec * 1000 / sample_.nEvents / event.p_C[0] * float(args.luminosity) * float(sample.event_factor)

        def reweightNoRef(event, sample):
            return event.lumiweight1fb * float(args.luminosity) * float(sample.event_factor) * float(sample.xsecScaleFactor)
#            return sample_.xsec * 1000 / sample_.nEvents * float(args.luminosity) * float(sample.event_factor)

        return reweightRef if checkReferencePoint( sample_ ) else reweightNoRef

    weight = [ [ get_reweight( allParams[i][j], sample_ ) for j, sample_ in enumerate(stackComponent) ] for i, stackComponent in enumerate(stack) ]

ttXWeightString = 'ref_lumiweight1fb*%s*%s*%s*%s'%(str(args.luminosity), str(ttXSample.event_factor), w.get_weight_string( **fisherInfo_WC ), str(ttXSample.xsecScaleFactor))

def drawObjects( hasData = False ):
    offset = 0.5 if args.addFisherInformation else 0.6
    titleAddon = ', FI: %s'%' '.join(args.variables) if args.addFisherInformation else ''
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    tex.SetTextFont(42)
    lines = [
      (0.15, 0.95, '#bf{CMS} #it{Simulation Preliminary}'),
      (offset, 0.95, '%i fb{}^{-1} (%s TeV)'% ( int(args.luminosity), '14' if args.scale14TeV else '13' ) )
#      (0.15, 0.95, ' '.join(args.processFile.split('_')[:2]) + '(' + args.detector + ')'),
#      (offset, 0.95, '%3.1f fb{}^{-1} @ 13 TeV%s'% ( float(args.luminosity), titleAddon) )
    ]
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots):

  # add nonInfo Signal to signal for stacked plot
  if args.noninfoSignal and len(nonInfoParams) != 0:
    indexNonInfo = len(params)
    for plot in plots:
      for nonInfo_histo in plot.histos[indexNonInfo]:
        for s, signal_histos in enumerate(plot.histos[:indexNonInfo]):
          signal_histos[0].Add(nonInfo_histo)

  # add bg to signal and nonInfo Signal for stacked plot
  if args.backgrounds and len(bgParams) != 0:
    indexBg = len(params) + len(nonInfoParams)
    for plot in plots:
      for bg_histo in plot.histos[indexBg]:
        for s, signal_histos in enumerate(plot.histos[:indexBg]):
            signal_histos[0].Add(bg_histo)

  for plot in plots:

    histoIndexSM      = len(params) - 1
    histoIndexNonInfo = len(params) if len(nonInfoParams) != 0 or args.addFisherInformationBackground else float('inf')
    histoIndexBg      = len(params) + len(nonInfoParams) if len(bgParams) != 0 or args.addFisherInformationBackground else float('inf')

    for i_h, h in enumerate(plot.histos):

#      if args.addFisherInformation and args.addFisherInformationBackground:
#        if not args.backgrounds and i_h == histoIndexBg: continue
#        if not args.noninfoSignal and i_h == histoIndexNonInfo: continue

      for j_hi, hi in enumerate(h):

        if i_h == histoIndexBg:
            # fill style for bg
            hi.style = styles.fillStyle(allParams[i_h][j_hi]['color'])
        elif i_h == histoIndexNonInfo:
            # fill style for nonInfo Signal
            hi.style = styles.fillStyle(allParams[i_h][j_hi]['color'])
            hi.SetFillStyle(3544)
            hi.SetFillColor(ROOT.kWhite)
            ROOT.gStyle.SetHatchesLineWidth(2)
        else:
            # fill style for signal and WC
            hi.style = styles.lineStyle(allParams[i_h][j_hi]['color'])
            hi.SetLineWidth(2)
            if i_h == histoIndexSM:
                hi.SetLineWidth(3)
            elif args.addFisherInformation and i_h == len(allParams)-2:
                hi.SetLineWidth(2)
            elif args.addFisherInformation and args.addFisherInformationBackground and i_h == len(allParams)-1:
                hi.SetLineStyle(2) #Fisher Info Bg
                hi.SetLineWidth(2)
        
  for log in [False, True]:
    # Directory structure
    WC_directory = '_'.join(args.parameters).rstrip('0').replace('-','m').replace('.','p') if args.parameters is not None else 'SM'
    plot_directory_ = os.path.join(\
        plot_directory,
        '%s_%s'%(args.level, args.version),
        args.detector,
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
        "log_%s"%('14TeV' if args.scale14TeV else '13TeV') if log else "lin_%s"%('14TeV' if args.scale14TeV else '13TeV'))

    # plot the legend
    l_plot = copy.deepcopy(plots[0])
    for i_h, h in enumerate(l_plot.histos):

#      if args.addFisherInformation and args.addFisherInformationBackground:
#        if not args.backgrounds and i_h == histoIndexBg: continue
#        if not args.noninfoSignal and i_h == histoIndexNonInfo: continue

      for j_hi, hi in enumerate(h):

        if args.addFisherInformation and args.addFisherInformationBackground:
          if not args.backgrounds and j_hi == histoIndexBg: continue
          if not args.noninfoSignal and j_hi == histoIndexNonInfo: continue

          hi.legendText = allParams[i_h][j_hi]['legendText']
          if i_h == histoIndexNonInfo or i_h == histoIndexBg: hi.style = styles.fillStyle(allParams[i_h][j_hi]['color'])
          else: hi.style = styles.lineStyle(allParams[i_h][j_hi]['color'])
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
        legend = ( (0.55,0.88-0.05*sum(map(len, l_plot.histos))/2,0.9,0.88), 2) if not args.wideLegendStyle else ( (0.22,0.88-0.05*sum(map(len, l_plot.histos))/4,0.9,0.88), 4),
        drawObjects = drawObjects(),
        copyIndexPHP = True,
    )

    # plot the plots
    for p, plot in enumerate(plots):
      histoIndexSM      = len(params) - 1
      histoIndexBg      = len(params) + len(nonInfoParams) if len(bgParams) != 0 or args.addFisherInformationBackground else float('inf')
      histoIndexNonInfo = len(params) if len(nonInfoParams) != 0 or args.addFisherInformationBackground else float('inf')
      histoIndexFI      = len(plot.histos) - 2 if args.addFisherInformationBackground else len(plot.histos) - 1

      if args.addFisherInformation:
        fi_Integral = plot.histos[histoIndexFI][0].Integral()
        fisherInfoScale = plot.histos[histoIndexSM][0].Integral() / fi_Integral if fi_Integral != 0 else 0 #scaling factor from SM histo
        fisherInfoScale *= 100 if log else 1.8 #offset
        if not log and not args.scaleLumi: fisherInfoScale *= 1.5 #offset

      for i_h, h in enumerate(plot.histos):

#        if args.addFisherInformation and args.addFisherInformationBackground:
#          if not args.backgrounds and i_h == histoIndexBg: continue
#          if not args.noninfoSignal and i_h == histoIndexNonInfo: continue

        for j_hi, hi in enumerate(h):
#          if args.addFisherInformation and args.addFisherInformationBackground:
#            if not args.backgrounds and j_hi == histoIndexBg: continue
#            if not args.noninfoSignal and j_hi == histoIndexNonInfo: continue

          if args.addFisherInformation and (i_h == histoIndexFI or (i_h == histoIndexFI+1 and args.addFisherInformationBackground)):
            if fisherInfoVariables[p] is not None: hi.Scale(fisherInfoScale) #signal
            else: hi.Scale(0) #no fisherInfo histo, only place holder
          hi.legendText = allParams[i_h][j_hi]['legendText']
          hi.GetXaxis().SetTickLength(0.04)
          hi.GetYaxis().SetTickLength(0.04)
#          hi.GetXaxis().SetLabelOffset(1.4)
#          hi.GetYaxis().SetLabelOffset(1.2)
#          hi.GetXaxis().SetTitleOffset(1.4)
          hi.GetYaxis().SetTitleOffset(1.)
          hi.GetXaxis().SetTitleSize(0.035)
          hi.GetYaxis().SetTitleSize(0.035)

#          ROOT.gStyle.SetLegendTextSize(2)

      if not max( max(li.GetMaximum() for li in l) for l in plot.histos): continue # Empty plot

      for i_h, h in enumerate(plot.histos[::-1]):
        if args.addFisherInformation and args.addFisherInformationBackground:
          if not args.backgrounds and i_h == len(plot.histos)-histoIndexBg-1:
            del plot.histos[i_h]
          if not args.noninfoSignal and i_h == len(plot.histos)-histoIndexNonInfo-1:
            del plot.histos[i_h]

#      for i_h, h in enumerate(plot.stack[::-1]):
#        if args.addFisherInformation and args.addFisherInformationBackground:
#          if not args.backgrounds and i_h == len(plot.histos)-histoIndexBg:
#            del plot.stack[i_h]
#          if not args.noninfoSignal and i_h == len(plot.histos)-histoIndexNonInfo:
#            del plot.stack[i_h]

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = None,
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0., "auto"),
#        scaling = {i:(len(params)-1) for i in range(len(params)-1)} if args.scaleLumi else {}, #Scale BSM shapes to SM (last in list)
        scaling = {i:(histoIndexSM) for i in range(histoIndexSM)} if args.scaleLumi else {}, #Scale BSM shapes to SM (last in list)
        legend = ( (0.55,0.88-0.05*sum(map(len, l_plot.histos))/2,0.9,0.88), 2)  if not args.wideLegendStyle else ( (0.22,0.88-0.05*sum(map(len, l_plot.histos))/4,0.9,0.88), 4),
#	    legend = ( (0.17,0.9-0.05*sum(map(len, plot.histos))/3,0.9,0.9), 3),
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

        fiHistoListIndex = len(plot.histos)-2 if args.addFisherInformationBackground else len(plot.histos)-1
        bins = plot.binning
        var = fisherInfoVariables[i] #replace with better solution than FI_plots!!

        hist = w.getFisherInformationHisto( ttXSample, var, bins, selectionString=cutInterpreter.cutString(args.selection), weightString=ttXWeightString, variables=args.variables, nEventsThresh = args.binThreshold, **fisherInfo_WC)
        plot.histos[fiHistoListIndex] = [plot.fromHisto( 'histo_%i'%i, hist ).histos]

        if args.addFisherInformationBackground:
            #get list of contents of bg histo
            bgContentList = histo_to_list( plot.histos_added[len(params)+1][0] )
            #get list of contents of bg histo
            nonInfoContentList = histo_to_list( plot.histos_added[len(params)][0] )
            #get list of contents of signal histo
            sigContentList = histo_to_list( plot.histos[len(params)-1][0] )
            #calculate purity for each bin
            purityBinList = [ ( sig + nonInfoContentList[j] ) / ( sig + nonInfoContentList[j] + bgContentList[j] ) if (sig + nonInfoContentList[j] + bgContentList[j])!=0 else 0 for j, sig in enumerate( sigContentList ) ]
            # multiply ideal fisher info histo purity
            for l in range(bins[0]):
                plot.histos[fiHistoListIndex+1][0].SetBinContent(l+1, plot.histos[fiHistoListIndex][0].GetBinContent(l+1) * purityBinList[l])

drawPlots(plots)

