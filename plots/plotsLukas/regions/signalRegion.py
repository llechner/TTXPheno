''' Plot script WC parameter LogLikelihood
'''

# Standard imports 
import sys
import ROOT
import imp
import pickle
import ctypes
import numpy as np
import copy

from math import sqrt
# turn off graphics
ROOT.gROOT.SetBatch( True )

# RootTools
from RootTools.core.standard import *
from plot_helpers            import getUncertaintyValue, getObservationValue
from TTXPheno.samples.color  import color

# Logger
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   'DEBUG', logFile = None)
logger_rt = logger_rt.get_logger('INFO', logFile = None)

# TTXPheno
from TTXPheno.samples.benchmarks import * 
from TTXPheno.Tools.user import plot_directory, cardfileLocation

# get the reweighting function
from TTXPheno.Tools.WeightInfo import WeightInfo

ROOT.gStyle.SetNumberContours(255)

# Arguments
import argparse

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--version',            action='store',     default='test', help='Appendix to plot directory')
argParser.add_argument('--process',            action='store',     default='ttZ', nargs='?', choices=['ttZ'], help="which process to calculate?")
argParser.add_argument('--sample',             action='store',     default='fwlite_ttZ_ll_LO_order2_15weights_ref', help='Sample name specified in sample/python/benchmarks.py, e.g. fwlite_ttZ_ll_LO_order2_15weights_ref')
argParser.add_argument('--order',              action='store',     default=2, help='Polynomial order of weight string (e.g. 2)')
argParser.add_argument('--selection',          action='store',     default='lepSel3-onZ-njet3p-nbjet1p-Zpt0', help="Specify cut.")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--level',              action='store',     default='reco', nargs='?', choices=['reco', 'gen'], help='Which level of reconstruction? reco, gen')
argParser.add_argument('--luminosity',         action='store',     default=150, type=int, help='Luminosity for weighting the plots')
argParser.add_argument('--detector',           action='store',     default='CMS', nargs='?', choices=['CMS', 'ATLAS', 'phase2_CMS'], help='Which Delphes detector simulation?')
argParser.add_argument('--scale14TeV',         action='store_true', help='scale 13 TeV cross-sections to 14 Tev?')
argParser.add_argument('--additionalCardFile', action='store',     default='TopEFTCardFile.txt', help='Cardfile where additional uncertainties are taken from')
argParser.add_argument('--addNonPrompt',       action='store_true', help='add nonPrompt?')
argParser.add_argument('--addOthers',          action='store_true', help='add rare, ttW, ZZ, Zgamma?')
argParser.add_argument('--addBinNumberShift',  action='store',     default = 0, type=int, help = "which bin number does the region start in the additional card file?")
argParser.add_argument('--combineTTX',         action='store_true', help='combine tWZ and ttgamma to ttX')
argParser.add_argument('--parameters',         action='store',     default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--topEFTcolors',       action='store_true', help='use comparable colors to TopEFT?')
argParser.add_argument('--nonInfoSignal',      action='store_true', help='add nonInfo Signal?')
argParser.add_argument('--addUncertainties',   action='store_true', help='add uncertainties?')
argParser.add_argument('--noTheoryUnc',        action='store_true', help='add uncertainties?')

args = argParser.parse_args()

if len(args.parameters) < 2: args.parameters = None

#colors = { 'WZ':ROOT.kAzure-3, 'tWZ':ROOT.kGreen-2, 'tZq':ROOT.kCyan-9, 'ttgamma':ROOT.kRed+2, 'ttZ':ROOT.kWhite, 'nonInfo':ROOT.kRed-7 }
colors = { 'WZ':ROOT.kAzure-3, 'tWZ':ROOT.kGreen-2, 'tZq':ROOT.kCyan-9, 'ttgamma':ROOT.kRed+2, 'ttZ':ROOT.kOrange, 'nonInfo':ROOT.kRed-7 }
#if args.nonInfoSignal: colors['ttZ']=ROOT.kOrange#41#ROOT.kWhite

if args.addUncertainties:
    if not os.path.isfile('TTXPhenoCardFile.txt'):
        raise Exception('Cardfile TTXPhenoCardFile.txt not found! Cannot add uncertainties!')
    if args.noTheoryUnc: uncertainties = [ 'lumi', 'JES', 'btagging', 'mistagging', 'muonId', 'electronId', 'trigger_2016', 'nonprompt' ]
    else: uncertainties = [ 'lumi', 'JES', 'btagging', 'mistagging', 'muonId', 'electronId', 'trigger_2016', 'scale', 'scale_sig', 'PDF', 'PartonShower', 'nonprompt', 'WZ_xsec', 'ttX' ]
    #uncertainties = []

params = {}
if args.parameters is not None:
    coeffs = args.parameters[::2]
    str_vals = args.parameters[1::2]
    vals = list( map( float, str_vals ) )
    str_vals = list( map( int, vals ))
    str_vals = list( map( str, str_vals ))
    for i_param, (coeff, val, str_val) in enumerate(zip(coeffs, vals, str_vals)):
#        params[coeff] = val
        params[coeff] = int(val)
    signalLegendText = ""
    for i, par in enumerate(params.keys()):
        signalLegendText += par.replace('gamma','#gamma').replace('c','C_{').replace('p','#phi').replace('M','').replace('I','}^{[Im]') + '} = ' + str(params[par])
        if i!=len(params.keys())-1: signalLegendText += ', '


if args.level == 'gen':
    if args.process == 'ttZ':
        from TTXPheno.Analysis.regions import genttZRegions as regions
    elif args.process == 'ttgamma':
        from TTXPheno.Analysis.regions import genttgammaRegions as regions

    # Import additional functions/classes specified for the level of reconstruction
    from TTXPheno.Tools.cutInterpreterGen import cutInterpreter

elif args.level == 'reco':
    if args.process == 'ttZ':
        from TTXPheno.Analysis.regions import recottZRegions as regions
    elif args.process == 'ttgamma':
        from TTXPheno.Analysis.regions import recottgammaRegions as regions

    # Import additional functions/classes specified for the level of reconstruction
    from TTXPheno.Tools.cutInterpreterReco import cutInterpreter

# Import samples
sample_file     = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
loadedSamples   = imp.load_source( "samples", os.path.expandvars( sample_file ) )

ttXSample       = getattr( loadedSamples, args.sample + '_%s' %args.detector )
WZSample        = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights_%s' %args.detector )
#ttSample        = getattr( loadedSamples, 'fwlite_tt_full_LO_order2_15weights_%s' %args.detector )
#tWSample        = getattr( loadedSamples, 'fwlite_tW_LO_order2_15weights_%s' %args.detector )
tWZSample       = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights_%s' %args.detector )
tZqSample       = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights_%s' %args.detector )
#ZgammaSample    = getattr( loadedSamples, 'fwlite_Zgamma_LO_order2_15weights_%s' %args.detector )
ttgammaSample   = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights_%s' %args.detector )

bg = [ WZSample, tZqSample, tWZSample, ttgammaSample ]
nonInfo = []
if args.nonInfoSignal:
    ttZISRSample  = copy.deepcopy( ttXSample ) #select ttgamma events with isolated gamma from ISR (cat a2)
    ttZISRSample.name = 'fwlite_nonInfo_LO_order2_15weights_ref'
    nonInfo.append( ttZISRSample)

def checkReferencePoint( sample ):
    ''' check if sample is simulated with a reference point
    '''
    return pickle.load(file(sample.reweight_pkl))['ref_point'] != {}

# set selection string
selectionString      = cutInterpreter.cutString(args.selection)
selectionString_up   = selectionString.replace('nBTag','nBTag_JEC_up').replace('nrecoJet','nrecoJets_JEC_up')
selectionString_down = selectionString.replace('nBTag','nBTag_JEC_down').replace('nrecoJet','nrecoJets_JEC_down')

# somehow has to be separate from the next loop
if args.small:
    for s in [ttXSample] + bg + nonInfo:
        s.reduceFiles( to = 30 )

# configure samples
for i, s in enumerate([ttXSample] + bg + nonInfo):

    s.shortname = s.name.split('_')[1]
    s.event_factor = s.nEvents / float( s.chain.GetEntries() )
    s.xsecScaleFactor = s.xsec14 / s.xsec if args.scale14TeV else 1.
    s.weightInfo = WeightInfo( s.reweight_pkl )
    s.weightInfo.set_order( args.order )
    s.setSelectionString( selectionString )
    s.color = getattr(color, s.shortname) if args.topEFTcolors else colors[s.shortname] #color[i]

    if checkReferencePoint( s ):
        s.setWeightString( 'ref_lumiweight1fb*(%s)*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor), str(s.xsecScaleFactor) ) )
    else:
        s.setWeightString( 'lumiweight1fb*(%s)*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor), str(s.xsecScaleFactor) ) )

if args.nonInfoSignal:
    ttXSample.addSelectionString( "signalZ==1" ) #Z from gluon or top
    ttZISRSample.addSelectionString( "signalZ==0" ) #Z from ISR or else

hists = {}
Nbins = len(regions)
for s in [ttXSample] + bg + nonInfo:
    hists[s.shortname] = ROOT.TH1F(s.shortname,"", Nbins, 0, Nbins)
if args.nonInfoSignal:
    hists['nonInfo'] = ROOT.TH1F('nonInfo',"", Nbins, 0, Nbins)
if args.addNonPrompt:
    hists['nonPrompt'] = ROOT.TH1F('nonprompt',"", Nbins, 0, Nbins)
if args.addOthers:
    hists['rare'] = ROOT.TH1F('rare',"", Nbins, 0, Nbins)
    hists['ttW'] = ROOT.TH1F('ttW',"", Nbins, 0, Nbins)
    hists['ZZ'] = ROOT.TH1F('ZZ',"", Nbins, 0, Nbins)
    hists['ZG'] = ROOT.TH1F('ZG',"", Nbins, 0, Nbins)
if args.combineTTX:
    hists['ttX'] = ROOT.TH1F('ttX',"", Nbins, 0, Nbins)
if args.parameters is not None:
    hists['signal'] = ROOT.TH1F('signal',"", Nbins, 0, Nbins)
hists['SM'] = ROOT.TH1F('SM',"", Nbins, 0, Nbins)

hists['empty'] = ROOT.TH1F('empty',"", Nbins, 0, Nbins)
hists['empty'].legendText = ''
hists['empty'].style = styles.fillStyle( ROOT.kWhite, lineColor=ROOT.kWhite, errors=False, width=0 )

rate = {}

for i_region, region in enumerate(regions):
    # compute signal yield for this region (this is the final code)

    logger.info( "At region %s", region )

    rate[region] = {}

    ttX_coeffList                     = ttXSample.weightInfo.getCoeffListFromDraw( ttXSample, selectionString = region.cutString() )
    rate[region][ttXSample.shortname] = ttXSample.weightInfo.get_weight_yield( ttX_coeffList )

    for i_background, background in enumerate(bg):
        rate[region][background.shortname] = background.getYieldFromDraw( selectionString=region.cutString() )['val']

    if args.nonInfoSignal:
        nonInfo_coeffList       = ttZISRSample.weightInfo.getCoeffListFromDraw( ttZISRSample, selectionString = region.cutString() )
        rate[region]['nonInfo'] = ttZISRSample.weightInfo.get_weight_yield( nonInfo_coeffList )

    if args.addNonPrompt:
        # scale nonprompt observation value from Run2 to args.luminosity
        rate[region]['nonPrompt'] = getObservationValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'nonPromptDD' ) * float(args.luminosity) / 35.9

    if args.addOthers:
        rate[region]['rare'] = getObservationValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'rare' ) * float(args.luminosity) / 35.9
        rate[region]['ttW'] = getObservationValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'TTW' ) * float(args.luminosity) / 35.9
        rate[region]['ZZ'] = getObservationValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'ZZ' ) * float(args.luminosity) / 35.9
        rate[region]['ZG'] = getObservationValue( args.additionalCardFile, args.addBinNumberShift + i_region, 'ZG' ) * float(args.luminosity) / 35.9

    if args.parameters is not None:
        rate[region]['signal'] = ttXSample.weightInfo.get_weight_yield( ttX_coeffList, **params)
        if args.nonInfoSignal: rate[region]['signal'] += ttZISRSample.weightInfo.get_weight_yield( nonInfo_coeffList, **params )

for i_region, region in enumerate(regions):

    totalUncertainty = 0

    for i, s in enumerate( [ttXSample] + bg ):
        hists[s.shortname].SetBinContent(i_region+1, rate[region][s.shortname])
        hists[s.shortname].SetBinError(i_region+1,0)
        hists[s.shortname].legendText = s.shortname.replace('gamma', '#gamma')
#        hists[s.shortname].style = styles.fillStyle( s.color, lineColor=s.color if i!=0 else ROOT.kBlack, errors=False, width=0 if i!=0 else 2 )
        hists[s.shortname].style = styles.fillStyle( s.color, lineColor=ROOT.kBlack, errors=False, width=1 if i!=0 else 2 )

        if args.addUncertainties and s.name != ttXSample.name:
            for unc in uncertainties:
                totalUncertainty += (abs(getUncertaintyValue( 'TTXPhenoCardFile.txt', i_region, '_'.join(s.name.split('_')[1:3]), unc ) - 1) * rate[region][s.shortname])**2

    if args.nonInfoSignal:
        hists['nonInfo'].SetBinContent(i_region+1, rate[region]['nonInfo'])
        hists['nonInfo'].SetBinError(i_region+1,0)
        hists['nonInfo'].legendText = ttXSample.shortname.replace('gamma', '#gamma') + ' (non-info)'
        hists['nonInfo'].style = styles.fillStyle( colors['nonInfo'], lineColor=ROOT.kBlack, errors=False, width=1, fillStyle=3345, hatchesWidth=1, hatchesSpacing=None )
#        hists['nonInfo'].style = styles.fillStyle( colors['nonInfo'], lineColor=colors['nonInfo'], errors=False )
#        hists['nonInfo'].SetFillStyle(3345)
#        hists['nonInfo'].SetFillStyle(3005)
#        hists['nonInfo'].SetFillStyle(3544)
#        hists['nonInfo'].SetFillColor(ROOT.kWhite)
#        hists['nonInfo'].SetLineWidth(1)
#        hists['nonInfo'].SetLineColor(ROOT.kBlack)
#        ROOT.gStyle.SetHatchesLineWidth(2)
#        ROOT.gStyle.SetHatchesSpacing(1.1)

    if args.addNonPrompt:
        hists['nonPrompt'].SetBinContent(i_region+1, rate[region]['nonPrompt'])
        hists['nonPrompt'].SetBinError(i_region+1,0)
        hists['nonPrompt'].legendText = 'non-prompt     '
#        hists['nonPrompt'].style = styles.fillStyle( getattr(color, 'nonprompt'), lineColor=getattr(color, 'nonprompt'), errors=False )
        hists['nonPrompt'].style = styles.fillStyle( getattr(color, 'nonprompt'), lineColor=ROOT.kBlack, errors=False, width=1 )

        if args.addUncertainties:
            for unc in uncertainties:
                totalUncertainty += (abs(getUncertaintyValue( 'TTXPhenoCardFile.txt', i_region, 'nonPrompt', unc ) - 1) * rate[region]['nonPrompt'])**2

    if args.addOthers:
        hists['rare'].SetBinContent(i_region+1, rate[region]['rare'])
        hists['rare'].SetBinError(i_region+1,0)
        hists['rare'].legendText = 'rare'
#        hists['rare'].style = styles.fillStyle( getattr(color, 'rare'), lineColor=getattr(color, 'rare'), errors=False )
        hists['rare'].style = styles.fillStyle( getattr(color, 'rare'), lineColor=ROOT.kBlack, errors=False, width=1 )

        hists['ttW'].SetBinContent(i_region+1, rate[region]['ttW'])
        hists['ttW'].SetBinError(i_region+1,0)
        hists['ttW'].legendText = 'ttW'
#        hists['ttW'].style = styles.fillStyle( getattr(color, 'ttW'), lineColor=getattr(color, 'ttW'), errors=False )
        hists['ttW'].style = styles.fillStyle( getattr(color, 'ttW'), lineColor=ROOT.kBlack, errors=False, width=1 )

        hists['ZZ'].SetBinContent(i_region+1, rate[region]['ZZ'])
        hists['ZZ'].SetBinError(i_region+1,0)
        hists['ZZ'].legendText = 'ZZ'
#        hists['ZZ'].style = styles.fillStyle( getattr(color, 'ZZ'), lineColor=getattr(color, 'ZZ'), errors=False )
        hists['ZZ'].style = styles.fillStyle( getattr(color, 'ZZ'), lineColor=ROOT.kBlack, errors=False, width=1 )

        hists['ZG'].SetBinContent(i_region+1, rate[region]['ZG'])
        hists['ZG'].SetBinError(i_region+1,0)
        hists['ZG'].legendText = 'Z#gamma'
#        hists['ZG'].style = styles.fillStyle( getattr(color, 'ZG'), lineColor=getattr(color, 'ZG'), errors=False )
        hists['ZG'].style = styles.fillStyle( getattr(color, 'ZG'), lineColor=ROOT.kBlack, errors=False, width=1 )

        if args.addUncertainties:
            for unc in uncertainties:
                for p in ['rare', 'ttW', 'ZZ', 'ZG']:
                    totalUncertainty += (abs(getUncertaintyValue( 'TTXPhenoCardFile.txt', i_region, p, unc ) - 1) * rate[region][p])**2

    if args.combineTTX:
        hists['ttX'].SetBinContent(i_region+1, rate[region]['tWZ'] + rate[region]['ttgamma'] + rate[region]['tZq'])
        hists['ttX'].SetBinError(i_region+1,0)
        hists['ttX'].legendText = 't(t)X'#=tZq,tWZ,tt#gamma'
#        hists['ttX'].style = styles.fillStyle( getattr(color, 'ttX'), lineColor=getattr(color, 'ttX'), errors=False )
        hists['ttX'].style = styles.fillStyle( getattr(color, 'ttX'), lineColor=ROOT.kBlack, errors=False, width=1 )

    if args.parameters is not None:
        hists['signal'].SetBinContent(i_region+1, rate[region]['signal'])
        hists['signal'].SetBinError(i_region+1,0)
        hists['signal'].legendText = signalLegendText
#        hists['signal'].style = styles.lineStyle( ROOT.kRed, width=2, dashed=True )
        hists['signal'].style = styles.lineStyle( ROOT.kRed, width=2, dashed=False )

    hists['SM'].SetBinContent(i_region+1, rate[region][ttXSample.shortname] + rate[region]['nonInfo'] if args.nonInfoSignal else rate[region][ttXSample.shortname])
    hists['SM'].SetBinError(i_region+1,0)
    hists['SM'].legendText = 'SM'
    hists['SM'].style = styles.lineStyle( ROOT.kBlack, width=2 )

    if args.addUncertainties:
        for unc in uncertainties:
            totalUncertainty += (abs(getUncertaintyValue( 'TTXPhenoCardFile.txt', i_region, 'signal', unc ) - 1) * (rate[region][ttXSample.shortname] + rate[region]['nonInfo'] if args.nonInfoSignal else rate[region][ttXSample.shortname]))**2

    totalUncertainty = sqrt(totalUncertainty)
    hists['SM'].SetBinError(i_region+1,totalUncertainty)

def drawDivisions(regions):
    min = 0.15
    max = 0.95
    diff = (max-min) / len(regions)
    lines = []
    lines2 = []
    line = ROOT.TLine()
#   line.SetLineColor(38)
    line.SetLineWidth(1)
    line.SetLineStyle(9)
#    lines = [ (min+3*i*diff,  0.023, min+3*i*diff, 0.93) if min+3*i*diff<0.74 else (min+3*i*diff,  0.023, min+3*i*diff, 0.52) for i in range(1,10) ]
    lines = [ (min+3*i*diff,  0.08, min+3*i*diff, 0.93) if min+3*i*diff<0.74 else (min+3*i*diff,  0.08, min+3*i*diff, 0.93) for i in range(1,10) ]
    return [line.DrawLineNDC(*l) for l in lines[:3]] + [tex.DrawLatex(*l) for l in []] + [tex2.DrawLatex(*l) for l in lines2]

def setBinLabels( hist ):
    for i in range(1, hist.GetNbinsX()+1):
        if i < 16:
            hist.GetXaxis().SetBinLabel(i, "%s"%i)
        else:
            hist.GetXaxis().SetBinLabel(i, "%s"%(i-15))

def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    tex.SetTextFont(42)
    lines = [
      (0.15, 0.95, '#bf{CMS Phase-2} #it{Simulation Preliminary}'),
      (0.76, 0.95, '%i ab{}^{-1} (%s TeV)'% ( int(args.luminosity/1000.), '14' if args.scale14TeV else '13' ) )
#      (0.15, 0.95, ' '.join(args.processFile.split('_')[:2]) + '(' + args.detector + ')'),
#      (offset, 0.95, '%3.1f fb{}^{-1} @ 13 TeV%s'% ( float(args.luminosity), titleAddon) )
    ]
    return [tex.DrawLatex(*l) for l in lines]

def drawLabels( regions ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.032)
    if len(regions)>12:
        tex.SetTextAngle(90)
    else:
        tex.SetTextAngle(0)
    tex.SetTextAlign(12) # align right
    min = 0.15
    max = 0.95
    diff = (max-min) / len(regions)
    y_pos = 0.6 if len(regions)>12 else 0.87
    x_pos = 0.90 if len(regions)>12 else 0.18
    if len(regions) > 12:
        lines =  [(min+(3*i+x_pos)*diff, y_pos,  r.texStringForVar('recoZ_pt'))   for i, r in enumerate(regions[:-3][::3])]
    else:
        lines =  [(min+(3*i+x_pos)*diff + 0.017 if i==0 or i==len(regions[::3])-1 else min+(3*i+x_pos)*diff, y_pos,  r.texStringForVar('recoZ_pt'))   for i, r in enumerate(regions[::3])]
    return [tex.DrawLatex(*l) for l in lines]


for hist in hists.keys():
    setBinLabels(hists[hist])

bkgHists = []
for s in bg:
    ind = s.shortname
    if args.combineTTX and (ind=='tWZ' or ind=='ttgamma' or ind=='tZq'): continue
    if ind=="ttgamma": bkgHists.append(hists['empty'])
    bkgHists.append(hists[ind])
if args.addOthers:
    bkgHists.append(hists['ZZ'])
if args.addNonPrompt:
    if args.combineTTX: bkgHists.append(hists['nonPrompt'])
    else: bkgHists = [hists['nonPrompt']] + bkgHists
if args.combineTTX:
    bkgHists.append(hists['ttX'])
#if args.addNonPrompt:
#    bkgHists = [hists['nonPrompt']] + bkgHists
if args.addOthers:
    bkgHists += [hists['ttW'], hists['ZG'], hists['rare']]


SM = [ [ hists[ttXSample.shortname], hists['nonInfo'] ] ] if args.nonInfoSignal else [[hists[ttXSample.shortname]]]

if args.parameters is not None:
    plots = [ [hists['signal']], [hists['SM']] ] + SM + [ bkgHists ]
#    plots = [ [hists['signal']] ] + SM + [ bkgHists ]
else:
    plots = [ [hists['SM']] ] + SM + [ bkgHists ]
#    plots = SM + [ bkgHists ]

#histos =  bkgHists  + [hists["total"]]
#if options.signal:
#    plots = [ [hists['BSM']], bkgHists, [hists['observed']] ]
#else:
#    plots = [ bkgHists, [hists['observed']]]
#if subDir:
#    subDir = "%s_"%subDir

plot_directory_ = os.path.join(\
    plot_directory,
    '%s_%s'%(args.level, args.version),
    args.detector,
    args.sample,
    'backgrounds',
    'regions_small' if args.small else 'regions',
    args.selection.replace('nRun2bjet', 'nbjet'),
    '_'.join(args.parameters).rstrip('0').replace('-','m').replace('.','p') if args.parameters is not None else 'SM'
)

if not os.path.isdir( plot_directory_ ):
    os.makedirs( plot_directory_ )

plotNameList = ['regions', 'lumi'+str(args.luminosity), "14TeV" if args.scale14TeV else "13TeV"]
if args.combineTTX: plotNameList.append('ttXcombined')
if args.addNonPrompt: plotNameList.append('nonPrompt')
if args.addOthers: plotNameList.append('others')
if args.noTheoryUnc: plotNameList.append('noTheoryUnc')
plotName = '_'.join(plotNameList)

plot = Plot.fromHisto(plotName, plots, texX = "Signal Regions", texY = "Number of Events" )

#bgHistIndex = 1 if args.nonInfoSignal else 1
for bgHisto in plot.histos[-1]:
    for signalHisto in plot.histos[:-1]:# + plot.histos[bgHistIndex+1:]:
        signalHisto[-1].Add(bgHisto)

boxes = []
ratio_boxes = []
tot = 0
tot_val = 0
for ib in range(1, 1 + hists['SM'].GetNbinsX() ):
    val = hists['SM'].GetBinContent(ib)
    if val<0: continue
    sys = hists['SM'].GetBinError(ib)
    sys_rel = sys/val
    # uncertainty box in main histogram
    box = ROOT.TBox( hists['SM'].GetXaxis().GetBinLowEdge(ib),  max([0.006, val-sys]), hists['SM'].GetXaxis().GetBinUpEdge(ib), max([0.006, val+sys]) )
    box.SetLineColor(ROOT.kBlack)
    box.SetFillStyle(3444)
    box.SetFillColor(ROOT.kBlack)

    print 'region', ib, 'yield', val, 'unc', sys, 'rel unc', sys_rel, 'stat+sys', sqrt(sys**2+val), 'stat+sys rel', sqrt(sys**2+val)/val
    tot += sys**2+val 
    tot_val += val
    # uncertainty box in ratio histogram
    r_box = ROOT.TBox( hists['SM'].GetXaxis().GetBinLowEdge(ib),  max(0.1, 1-sys_rel), hists['SM'].GetXaxis().GetBinUpEdge(ib), min(1.9, 1+sys_rel) )
    r_box.SetLineColor(ROOT.kBlack)
    r_box.SetFillStyle(3444)
    r_box.SetFillColor(ROOT.kBlack)

    boxes.append( box )
    hists['SM'].SetBinError(ib, 0)
    ratio_boxes.append( r_box )

print 'tot', sqrt(tot), 'tot_rel', sqrt(tot)/tot_val

def histmodification(log):
    def histmod(h):
        h.GetXaxis().SetTitleOffset( 1.06 )
        h.GetYaxis().SetTitleOffset( 1.08 if log else 1.5 )

        h.GetXaxis().SetTitleSize( 0.042 )
        h.GetYaxis().SetTitleSize( 0.042 )

        h.GetXaxis().SetLabelSize( 0.06 )
        h.GetYaxis().SetLabelSize( 0.04 )
    return histmod

def legendmodification(l):
    l.SetTextSize(.032)

for logY in [True, False]:

    if logY:
        ymax = 2000000
    else:
        ymax = 12000

    plotting.draw(
        plot,
        plot_directory = os.path.join(plot_directory_, 'log' if logY else 'lin'),
        logX = False, logY = logY, sorting = False,
        legend = ( (0.55,0.6, 0.925, 0.85), 2 ),
#        legend = (0.75,0.49, 0.95, 0.85),
        widths = {'x_width':750, 'y_width':600},
        yRange = (0.7,ymax),
#       ratio = {'yRange': (0.6, 1.4), 'drawObjects':boxes},
        drawObjects = drawObjects() + drawDivisions( regions ) + drawLabels( regions ) + boxes ,
        histModifications = [histmodification(logY)],
        legendModifications = [legendmodification],
        copyIndexPHP = True,
    )


