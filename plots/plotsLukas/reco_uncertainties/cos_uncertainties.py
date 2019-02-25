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

from TTXPheno.Analysis.regions import recottZRegionsCosPTZ200 as regions
from TTXPheno.Tools.cutInterpreterReco import cutInterpreter

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
argParser.add_argument('--selection',          action='store',     default='lepSel3-onZ-njet3p-nbjet1p-Zpt200-leptonIso3', help="Specify cut.")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--level',              action='store',     default='reco', nargs='?', choices=['reco', 'gen'], help='Which level of reconstruction? reco, gen')
argParser.add_argument('--luminosity',         action='store',     default=3000, type=int, help='Luminosity for weighting the plots')
argParser.add_argument('--detector',           action='store',     default='phase2_CMS', nargs='?', choices=['CMS', 'ATLAS', 'phase2_CMS'], help='Which Delphes detector simulation?')
argParser.add_argument('--scale14TeV',         action='store_true', help='scale 13 TeV cross-sections to 14 Tev?')
argParser.add_argument('--cardFile',           action='store',     default='Cos_SM_cardfile.txt', help='Cardfile where additional uncertainties are taken from')
argParser.add_argument('--parameters',         action='store',     default = [], type=str, nargs='+', help = "argument parameters")
argParser.add_argument('--nonInfoSignal',      action='store_true', help='add nonInfo Signal?')
argParser.add_argument('--addUncertainties',   action='store_true', help='add uncertainties?')
argParser.add_argument('--shape',              action='store_true', help='Scale lumi only?')

args = argParser.parse_args()

if len(args.parameters) < 2: args.parameters = None

var = "Z_cosThetaStar5"
colors = { 'WZ':ROOT.kAzure-3, 'tWZ':ROOT.kGreen-2, 'tZq':ROOT.kCyan-9, 'ttgamma':ROOT.kRed+2, 'ttZ':ROOT.kWhite, 'nonInfo':ROOT.kRed-7 }
bsmColors = [ ROOT.kRed+1, ROOT.kGreen+2, ROOT.kOrange+1, ROOT.kViolet+9, ROOT.kSpring-7, ROOT.kRed+2,  ROOT.kPink-9, ROOT.kBlue,  ROOT.kRed-7, ROOT.kRed-10, ROOT.kRed+3,  ROOT.kGreen-7, ROOT.kGreen-10]

if args.addUncertainties:
    if not os.path.isfile(args.cardFile):
        raise Exception('Cardfile not found! Cannot add uncertainties!')
    uncertainties = [ 'lumi', 'JES', 'btagging', 'mistagging', 'muonId', 'electronId', 'trigger_2016', 'scale', 'scale_sig', 'PDF', 'PartonShower', 'WZ_xsec', 'ttX' ]

params = []
if args.parameters is not None:
    coeffs = args.parameters[::2]
    vals = list( map( float, args.parameters[1::2] ))
#    str_vals = list( map( "{:3.1f}".format, vals ))
    str_vals = list( map( int, vals ))
    str_vals = list( map( str, str_vals ))
    for i_param, (coeff, val, str_val) in enumerate(zip(coeffs, vals, str_vals)):
        params.append( {
            'legendText': ' '.join([coeff,str_val]).replace('gamma','#gamma').replace('c','C_{').replace(' ','} = ').replace('p','#phi').replace('M','').replace('I','}^{[Im]') + '  ',
#            'legendText': ' '.join([coeff]).replace('gamma','#gamma').replace('c','C_{').replace('p','#phi').replace('M','').replace('I','}^{[Im]') + '}  ',
            'WC' : { coeff:val },
            'name' : coeff ,
            'color' : bsmColors[i_param],
            } )

# Import samples
sample_file     = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
loadedSamples   = imp.load_source( "samples", os.path.expandvars( sample_file ) )

ttXSample       = getattr( loadedSamples, args.sample + '_%s' %args.detector )
WZSample        = getattr( loadedSamples, 'fwlite_WZ_lep_LO_order2_15weights_%s' %args.detector )
tWZSample       = getattr( loadedSamples, 'fwlite_tWZ_LO_order2_15weights_%s' %args.detector )
tZqSample       = getattr( loadedSamples, 'fwlite_tZq_LO_order2_15weights_%s' %args.detector )
ttgammaSample   = getattr( loadedSamples, 'fwlite_ttgamma_bg_LO_order2_15weights_%s' %args.detector )

bg = [ tZqSample, WZSample, tWZSample, ttgammaSample ]
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
        s.reduceFiles( to = 5 )

# configure samples
for i, s in enumerate([ttXSample] + bg + nonInfo):

    s.shortname = s.name.split('_')[1]
    s.event_factor = s.nEvents / float( s.chain.GetEntries() )
    s.xsecScaleFactor = s.xsec14 / s.xsec if args.scale14TeV else 1.
    s.weightInfo = WeightInfo( s.reweight_pkl )
    s.weightInfo.set_order( args.order )
    s.setSelectionString( selectionString )
    s.color = colors[s.shortname] #color[i]

    if checkReferencePoint( s ):
        s.setWeightString( 'ref_lumiweight1fb*(%s)*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor), str(s.xsecScaleFactor) ) )
    else:
        s.setWeightString( 'lumiweight1fb*(%s)*(%s)*(%s)'%( str(args.luminosity), str(s.event_factor), str(s.xsecScaleFactor) ) )

if args.nonInfoSignal:
    ttXSample.addSelectionString( "signalZ==1" ) #Z from gluon or top
    ttZISRSample.addSelectionString( "signalZ==0" ) #Z from ISR or else

hists = {}
Nbins = len(regions)
maxval = 1
minval = -1
for s in [ttXSample] + bg + nonInfo:
    hists[s.shortname] = ROOT.TH1F(s.shortname,"", Nbins, minval, maxval)
if args.nonInfoSignal:
    hists['nonInfo'] = ROOT.TH1F('nonInfo',"", Nbins, minval, maxval)
if args.parameters is not None:
    for param in params:
        hists['signal_%s'%param['name']] = ROOT.TH1F('signal_' + param['name'],"", Nbins, minval, maxval)
#hists['SM'] = ROOT.TH1F('SM',"", Nbins, 0, Nbins)

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

    if args.parameters is not None:
        for param in params:
            rate[region]['signal_' + param['name']] = ttXSample.weightInfo.get_weight_yield( ttX_coeffList, **param['WC'])
            if args.nonInfoSignal: rate[region]['signal_' + param['name']] += ttZISRSample.weightInfo.get_weight_yield( nonInfo_coeffList, **param['WC'] )

for i_region, region in enumerate(regions):

    totalUncertainty = 0

    for i_s, s in enumerate([ttXSample] + bg):
        hists[s.shortname].SetBinContent(i_region+1, rate[region][s.shortname])
        hists[s.shortname].SetBinError(i_region+1,0)
        hists[s.shortname].legendText = s.shortname.replace('gamma', '#gamma')
        hists[s.shortname].style = styles.fillStyle( s.color, lineColor=ROOT.kBlack, errors=False, width=1 if i_s != 0 else 3 )

        if args.addUncertainties and s.name != ttXSample.name:
            for unc in uncertainties:
                totalUncertainty += (abs(getUncertaintyValue( args.cardFile, i_region, '_'.join(s.name.split('_')[1:3]), unc ) - 1) * rate[region][s.shortname])**2

    if args.nonInfoSignal:
        hists['nonInfo'].SetBinContent(i_region+1, rate[region]['nonInfo'])
        hists['nonInfo'].SetBinError(i_region+1,0)
        hists['nonInfo'].legendText = ttXSample.shortname.replace('gamma', '#gamma') + ' (non-info)'
        hists['nonInfo'].style = styles.fillStyle( colors['nonInfo'], lineColor=ROOT.kBlack, errors=False, width=1, fillStyle=3645, hatchesWidth=1, hatchesSpacing=None )
#        hists['nonInfo'].style = styles.fillStyle( colors['nonInfo'], lineColor=ROOT.kBlack, errors=False, width=1 )
#        hists['nonInfo'].SetFillStyle(3345)
#        hists['nonInfo'].SetFillStyle(3544)
#        hists['nonInfo'].SetFillColor(ROOT.kWhite)
#        ROOT.gStyle.SetHatchesLineWidth(3)

    if args.parameters is not None:
        for i_param, param in enumerate(params):
            hists['signal_' + param['name']].SetBinContent(i_region+1, rate[region]['signal_' + param['name']])
            hists['signal_' + param['name']].SetBinError(i_region+1,0)
            hists['signal_' + param['name']].legendText = param['legendText']
            hists['signal_' + param['name']].style = styles.lineStyle( param['color'], width=2, dashed=False )

    if args.addUncertainties:
        for unc in uncertainties:
            totalUncertainty += (abs(getUncertaintyValue( args.cardFile, i_region, 'signal', unc ) - 1) * (rate[region][ttXSample.shortname] + rate[region]['nonInfo'] if args.nonInfoSignal else rate[region][ttXSample.shortname]))**2

    totalUncertainty = sqrt(totalUncertainty)
    hists[ttXSample.shortname].SetBinError(i_region+1,totalUncertainty)

def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.042)
    tex.SetTextAlign(11) # align right
    tex.SetTextFont(42)
    lines = [
      (0.02, 0.95, '#bf{CMS Phase-2} #it{Simulation Preliminary}'),
      (0.69, 0.95, '%i ab{}^{-1} (%s TeV)'% ( int(args.luminosity/1000.), '14' if args.scale14TeV else '13' ) )
    ]
    return [tex.DrawLatex(*l) for l in lines]


bkgHists = []
for s in bg:
    ind = s.shortname
    bkgHists.append(hists[ind])

SM = [ [ hists[ttXSample.shortname], hists['nonInfo'] ] ] if args.nonInfoSignal else [[hists[ttXSample.shortname]]]

if args.parameters is not None:
    plots = [ [hists['signal_' + param['name']]] for param in list(reversed(params)) ] + SM + [ bkgHists ]
else:
    plots = SM + [ bkgHists ]

plot_directory_ = os.path.join(\
    plot_directory,
    '%s_%s'%(args.level, args.version),
    args.detector,
    args.sample,
    'backgrounds',
    'uncertainties_small' if args.small else 'uncertainties',
    args.selection.replace('nRun2bjet', 'nbjet'),
    '_'.join(args.parameters).rstrip('0').replace('-','m').replace('.','p') if args.parameters is not None else 'SM'
)

if not os.path.isdir( plot_directory_ ):
    os.makedirs( plot_directory_ )

plotNameList = [var]
plotName = '_'.join(plotNameList)

plot = Plot.fromHisto(plotName, plots, texX = "cos(#theta_{Z}*)", texY = "Number of Events" )

#bgHistIndex = 1 if args.nonInfoSignal else 1
for bgHisto in plot.histos[-1]:
    for signalHisto in plot.histos[:-1]:# + plot.histos[bgHistIndex+1:]:
        signalHisto[-1].Add(bgHisto)

boxes = []
ratio_boxes = []
tot = 0
tot_val = 0
for ib in range(1, 1 + hists[ttXSample.shortname].GetNbinsX() ):
    val = hists[ttXSample.shortname].GetBinContent(ib) + hists['nonInfo'].GetBinContent(ib)
    if val<0: continue
    sys = hists[ttXSample.shortname].GetBinError(ib)
    sys_rel = sys/val
    # uncertainty box in main histogram
    box = ROOT.TBox( hists[ttXSample.shortname].GetXaxis().GetBinLowEdge(ib),  max([0.006, val-sys]), hists[ttXSample.shortname].GetXaxis().GetBinUpEdge(ib), max([0.006, val+sys]) )
    box.SetLineColor(ROOT.kBlack)
    box.SetFillStyle(3444)
    box.SetFillColor(ROOT.kBlack)

    # uncertainty box in ratio histogram
    r_box = ROOT.TBox( hists[ttXSample.shortname].GetXaxis().GetBinLowEdge(ib),  max(0.1, 1-sys_rel), hists[ttXSample.shortname].GetXaxis().GetBinUpEdge(ib), min(1.9, 1+sys_rel) )
    r_box.SetLineColor(ROOT.kBlack)
    r_box.SetFillStyle(3444)
    r_box.SetFillColor(ROOT.kBlack)

    boxes.append( box )
    hists[ttXSample.shortname].SetBinError(ib, 0)
    ratio_boxes.append( r_box )

def histmodification(log):
    def histmod(h):
        h.GetXaxis().SetTitleOffset( 1.06 )
        h.GetYaxis().SetTitleOffset( 1.08 if log else 1.6 )

        h.GetXaxis().SetTitleSize( 0.045 )
        h.GetYaxis().SetTitleSize( 0.045 )

        h.GetXaxis().SetLabelSize( 0.04 )
        h.GetYaxis().SetLabelSize( 0.04 )

    return histmod

def ratiomodification(h):
    h.GetXaxis().SetTitleOffset( 1.25 )
    h.GetYaxis().SetTitleOffset( 0.65 )

    h.GetXaxis().SetTitleSize( 0.11 )
    h.GetYaxis().SetTitleSize( 0.11 )

    h.GetXaxis().SetLabelSize( 0.10 )
    h.GetYaxis().SetLabelSize( 0.10 )

def legendmodification(l):
    l.SetTextSize(.04)

#ROOT.gStyle.SetPadLeftMargin(0.14)
#ROOT.gStyle.SetPadRightMargin(0.1)

for logY in [True, False]:

    plotting.draw(
        plot,
        plot_directory = os.path.join(plot_directory_, 'log' if logY else 'lin'),
        logX = False, logY = logY, sorting = False,
        yRange = (0.01,'auto'),
        legend = ((0.18,0.75,0.91,0.9), 4),
        widths = {'x_width':500, 'y_width':500},
        scaling = {i:len(params) for i in range(0,len(params))} if args.shape else {},
        drawObjects = drawObjects() + boxes, #+ ratio_boxes,
        ratio = {'yRange': (0.3, 1.7), 'histos':[(i,len(params)) for i in range(0,len(params))], 'texY':'BSM / SM', 'drawObjects':ratio_boxes},
        histModifications = [histmodification(logY)],
        ratioModifications = [ratiomodification],
        legendModifications = [legendmodification],
        copyIndexPHP = True,
    )

