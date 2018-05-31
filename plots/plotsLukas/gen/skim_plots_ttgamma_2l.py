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
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR 
from TTXPheno.Tools.WeightInfo           import WeightInfo
from TTXPheno.Tools.cutInterpreter       import cutInterpreter

# Import samples
from TTXPheno.samples.benchmarks         import *

# Import helpers
from plot_helpers                        import *

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='gen')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttgamma_LO_order3_8weights')
argParser.add_argument('--order',              action='store',      default=3)
argParser.add_argument('--selection',          action='store',      default='gammapt40-nlep2p-njet2p-nbjet1p', help="Specify cut.")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--scaleLumi',          action='store_true', help='Scale lumi only??')
argParser.add_argument('--parameters',         action='store',      default = ['ctZI', '2'], type=str, nargs='+', help = "argument parameters")

args = argParser.parse_args()

#
# Logger
#
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(    args.logLevel, logFile = None )
logger_rt = logger_rt.get_logger( args.logLevel, logFile = None )

# Make subdirectory
subDirectory = []
if args.scaleLumi:  subDirectory.append("shape")
else:               subDirectory.append("lumi")
if args.small:      subDirectory.append("small")
subDirectory = '_'.join( subDirectory )

# Import samples
#sample = fwlite_ttZ_ll_LO_order3_8weights 
sample_file = "$CMSSW_BASE/python/TTXPheno/samples/benchmarks.py"
samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
sample = getattr( samples, args.sample )

# Polynomial parametrization
w = WeightInfo(sample.reweight_pkl)
w.set_order(int(args.order))

# Parameters
params = [  
    {'legendText':'SM', 'WC':{}, 'color':ROOT.kBlack},
   ] 

colors = [ ROOT.kMagenta+1, ROOT.kOrange, ROOT.kBlue, ROOT.kCyan+1, ROOT.kGreen+1, ROOT.kRed, ROOT.kViolet, ROOT.kYellow+2 ]

coeffs = args.parameters[::2]
str_vals = args.parameters[1::2]
vals   = list( map( float, str_vals ) )
for i_param, (coeff, val, str_val) in enumerate(zip(coeffs, vals, str_vals)):
    params.append( { 
        'legendText': ' '.join([coeff,str_val]),
        'WC'        : { coeff:val },
        'color'     : colors[i_param], 
        })

# Make stack and weight
stack  = Stack(*[ [ sample ] for param in params ] )
weight = [ [ w.get_weight_func( **param['WC'] ) ] for param in params ]

def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'data' if hasData else ""), 
      #(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
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
        args.plot_directory, 
        sample.name, 
        subDirectory, 
        args.selection, 
        '_'.join(args.parameters).rstrip('0').replace('-','m').replace('.','p'),
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
        ratio = None, #{'yRange':(0.1,1.9)} if not args.noData else None,
        logX = False, logY = log, sorting = True,
        #yRange = (0.03, "auto") if log else (0., "auto"),
        #scaling = {i:0 for i in range(1, len(params))} if args.scaleLumi else {}, #Scale BSM shapes to SM (first in list)
        legend =  ( (0.17,0.9-0.05*sum(map(len, l_plot.histos))/3,1.,0.9), 3),
        drawObjects = drawObjects( ),
        copyIndexPHP = True,
    )

    # plot the plots
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = None, #{'yRange':(0.1,1.9)} if not args.noData else None,
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0., "auto"),
	    scaling = {i:0 for i in range(1, len(params))} if args.scaleLumi else {}, #Scale BSM shapes to SM (first in list)
	    legend = ( (0.17,0.9-0.05*sum(map(len, plot.histos))/3,1.,0.9), 3),
	    drawObjects = drawObjects( ),
        copyIndexPHP = True,
      )

#
# Read variables and sequences
#

read_variables = [
    "GenMet_pt/F", "GenMet_phi/F", 
    "nGenJet/I", "GenJet[pt/F,eta/F,phi/F,matchBParton/I]", 
    "nGenLep/I", "GenLep[pt/F,eta/F,phi/F,pdgId/I,motherPdgId/I]", 
    "ntop/I", "top[pt/F,eta/F,phi/F]", 
    "gamma_pt/F", "gamma_eta/F", "gamma_phi/F", "gamma_mass/F",
]
read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

logger.info( "Translating cut %s to %s", args.selection, cutInterpreter.cutString(args.selection) )
sample.setSelectionString( cutInterpreter.cutString(args.selection) )
sample.style = styles.lineStyle(ROOT.kBlue)

if args.small:
    for sample in stack.samples:
        sample.reduceFiles( to = 1 )


#sequence functions
sequence = []

def makeJets( event, sample ):
    ''' Add a list of filtered jets to the event
    '''
    # load jets
    event.jets = getCollection( event, 'GenJet', ['pt', 'eta', 'phi', 'matchBParton'], 'nGenJet' )

    # filter, pre-selection requires 3 leptons (no default leptons necessary)
    event.jets = list( filter( lambda j:isGoodJet( j ), event.jets ) )

    # sort
    event.jets = sorted( event.jets, key=lambda k: -k['pt'] )

    # Add extra vectors
    for p in event.jets:
        addTransverseVector( p )
        addTLorentzVector( p )

    # True B's
    event.trueBjets = list( filter( lambda j: j['matchBParton'], event.jets ) ) 
    event.trueNonBjets = list( filter( lambda j: not j['matchBParton'], event.jets ) )

    # Mimick b reconstruction ( if the trailing b fails acceptance, we supplement with the leading non-b jet ) 
    event.bj0, event.bj1 = ( event.trueBjets + event.trueNonBjets + [NanJet(), NanJet()] )[:2] 
    
sequence.append( makeJets )


def makeGamma( event, sample ):
    ''' Make a gamma vector to facilitate further calculations
    '''
    event.gamma_unitVec2D = UnitVectorT2( event.gamma_phi )
    event.gamma_vec4D     = ROOT.TLorentzVector()
    event.gamma_vec4D.SetPtEtaPhiM( event.gamma_pt, event.gamma_eta, event.gamma_phi, event.gamma_mass )
    event.gamma_unitVec3D = event.gamma_vec4D.Vect().Unit()

sequence.append( makeGamma )


def makeLeps( event, sample ):
    ''' Add a list of filtered leptons to the event
    '''
    # load leps
    event.leps = getCollection( event, 'GenLep', ['pt', 'eta', 'phi', 'pdgId', 'motherPdgId'], 'nGenLep' )

    # filter, pre-selection requires 3 leptons (no default leptons necessary)
    event.leps = list( filter( lambda l:isGoodLepton( l ), event.leps ) )

    # Cross-cleaning: remove leptons that overlap with a jet within 0.4
    event.leps = list(filter( lambda l: min( [ deltaR(l, j) for j in event.jets ] + [999] ) > 0.4 , event.leps ))

    # sort
    event.leps = sorted( event.leps, key=lambda l: -l['pt'] )

    # Add extra vectors
    for p in event.leps:
        addTransverseVector( p )
        addTLorentzVector( p )

    # 2l 
    event.l0, event.l1 = ( event.leps + [NanLepton(), NanLepton()] )[:2] 
    
    # We may loose some events by cross-cleaning or by thresholds.
    event.passing_2lep = len(event.leps)>=2 and event.l0['pdgId']*event.l1['pdgId'] > 0.

sequence.append( makeLeps )


def makeObservables( event, sample):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.deltaPhi_bb = deltaPhi( event.bj0['phi'], event.bj1['phi'] )
    event.deltaR_bb = deltaR( event.bj0, event.bj1 )

    # double l kinematic
    event.deltaPhi_ll = deltaPhi( event.l0['phi'], event.l1['phi'] )
    event.deltaR_ll = deltaR( event.l0, event.l1 )

    # signed lepton pt, Nan if len(event.lepsNotFromZ == 0)
    event.getl0chargept = event.l0['pt'] if event.l0['pdgId']>0 else -event.l0['pt']
    event.getl1chargept = event.l1['pt'] if event.l1['pdgId']>0 else -event.l1['pt']

sequence.append( makeObservables )


# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, addOverFlowBin=None)

  
plots = []

plots.append(Plot( name = "gamma_pt",
  texX = 'p_{T}(#gamma) [GeV]', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.gamma_pt if event.passing_2lep else float('nan'),
  binning=[20,0,400],
))

plots.append(Plot( name = "gamma_mass",
  texX = 'm(#gamma) [GeV]', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.gamma_mass if event.passing_2lep else float('nan'),
  binning=[20,-5,5],
))

plots.append(Plot( name = "b0_pt",
  texX = 'p_{T}(b_{0}) [GeV]', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.bj0['pt'] if event.passing_2lep else float('nan'),
  binning=[20,0,400],
))

plots.append(Plot( name = "b1_pt",
  texX = 'p_{T}(b_{1}) [GeV]', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.bj1['pt'] if event.passing_2lep else float('nan'),
  binning=[20,0,400],
))

plots.append(Plot( name = "b0_eta",
  texX = '#eta(b_{0})', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.bj0['eta'] if event.passing_2lep else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = "b1_eta",
  texX = '#eta(b_{1})', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.bj1['eta'] if event.passing_2lep else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = "b0_phi",
  texX = '#phi(b_{0})', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.bj0['phi'] if event.passing_2lep else float('nan'),
  binning=[20,pi,pi],
))

plots.append(Plot( name = "b1_phi",
  texX = '#phi(b_{1})', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.bj1['phi'] if event.passing_2lep else float('nan'),
  binning=[20,pi,pi],
))

plots.append(Plot( name = "l0_pt",
  texX = 'p_{T}(l_{0}) [GeV]', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.l0['pt'] if event.passing_2lep else float('nan'),
  binning=[20,0,300],
))

plots.append(Plot( name = "l1_pt",
  texX = 'p_{T}(l_{1}) [GeV]', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.l1['pt'] if event.passing_2lep else float('nan'),
  binning=[20,0,200],
))

plots.append(Plot( name = "l0_eta",
  texX = '#eta(l_{0})', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.l0['eta'] if event.passing_2lep else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = "l1_eta",
  texX = '#eta(l_{1})', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.l1['eta'] if event.passing_2lep else float('nan'),
  binning=[20,-3,3],
))

plots.append(Plot( name = "l0_phi",
  texX = '#phi(l_{0})', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.l0['phi'] if event.passing_2lep else float('nan'),
  binning=[20,pi,pi],
))

plots.append(Plot( name = "l1_phi",
  texX = '#phi(l_{1})', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.l1['phi'] if event.passing_2lep else float('nan'),
  binning=[20,pi,pi],
))

plots.append(Plot( name = 'deltaPhi_bb',
  texX = '#Delta#phi(bb)', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.deltaPhi_bb if event.passing_2lep else float('nan'),
  binning=[20,0,pi],
))

plots.append(Plot( name = 'deltaR_bb',
  texX = '#DeltaR(bb)', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.deltaR_bb if event.passing_2lep else float('nan'),
  binning=[20,0,6],
))

plots.append(Plot( name = 'deltaPhi_ll',
  texX = '#Delta#phi(ll)', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.deltaPhi_ll if event.passing_2lep else float('nan'),
  binning=[20,0,pi],
))

plots.append(Plot( name = 'deltaR_ll',
  texX = '#DeltaR(ll)', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.deltaR_ll if event.passing_2lep else float('nan'),
  binning=[20,0,4],
))

plots.append(Plot( name = '2nu_Met_pt',
  texX = 'E_{T}^{miss}(2#nu) [GeV]', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.GenMet_pt if event.passing_2lep else float('nan'),
  binning=[20,0,400],
))

plots.append(Plot( name	= '2nu_Met_phi',
  texX = '#phi(E_{T}^{miss}(2#nu))', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.GenMet_phi if event.passing_2lep else float('nan'),
  binning=[20,-pi,pi],
))

plots.append(Plot( name = 'nbjets',
  texX = 'Number of b-Jets', texY = 'Number of Events',
  attribute = lambda event, sample: len( event.trueBjets ) if event.passing_2lep else float('nan'),
  binning=[4,0,4],
))

plots.append(Plot( name = 'njets',
  texX = 'Number of Jets', texY = 'Number of Events',
  attribute = lambda event, sample: len( event.jets ) if event.passing_2lep else float('nan'),
  binning=[10,0,10],
))

plots.append(Plot( name = 'nleps',
  texX = 'Number of Leptons', texY = 'Number of Events',
  attribute = lambda event, sample: len( event.leps ) if event.passing_2lep else float('nan'),
  binning=[8,0,8],
))

plots.append(Plot( name = 'l0_pt_charge',
  texX = 'p_{T}(l_{0}) [GeV] signed with lepton charge', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.getl0chargept if event.passing_2lep else float('nan'),
  binning=[20,-200,200],
))

plots.append(Plot( name = 'l1_pt_charge',
  texX = 'p_{T}(l_{1}) [GeV] signed with lepton charge', texY = 'Number of Events / bin',
  attribute = lambda event, sample: event.getl1chargept if event.passing_2lep else float('nan'),
  binning=[20,-200,200],
))


plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

drawPlots(plots)

#logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )