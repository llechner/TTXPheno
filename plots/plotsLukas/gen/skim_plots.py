#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os, itertools
ROOT.gROOT.SetBatch(True)

from math                                import sqrt, cos, sin, pi, isnan, sinh
from RootTools.core.standard             import *
from TTXPheno.Tools.user                 import plot_directory, deltaPhi, getCollection
from TTXPheno.Tools.helpers              import deltaR, mz, sign, TransMT2, TransMT
from TTXPheno.Tools.helpers              import TransVecSum, VecSum, returnNanDict, returnNan
from TTXPheno.Tools.helpers              import createVec2, createLVec, UnitVectorT2, NUnitVectorT2
from TTXPheno.Tools.WeightInfo           import WeightInfo

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--samples',            action='store',      nargs='*',               help="Which samples?")
argParser.add_argument('--plot_directory',     action='store',      default='gen')
#argParser.add_argument('--selection',          action='store',      default='njet2p-btag1p-relIso0.12-looseLeptonVeto-mll20-met80-metSig5-dPhiJet0-dPhiJet1')
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

#
# Logger
#
import TTXPheno.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small: args.plot_directory += "_small"

# Import samples
from TTXPheno.samples.benchmarks import *

samples = map( eval, ["fwlite_ttZ_ll_LO_order3_8weights"] ) 

##
## Text on the plots
##
def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if hasData else "Robert\'s CMS Simulation"), 
      #(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

def drawPlots(plots, subDirectory=''):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, subDirectory)
    plot_directory_ = os.path.join(plot_directory_, "log") if log else os.path.join(plot_directory_, "lin")
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      plotting.draw(plot,
	    plot_directory = plot_directory_,
	    ratio = None, #{'yRange':(0.1,1.9)} if not args.noData else None,
	    logX = False, logY = log, sorting = True,
	    yRange = (0.03, "auto") if log else (0.001, "auto"),
	    scaling = {},
	    legend =  ( (0.17,0.9-0.05*sum(map(len, plot.histos))/2,1.,0.9), 2),
	    drawObjects = drawObjects( ),
      )

#
# Read variables and sequences
#

read_variables = [
    "GenMet_pt/F", "GenMet_phi/F", 
    "nGenJet/I", "GenJet[pt/F,eta/F,phi/F,matchBParton/I]", 
    "nGenLep/I", "GenLep[pt/F,eta/F,phi/F,pdgId/I,motherPdgId/I]", 
    "ntop/I", "top[pt/F,eta/F,phi/F]", 
    "Z_pt/F", "Z_eta/F", "Z_phi/F", "Z_mass/F", "Z_cosThetaStar/F", "Z_daughterPdg/I",
    "p[C/F]"
]

selection = [ 
##     ('all','1')
##    ("nlep4p", "Sum$(GenLep_pt>10&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.5)>=4&&Sum$(GenLep_pt>40&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.5)>=1"),
##    ("njet1p", "Sum$(GenJet_pt>10&&abs(GenJet_eta)<2.5)>=1"),
##    ("met40", "GenMet_phi>40"),
    ("hasZ", "Z_pt>0"),
    ("onZ", "abs(Z_mass-"+str(mZ)+")<=15"),
    ("njet3p", "Sum$(GenJet_pt>30&&abs(GenJet_eta)<2.4)>=3"),
    ("nlep3p", "Sum$(GenLep_pt>10&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.4)>=3&&Sum$(GenLep_pt>20&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.4)>=2&&Sum$(GenLep_pt>40&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.4)>=1"),
    ("nbjet1p", "Sum$(GenJet_pt>30&&GenJet_matchBParton>=1&&abs(GenJet_eta)<2.4)>=1"),
]

selectionString = "&&".join( c[1] for c in selection )
subDirectory    = '_'.join(  c[0] for c in selection )

for sample in samples:
    sample.setSelectionString( selectionString )
    sample.style = styles.lineStyle(ROOT.kBlue)


# Initialize weightstring
weight_         = None
#weight_         = xsec-factor

#WC = {'cpt':-7, 'cpQM':10, 'ctZ':2, 'ctZI':2}
#weight_ = []
#for sample in samples:
#    w = WeightInfo(sample.reweight_pkl)
#    w.set_order( 3 )
#    weightstring = '(' + w.arg_weight_string(**WC) + ')/p_C[0]'
#    print(weightstring.replace('p','event.p'))
#    weight_.append( weightstring.replace('p','event.p') )

stack = Stack(*[ [ sample ] for sample in samples] )

if args.small:
    for sample in stack.samples:
        sample.reduceFiles( to = 5 )

sequence = []

#reconstruction functions
def getblep( event, sample ):
    b = (event.bjets + event.jets)[:2]
    if TransVecSum( b[0], event.LepsFromNonZ[0], event.MetVec ).X() > TransVecSum( b[1], event.LepsFromNonZ[0], event.MetVec ).X(): return b[0]
    else: return b[1]

def getLp( event, sample ):
    WlepTransVec = TransVecSum( event.MetVec, event.LepsFromNonZ[0] )
    return ( WlepTransVec*event.LepsFromNonZ[0]['transverse_vector'] ) / ( WlepTransVec*WlepTransVec )

def getbleplepDoteZ3( event, sample ):
    return VecSum( getblep( event, sample ), event.LepsFromNonZ[0] ).Vect() * createLVec( event.Z_pt, event.Z_phi, event.Z_eta ).Vect().Unit() )

def getbleplepDoteZ2( event, sample ):
    return TransVecSum( getblep( event, sample ), event.LepsFromNonZ[0] ) * UnitVectorT2( event.Z_phi )

def gettDoteZ( event, sample ):
    tTransVec = TransVecSum( event.MetVec, event.LepsFromNonZ[0], getblep( event, sample ) )
    e_para = tTransVec * UnitVectorT2( event.ZVec_phi )
    e_norm = tTransVec * NUnitVectorT2( event.ZVec_phi )
    return e_para, e_norm

def getTransMTW( event, sample ):
    return sqrt( TransMT2( event.LepsFromNonZ[0]['pt'], event.LepsFromNonZ[0]['phi'], event.GenMet_pt, event.GenMet_phi ) )

def getTransMTt( event, sample ):
    return TransMT( event.LepsFromNonZ[0], event.MetVec, getblep( event, sample ) )

def getDeltaPhibb( event, sample ):
    b = (event.bjets + event.jets)[:2]
    return deltaPhi(b[0]['phi'], b[1]['phi'])

def getDeltaRbb( event, sample ):
    b = (event.bjets + event.jets)[:2]
    return deltaR(b[0], b[1])

def getnonZlepchargept( event, sample ):
    return -event.LepsFromNonZ[0]['pt'] * sign( event.LepsFromNonZ[0]['pdgId'] )


#sequence functions
def makeMETVec( event, sample ):
    event.MetVec = {'pt':event.GenMet_pt, 'phi':event.GenMet_phi, 'transverse_vector':createVec2(event.GenMet_pt, event.GenMet_phi)}

def makeLeps( event, sample ):
    event.leps = getCollection( event, 'GenLep', ['pt', 'eta', 'phi', 'pdgId', 'motherPdgId'], 'nGenLep' )
    event.leps = sorted( event.leps, key=lambda k: k['pt'] )
    for p in event.leps:
        p['transverse_vector'] = createVec2( p['pt'], p['phi'] )
        p['vector'] = createLVec( p['pt'], p['phi'], p['eta'] )
    event.LepsFromZ = list( filter( lambda j: j['motherPdgId'] == 23 , event.leps ) )
    if len( event.LepsFromZ ) != 2 or event.LepsFromZ[0]['pdgId'] * event.LepsFromZ[1]['pdgId'] > 0: event.LepsFromZ = []
    event.LepsFromNonZ = list( filter( lambda j: j['motherPdgId'] != 23, event.leps ) )
    
def makeJets( event, sample ):
    event.jets = getCollection( event, 'GenJet', ['pt', 'eta', 'phi', 'matchBParton'], 'nGenJet' )
    event.jets = sorted( event.jets, key=lambda k: k['pt'] )
    for p in event.jets:
        p['transverse_vector'] = createVec2( p['pt'], p['phi'] )
        p['vector'] = createLVec( p['pt'], p['phi'], p['eta'] )
    event.bjets = list( filter( lambda j: j['matchBParton'], event.jets ) )
    event.bjets = list( filter( lambda j: not j['matchBParton'], event.jets ) )


sequence.append( makeJets )
sequence.append( makeLeps )
sequence.append( makeMETVec )

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight_, selectionString = selectionString, addOverFlowBin=None)
  
plots = []

plots.append(Plot( name = "Z_pt",
  texX = 'p_{T}(Z) (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.Z_pt,
  binning=[400/20,0,400],
))

plots.append(Plot( name = "Z_mass",
  texX = 'm(ll) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.Z_mass,
  binning=[50,50,130],
))

plots.append(Plot( name = "Z_cosThetaStar",
  texX = 'cos(#theta*)', texY = 'Number of Events',
  attribute = lambda event, sample: event.Z_cosThetaStar,
  binning=[50,-1.2,1.2],
))

plots.append(Plot( name = "b0_pt",
  texX = 'p_{T}(b0) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: (event.bjets + event.jets)[0]['pt'],
  binning=[400/20,0,400],
))

plots.append(Plot( name = "b1_pt",
  texX = 'p_{T}(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: (event.bjets + event.jets)[1]['pt'],
  binning=[400/20,0,400],
))

plots.append(Plot( name = "b0_eta",
  texX = '#eta(b0) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: (event.bjets + event.jets)[0]['eta'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = "b1_eta",
  texX = '#eta(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: (event.bjets + event.jets)[1]['eta'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = "b0_phi",
  texX = '#Phi(b0) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: (event.bjets + event.jets)[0]['phi'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = "b1_phi",
  texX = '#Phi(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: (event.bjets + event.jets)[1]['phi'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaPhi_bb',
  texX = '#Delta#phi(bb)', texY = 'Number of Events',
  attribute = lambda event, sample: getDeltaPhibb( event, sample ),
  binning=[32,-0.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaR_bb',
  texX = '#DeltaR(bb)', texY = 'Number of Events',
  attribute = lambda event, sample: getDeltaRbb( event, sample ),
  binning=[32,-0.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaPhi_ll',
  texX = '#Delta#phi(ll)', texY = 'Number of Events',
  attribute = lambda event, sample: deltaPhi( event.LepsFromZ[0]['phi'], event.LepsFromZ[1]['phi'] ),
  binning=[32,-0.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaR_ll',
  texX = '#DeltaR(ll)', texY = 'Number of Events',
  attribute = lambda event, sample: deltaR( event.LepsFromZ[0], event.LepsFromZ[1] ),
  binning=[32,-0.2*pi,1.2*pi],
))

plots.append(Plot( name = 'Met_pt',
  texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.GenMet_pt,
  binning=[400/20,0,400],
))

plots.append(Plot( name	= 'Met_phi',
  texX = '#phi(E_{T}^{miss})', texY = 'Number of Events / 20 GeV',
  attribute = lambda event, sample: event.GenMet_phi,
  binning=[10,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = 'nbjets',
  texX = 'number of gen bjets', texY = 'Number of Events',
  attribute = lambda event, sample: len( event.bjets ),
  binning=[6,0,6],
))

plots.append(Plot( name = 'njets',
  texX = 'number of gen jets', texY = 'Number of Events',
  attribute = lambda event, sample: len( event.jets ),
  binning=[10,0,10],
))

plots.append(Plot( name = 'nleps',
  texX = 'number of gen leps', texY = 'Number of Events',
  attribute = lambda event, sample: len( event.leps ),
  binning=[10,0,10],
))

plots.append(Plot( name = 'Lp',
  texX = 'L_{p}', texY = 'Number of Events',
  attribute = lambda event, sample: getLp( event, sample ),
  binning=[50,-1,1],
))

plots.append(Plot( name = 'W_dot_eZ2',
  texX = '(b_{lep}+l) * e(Z) (2D)', texY = 'Number of Events',
  attribute = lambda event, sample: getbleplepDoteZ2( event, sample ),
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'W_dot_eZ3',
  texX = '(b_{lep}+l) * e(Z) (3D)', texY = 'Number of Events',
  attribute = lambda event, sample: getbleplepDoteZ3( event, sample ),
  binning=[400/20,0,400],
))

plots.append(Plot( name = 't_dot_eZ',
  texX = 't_{lep}^{rec} * e(Z)', texY = 'Number of Events',
  attribute = lambda event, sample: gettDoteZ( event, sample )[0],
  binning=[400/20,0,400],
))

plots.append(Plot( name = 't_dot_enZ',
  texX = 't_{lep}^{rec} * e_{n}(Z)', texY = 'Number of Events',
  attribute = lambda event, sample: gettDoteZ( event, sample )[1],
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'l_pt_charge',
  texX = 'p_{T}(l) (GeV) * sign(q(l))', texY = 'Number of Events',
  attribute = lambda event, sample: getnonZlepchargept( event, sample ),
  binning=[50,-100,100],
))

plots.append(Plot( name = 'mT_W',
  texX = 'm_{T}(W) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getTransMTW( event, sample ),
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'mT_t',
  texX = 'm_{T}(t) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getTransMTt( event, sample ),
  binning=[400/20,0,400],
))



plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = 100 if args.small else -1)

drawPlots(plots, subDirectory = subDirectory)

#logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
