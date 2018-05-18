#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)

from math                                import sqrt, cos, sin, pi, isnan, sinh
from RootTools.core.standard             import *
from TTXPheno.Tools.user                 import plot_directory
from TTXPheno.Tools.helpers              import deltaPhi
from TTXPheno.Tools.helpers              import getCollection
from TTXPheno.Tools.helpers              import deltaR
from TTXPheno.Tools.helpers              import mZ
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
      (0.15, 0.95, 'CMS Preliminary' if hasData else 'CMS Simulation'), 
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
    ("nlep3p", "Sum$(GenLep_pt>10&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.5)>=3&&Sum$(GenLep_pt>20&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.5)>=2&&Sum$(GenLep_pt>40&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.5)>=1"),
    ("onZ", "abs(Z_mass-"+str(mZ)+")<=15"),
    ("hasZ", "Z_pt>0"),
    ("lepSelTTZ", "Sum$(GenLep_pt>40)>=1 && Sum$(GenLep_pt>20)>=2 && Sum$(GenLep_pt>10)>=3"),
    ("njet3p", "Sum$(GenJet_pt>30&&abs(GenJet_eta)<2.5)>=3"),
    ("nbjet1p", "Sum$(GenJet_pt>30&&GenJet_matchBParton>=1&&abs(GenJet_eta)<2.5)>=1"),
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

#additional functions
def sign( x ):
    return x and (1, -1)[x < 0]

def TransMT( Vec2 ):
    return sqrt( Vec2.Mod2() )

def TransVecSum( *args ):
    return sum( [p['transverse_vector'] for p in args], ROOT.TVector2() )

def VecSum( *args ):
    return sum( [p['vector'] for p in args], ROOT.TLorentzVector() )

def returnNanDict():
    return {'pt':float('nan'), 'phi':float('nan'), 'pdgId':float('nan'), 'eta':float('nan'), 'motherPdgId':float('nan'), 'matchBParton':float('nan'), 'transverse_vector':float('nan'), 'vector':float('nan')}

def returnNan():
    return float('nan')

def createVec2( pt, phi ):
    return ROOT.TVector2( pt*cos(phi), pt*sin(phi) )

def createLVec( pt, phi, eta ):
    return ROOT.TLorentzVector( pt*cos(phi), pt*sin(phi), pt*sinh(eta), 0 )

def UnitVectorT3( Vec3 ):
    return Vec3.Unit()

def UnitVectorT2( phi ):
    return ROOT.TVector2( cos(phi), sin(phi) )

def NUnitVectorT2( phi ):
    return ROOT.TVector2( -sin(phi), cos(phi) )

#reconstruction functions
def getb( event, sample, index ):
    if event.nbjets <= index: return returnNanDict()
    return sorted( event.bjets, key=lambda k: k['pt'] )[index]

def getBorJet( event, sample, index ):
    if event.nbjets <= index: return getjet( event, sample, index - event.nbjets )
    else: return getb( event, sample, index )

def getjet( event, sample, index ):
    if event.nNonbjets <= index: return returnNanDict()
    return sorted( event.Nonbjets, key=lambda k: k['pt'] )[index]

def getNonZLep( event, sample, index ):
    if event.nLepsFromNonZ <= index: return returnNanDict()
    return sorted( event.LepsFromNonZ, key=lambda k: k['pt'] )[index]

def getblep( event, sample ):
    b1 = getBorJet( event, sample, 0 )
    b2 = getBorJet( event, sample, 1 )
    lep = getNonZLep( event, sample, 0 ) #really the right lepton?
    if TransVecSum( b1, lep, event.MetVec ).X() > TransVecSum( b2, lep, event.MetVec ).X(): return b1
    else: return b2

def getWlepTransVec( event, sample ):
    lep = getNonZLep( event, sample, 0 ) #really the right lepton?
    return TransVecSum( event.MetVec, lep )

def gettlepTransVec( event, sample ):
    blepTransVec = getblep( event, sample )['transverse_vector']
    WlepTransVec = getWlepTransVec( event, sample )
    return blepTransVec + WlepTransVec

def getbleplepLVec( event, sample ):
    blep = getblep( event, sample )
    lep = getNonZLep( event, sample, 0 ) #really the right lepton?
    return VecSum( blep, lep )

def getbleplepTransVec( event, sample ):
    blep = getblep( event, sample )
    lep = getNonZLep( event, sample, 0 ) #really the right lepton?
    return TransVecSum( blep, lep )

def getLp( event, sample ):
    lepTransVec = getNonZLep( event, sample, 0 )['transverse_vector'] #really the right lepton?
    WlepTransVec = getWlepTransVec( event, sample )
    return ( WlepTransVec*lepTransVec ) / ( WlepTransVec*WlepTransVec )

def getbleplepDoteZ3( event, sample ):
    return getbleplepLVec( event, sample ).Vect()*UnitVectorT3( event.ZVec['vector'].Vect() )

def getbleplepDoteZ2( event, sample ):
    return getbleplepTransVec( event, sample )*UnitVectorT2( event.ZVec['phi'] )

def gettDoteZ2( event, sample ):
    return gettlepTransVec( event, sample )*UnitVectorT2( event.ZVec['phi'] )

def gettDotneZ2( event, sample ):
    return gettlepTransVec( event, sample )*NUnitVectorT2( event.ZVec['phi'] )

def getTransMTW( event, sample ):
    return TransMT( getWlepTransVec( event, sample ) )

def getTransMTt( event, sample ):
    return TransMT( gettlepTransVec( event, sample ) )

#def getnonZlepchargept( event, sample ):
#    lep = getNonZLep( event, sample, 0 ) #really the right lepton?
#    return -lep['pt'] * sign( lep['pdgId'] )


#sequence functions
def makeMETVec( event, sample ):
    event.MetVec = {'pt':event.GenMet_pt, 'phi':event.GenMet_phi, 'transverse_vector':createVec2(event.GenMet_pt, event.GenMet_phi)}

def makeZVec( event, sample ):
    event.ZVec = {'pt':event.Z_pt, 'phi':event.Z_phi, 'eta':event.Z_eta, 'transverse_vector':createVec2(event.Z_pt, event.Z_phi), 'vector':createLVec(event.Z_pt, event.Z_phi, event.Z_eta)}

def makeLeps( event, sample ):
    event.leps = getCollection( event, 'GenLep', ['pt', 'eta', 'phi', 'pdgId', 'motherPdgId'], 'nGenLep' )
    for p in event.leps:
        p['transverse_vector'] = createVec2( p['pt'], p['phi'] )
        p['vector'] = createLVec( p['pt'], p['phi'], p['eta'] )
    event.nleps = len( event.leps )
    event.LepsFromZ = filter( lambda j: j['motherPdgId'] == 23, event.leps )
    event.nLepsFromZ = len( event.LepsFromZ )
    event.LepsFromNonZ = filter( lambda j: j['motherPdgId'] != 23, event.leps )
    event.nLepsFromNonZ = len( event.LepsFromNonZ )
    
def makeJets( event, sample ):
    event.jets = getCollection( event, 'GenJet', ['pt', 'eta', 'phi', 'matchBParton'], 'nGenJet' )
    for p in event.jets:
        p['transverse_vector'] = createVec2( p['pt'], p['phi'] )
        p['vector'] = createLVec( p['pt'], p['phi'], p['eta'] )
    event.njets = len( event.jets )
    event.bjets = filter( lambda j: j['matchBParton'], event.jets )
    event.nbjets = len( event.bjets )
    event.Nonbjets = filter( lambda j: not j['matchBParton'], event.jets )
    event.nNonbjets = len( event.Nonbjets )

def makeDeltaPhill( event, sample ):
    if event.nLepsFromZ < 2 or event.LepsFromZ[0]['pdgId'] * event.LepsFromZ[1]['pdgId'] > 0: event.dPhi_ll = returnNan()
    else: event.dPhi_ll = deltaPhi( event.LepsFromZ[0]['phi'], event.LepsFromZ[1]['phi'] )

def makeDeltaRll( event, sample ):
    if event.nLepsFromZ < 2 or event.LepsFromZ[0]['pdgId'] * event.LepsFromZ[1]['pdgId'] > 0: event.dR_ll = returnNan()
    else: event.dR_ll = deltaR( event.LepsFromZ[0], event.LepsFromZ[1] )

def makeDeltaPhibb( event, sample ):
    b1 = getBorJet( event, sample, 0 )
    b2 = getBorJet( event, sample, 1 )
    event.dPhi_bb = deltaPhi(b1['phi'], b2['phi'])

def makeDeltaRbb( event, sample ):
    b1 = getBorJet( event, sample, 0 )
    b2 = getBorJet( event, sample, 1 )
    event.dR_bb = deltaR(b1, b2)

sequence.append( makeJets )
sequence.append( makeLeps )
sequence.append( makeDeltaPhill )
sequence.append( makeDeltaRll )
sequence.append( makeDeltaPhibb )
sequence.append( makeDeltaRbb )
sequence.append( makeMETVec )
sequence.append( makeZVec )

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

plots.append(Plot( name = "b1_pt",
  texX = 'p_{T}(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getBorJet( event, sample, 0 )['pt'],
  binning=[400/20,0,400],
))

plots.append(Plot( name = "b2_pt",
  texX = 'p_{T}(b2) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getBorJet( event, sample, 1 )['pt'],
  binning=[400/20,0,400],
))

plots.append(Plot( name = "b1_eta",
  texX = '#eta(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getBorJet( event, sample, 0 )['eta'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = "b2_eta",
  texX = '#eta(b2) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getBorJet( event, sample, 1 )['eta'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = "b1_phi",
  texX = '#Phi(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getBorJet( event, sample, 0 )['phi'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = "b2_phi",
  texX = '#Phi(b2) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getBorJet( event, sample, 1 )['phi'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaPhi_bb',
  texX = '#Delta#phi(bb)', texY = 'Number of Events',
  attribute = lambda event, sample: event.dPhi_bb,
  binning=[32,-0.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaR_bb',
  texX = '#DeltaR(bb)', texY = 'Number of Events',
  attribute = lambda event, sample: event.dR_bb,
  binning=[32,-0.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaPhi_ll',
  texX = '#Delta#phi(ll)', texY = 'Number of Events',
  attribute = lambda event, sample: event.dPhi_ll,
  binning=[32,-0.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaR_ll',
  texX = '#DeltaR(ll)', texY = 'Number of Events',
  attribute = lambda event, sample: event.dR_ll,
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
  attribute = lambda event, sample: event.nbjets,
  binning=[6,0,6],
))

plots.append(Plot( name = 'njets',
  texX = 'number of gen jets', texY = 'Number of Events',
  attribute = lambda event, sample: event.njets,
  binning=[10,0,10],
))

plots.append(Plot( name = 'nleps',
  texX = 'number of gen leps', texY = 'Number of Events',
  attribute = lambda event, sample: event.nleps,
  binning=[10,0,10],
))

plots.append(Plot( name = 'Lp',
  texX = 'L_{p}', texY = 'Number of Events',
  attribute = lambda event, sample: getLp( event, sample ),
  binning=[50,-1,1], #??
))

plots.append(Plot( name = 'bleplepDoteZ2',
  texX = '(b_{lep}+l) * e(Z) (2D)', texY = 'Number of Events',
  attribute = lambda event, sample: getbleplepDoteZ2( event, sample ),
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'bleplepDoteZ3',
  texX = '(b_{lep}+l) * e(Z) (3D)', texY = 'Number of Events',
  attribute = lambda event, sample: getbleplepDoteZ3( event, sample ),
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'tDoteZ',
  texX = 't_{lep}^{rec} * e(Z)', texY = 'Number of Events',
  attribute = lambda event, sample: gettDoteZ2( event, sample ),
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'tDotenZ',
  texX = 't_{lep}^{rec} * e_{n}(Z)', texY = 'Number of Events',
  attribute = lambda event, sample: gettDotneZ2( event, sample ),
  binning=[400/20,0,400],
))

#plots.append(Plot( name = 'leppt_charge',
#  texX = 'p_{T}(l) (GeV) * sign(q(l))', texY = 'Number of Events',
#  attribute = lambda event, sample: getnonZlepchargept( event, sample ),
#  binning=[50,-100,100],
#))

plots.append(Plot( name = 'MissingMass_W',
  texX = 'm_{T}(W) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getTransMTW( event, sample ),
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'MissingMass_t',
  texX = 'm_{T}(t) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: getTransMTt( event, sample ),
  binning=[400/20,0,400],
))



plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = 100 if args.small else -1)

drawPlots(plots, subDirectory = subDirectory)

#logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
