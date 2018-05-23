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
from TTXPheno.Tools.user                 import plot_directory
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR, mZ
from TTXPheno.Tools.WeightInfo           import WeightInfo

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--samples',            action='store',      nargs='*',               help="Which samples?")
argParser.add_argument('--plot_directory',     action='store',      default='gen')
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

preselection = [ 
##     ('all','1')
##    ("nlep4p", "Sum$(GenLep_pt>10&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.5)>=4&&Sum$(GenLep_pt>40&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.5)>=1"),
##    ("njet1p", "Sum$(GenJet_pt>10&&abs(GenJet_eta)<2.5)>=1"),
##    ("met40", "GenMet_phi>40"),
    ("hasZ", "Z_pt>0"),
    ("onZ", "abs(Z_mass-%f)<=15"%mZ),
    ("njet3p", "Sum$(GenJet_pt>30&&abs(GenJet_eta)<2.4)>=3"),
    ("nlep3p", "Sum$(GenLep_pt>10&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.4)>=3&&Sum$(GenLep_pt>20&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.4)>=2&&Sum$(GenLep_pt>40&&(abs(GenLep_pdgId)==11||abs(GenLep_pdgId)==13)&&abs(GenLep_eta)<2.4)>=1"),
    ("nbjet1p", "Sum$(GenJet_pt>30&&GenJet_matchBParton>=1&&abs(GenJet_eta)<2.4)>=1"),
]

selectionString = "&&".join( c[1] for c in preselection )
subDirectory    = '_'.join(  c[0] for c in preselection )

for sample in samples:
    sample.setSelectionString( selectionString )
    sample.style = styles.lineStyle(ROOT.kBlue)


# Initialize weightstring

#WC = {'cpt':-7, 'cpQM':10, 'ctZ':2, 'ctZI':2}
#for sample in samples:
#    w = WeightInfo(sample.reweight_pkl)
#    w.set_order( 3 )
#    weightstring = '(' + w.arg_weight_string(**WC) + ')/p_C[0]'
#    print(weightstring.replace('p','event.p'))
#    weight_.append( weightstring.replace('p','event.p') )

stack = Stack(*[ [ sample ] for sample in samples] )

if args.small:
    for sample in stack.samples:
        sample.reduceFiles( to = 1 )

# Helpers
def addTransverseVector( p_dict ):
    ''' add a transverse vector for further calculations
    '''
    p_dict['vec2D'] = ROOT.TVector2( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']) )

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vec4D'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), 0 )

def NanJet():
    return {'pt':float('nan'), 'phi':float('nan'), 'pdgId':float('nan'), 'eta':float('nan'), 'motherPdgId':float('nan'), 'matchBParton':float('nan'), 'transverse_vector':float('nan'), 'vector':float('nan')}

def UnitVectorT2( phi ):
    return ROOT.TVector2( cos(phi), sin(phi) )

def isGoodJet( j ):
    ''' jet object selection
    '''
    return j['pt']>30 and abs(j['eta'])<2.4

def isGoodLepton( l ):
    ''' lepton object selection
    '''
    return l['pt']>10 and abs(l['eta'])<2.5

def MTSquared( p1, p2 ):
    ''' compute MT from 2 particles
    '''
    return 2*p1['pt']*p2['pt']*(1-cos(p1['phi']-p2['phi']) )

def MSquared( p1, p2 ):
    ''' compute MassSquared from 2 particles
    '''
    return 2*p1['pt']*p2['pt']*(cosh(p1['eta']-p2['eta'])-cos(p1['phi']-p2['phi']) )

#def vecSum( *args ):
#    if type(args[0])==type(ROOT.TLorentzVector()):
#        return sum( [p['vec4D'] for p in args], ROOT.TLorentzVector() )
#    else:
#        return sum( [p['vec2D'] for p in args], ROOT.TVector2() )

#sequence functions
sequence = []

def makeJets( event, sample ):
    ''' Add a list of filtered all jets to the event
    '''

    # Retrieve & filter
    event.jets = list( 
            filter( lambda j:isGoodJet( j ), \
            getCollection( event, 'GenJet', ['pt', 'eta', 'phi', 'matchBParton'], 'nGenJet' ) )
            )
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
    event.bj0, event.bj1 = ( event.trueBjets + event.trueNonBjets + [NanJet(), NanJet()])[:2] 
    
sequence.append( makeJets )

def makeMET( event, sample ):
    ''' Make a MET vector to facilitate further calculations
    '''
    event.MET = {'pt':event.GenMet_pt, 'phi':event.GenMet_phi}
    addTransverseVector( event.MET )

sequence.append( makeMET )

def makeLeps( event, sample ):
    event.leps = list( filter( 
        lambda l: isGoodLepton( l ),
        getCollection( event, 'GenLep', ['pt', 'eta', 'phi', 'pdgId', 'motherPdgId'], 'nGenLep' )
        ))

    # Add extra vectors
    for p in event.leps:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Sort
    event.leps = sorted( event.leps, key=lambda k: -k['pt'] )

    # Cross-cleaning: remove leptons that overlap with a jet within 0.4
    event.leps = list(filter( lambda l: min( [ deltaR(l, j) for j in event.jets ] + [999] ) > 0.4 , event.leps ))

    # find leptons from Z
    event.lepsFromZ = list( filter( lambda j: j['motherPdgId'] == 23 , event.leps ) )
    event.foundZ    = len( event.lepsFromZ )==2 and event.lepsFromZ[0]['pdgId'] * event.lepsFromZ[1]['pdgId'] < 0
    event.Z_deltaPhi_ll = deltaPhi( *event.lepsFromZ) if event.foundZ else float('nan')
    event.Z_deltaR_ll   = deltaR( *event.lepsFromZ) if event.foundZ else float('nan')
 
    # convinience
    event.Z_unitVec2D = UnitVectorT2( event.Z_phi )
    event.Z_vec4D     = ROOT.TLorentzVector()
    event.Z_vec4D.SetPtEtaPhiM( event.Z_pt, event.Z_eta, event.Z_phi, event.Z_mass )
    event.Z_unitVec3D = event.Z_vec4D.Vect()
    #event.Z_unitVec3D /= event.Z_unitVec3D.Mag() 

    # find leptons that are NOT from Z 
    event.leptonsNotFromZ = [lepton for lepton in event.leps if lepton not in event.lepsFromZ] 

    # We may loose some events by cross-cleaning or by thresholds. Make a bool for that.
    event.passing_3lep    = event.foundZ and len(event.leptonsNotFromZ)==1

    print event.foundZ, event.leptonsNotFromZ

sequence.append( makeLeps )

def makeObservables( event, sample):
    ''' Compute all relevant observables
    '''

    # Only continue if we pass the 3-lep selection at least
    if not event.passing_3lep: return

    # Resolve pairing ambiguity by maximizing resulting top-lep pt  
    if vecSum( event.bj0, event.leptonsNotFromZ[0], event.MET ).Mag2() > vecSum( event.bj1, event.leptonsNotFromZ[0], event.MET ).Mag2(): 
        event.b_lep = bj0, event.b_had = bj1
    else:
        event.b_lep = bj1, event.b_had = bj0

    # Make leptonic W
    event.Wlep_vec2D = vecSum( event.MET, event.leptonsNotFromZ[0] )
    # Lp
    event.W_Lp = (Wlep_vec2D*event.leptonsNotFromZ[0]['vec2D'] ) / ( event.Wlep_vec2D*event.Wlep_vec2D )
    # classic MT
    event.W_MT = sqrt( MTSquared( event.MET, event.leptonsNotFromZ[0] ) )
        
    # blep+lep subsystem
    event.bleplep_vec2D = vecSum( b_lep['vec2D'], event.leptonsNotFromZ[0]['vec2D'] )
    event.bleplep_vec4D = vecSum( b_lep['vec4D'], event.leptonsNotFromZ[0]['vec4D'] )

    # transverse mass of top
    event.t_MT = sqrt( MTSquared( event.leptonsNotFromZ[0], event.MET) + MSquared( event.leptonsNotFromZ[0], event.b_lep ) + MTSquared( event.MET, event.b_lep ) )
    event.t_vec2D = vecSum( event.leptonsNotFromZ[0]['vec2D'], event.MET['vec2D'], event.b_lep['vec2D'] )

    # double b kinematic
    event.deltaPhi_bb = deltaPhi(bj0['phi'], bj1['phi'])
    event.deltaR_bb = deltaR(bj0, bj1)

    # Well... signed lepton pt.
    event.getnonZlepchargept = leptonsNotFromZ[0]['pt'] if leptonsNotFromZ[0]['pdgId']>0 else -leptonsNotFromZ[0]['pt'] 

sequence.append( makeObservables )

# Weight <- Here we remove events where leptons fail the analysis selection despite passing the preselection
weight_ = None#lambda event, sample: event.passing_3lep
    
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
  attribute = lambda event, sample: event.bj0['pt'],
  binning=[400/20,0,400],
))

plots.append(Plot( name = "b1_pt",
  texX = 'p_{T}(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bj1['pt'],
  binning=[400/20,0,400],
))

plots.append(Plot( name = "b0_eta",
  texX = '#eta(b0) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bj0['eta'],
  binning=[60,-3,3],
))

plots.append(Plot( name = "b1_eta",
  texX = '#eta(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bj1['eta'],
  binning=[60,-3,3],
))

plots.append(Plot( name = "b0_phi",
  texX = '#Phi(b0) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bj0['phi'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = "b1_phi",
  texX = '#Phi(b1) (GeV)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bj1['phi'],
  binning=[20,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaPhi_bb',
  texX = '#Delta#phi(bb)', texY = 'Number of Events',
  attribute = lambda event, sample: event.deltaPhi_bb,
  binning=[32,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = 'deltaR_bb',
  texX = '#DeltaR(bb)', texY = 'Number of Events',
  attribute = lambda event, sample: event.deltaR_bb,
  binning=[32,-6,6],
))

plots.append(Plot( name = 'Z_deltaPhi_ll',
  texX = '#Delta#phi(ll) from Z', texY = 'Number of Events',
  attribute = lambda event, sample: event.Z_deltaPhi_ll,
  binning=[32,-0.2*pi,1.2*pi],
))

plots.append(Plot( name = 'Z_deltaR_ll',
  texX = '#DeltaR(ll)from Z', texY = 'Number of Events',
  attribute = lambda event, sample: event.Z_deltaR_ll,
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
  binning=[50,-1.2*pi,1.2*pi],
))

plots.append(Plot( name = 'gen_nbjets',
  texX = 'number of gen bjets', texY = 'Number of Events',
  attribute = lambda event, sample: len( event.trueBjets ),
  binning=[10,0,10],
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

plots.append(Plot( name = 'W_Lp',
  texX = 'L_{p} from W', texY = 'Number of Events',
  attribute = lambda event, sample: event.W_Lp,
  binning=[50,-1,1],
))

plots.append(Plot( name = 'W_dot_eZ2',
  texX = '(b_{lep}+l) #cdot e(Z) (2D)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bleplep_vec2D*event.Z_unitVec2D,
  binning=[400/20,0,400],
))

plots.append(Plot( name = 'W_dot_eZ3',
  texX = '(b_{lep}+l) * e(Z) (3D)', texY = 'Number of Events',
  attribute = lambda event, sample: event.bleplep_vec4D.Vect()*event.Z_unitVec3D,
  binning=[400/20,0,400],
))

plots.append(Plot( name = 't_dot_eZ',
  texX = 't_{lep}^{rec} * e(Z)', texY = 'Number of Events',
  attribute = lambda event, sample: event.t_vec2D*event.Z_unitVec2D,
  binning=[400/20,0,400],
))

#plots.append(Plot( name = 't_dot_enZ',
#  texX = 't_{lep}^{rec} * e_{n}(Z)', texY = 'Number of Events',
#  attribute = lambda event, sample: gettDoteZ( event, sample )[1],
#  binning=[400/20,0,400],
#))
#
#plots.append(Plot( name = 'l_pt_charge',
#  texX = 'p_{T}(l) (GeV) * sign(q(l))', texY = 'Number of Events',
#  attribute = lambda event, sample: getnonZlepchargept( event, sample ),
#  binning=[50,-100,100],
#))
#
#plots.append(Plot( name = 'mT_W',
#  texX = 'm_{T}(W) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: getTransMTW( event, sample ),
#  binning=[400/20,0,400],
#))
#
#plots.append(Plot( name = 'mT_t',
#  texX = 'm_{T}(t) (GeV)', texY = 'Number of Events',
#  attribute = lambda event, sample: getTransMTt( event, sample ),
#  binning=[400/20,0,400],
#))



plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

drawPlots(plots, subDirectory = subDirectory)

#logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
