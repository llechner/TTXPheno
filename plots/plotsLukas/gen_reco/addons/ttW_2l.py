#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
ROOT.gROOT.SetBatch(True)
from math                                import pi, sqrt

# RootTools
from RootTools.core.standard             import *

# TTXPheno
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR, nanJet, nanLepton, getObjDict

# Import samples
from TTXPheno.samples.benchmarks         import *

# Import helpers
from plot_helpers                        import *


def getVariableList( level ):
   
    # List of variables where gen is replaced by reco for reco
    read_variables_gen = [
#        "ref_lumiweight1fb/F",
        "lumiweight1fb/F",
        "genMet_pt/F", "genMet_phi/F",
    
        "ngenJet/I", #"genJet[pt/F,eta/F,phi/F]",
        "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I]",
    
        "genW_pt/F", "genW_eta/F", "genW_phi/F", "genW_mass/F",
    
        "genBj0_pt/F", "genBj0_phi/F", "genBj0_eta/F",
        "genBj1_pt/F", "genBj1_phi/F", "genBj1_eta/F",
    ]
     
    # List of variables where genLep is replaced by reco for reco
    read_variables_genLep = [
    ]
    
    if level == 'reco':
        read_variables_gen    = [ variable.replace('gen', 'reco') for variable in read_variables_gen ]
        read_variables_genLep = [ variable.replace('genLep', 'reco') for variable in read_variables_genLep ]
        read_variables_gen.append("recoJet[pt/F,eta/F,phi/F,bTag/F]")
    else:
        read_variables_gen.append("genJet[pt/F,eta/F,phi/F,matchBParton/F]")

    read_variables = read_variables_gen + read_variables_genLep
    read_variables = list( set( read_variables ) ) # remove double entries
#    read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

    return read_variables


def makeJets( event, sample, level ):
    ''' Add a list of filtered jets to the event (full list is required for lepton cross cleaning)
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    # load jets
    btag = 'bTag' if level == 'reco' else 'matchBParton'
    event.jets = getCollection( event, '%sJet'%preTag, ['pt', 'eta', 'phi', btag ], 'n%sJet'%preTag )
    event.bjets = list( filter( lambda j: j[btag], event.jets ) )

    # get (second) hardest bjets
    event.bj0 = {'pt':getattr( event, '%sBj0_pt'%preTag ), 'phi':getattr( event, '%sBj0_phi'%preTag ), 'eta':getattr( event, '%sBj0_eta'%preTag )}
    event.bj1 = {'pt':getattr( event, '%sBj1_pt'%preTag ), 'phi':getattr( event, '%sBj1_phi'%preTag ), 'eta':getattr( event, '%sBj1_eta'%preTag )}

    # Add extra vectors
    for p in [event.bj0, event.bj1]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoJet       as isGoodJet
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenJet      as isGoodJet

    # selection checks
    event.foundBj0 =  isGoodJet( event.bj0 )

    # choose your selection on b-jets
    event.passing_bjets = event.foundBj0


def makeMET( event, sample, level ):
    ''' Make a MET vector to facilitate further calculations
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    event.MET = {'pt':getattr(event, '%sMet_pt'%preTag), 'phi':getattr(event, '%sMet_phi'%preTag)}
    addTransverseVector( event.MET )

def makeW( event, sample, level ):
    ''' Make a W vector to facilitate further calculations (either recoZ, genLepZ or genZ)
    '''
    event.W_unitVec2D = UnitVectorT2( getattr( event, '%sW_phi'%level     ) )
    event.W_vec4D     = ROOT.TLorentzVector()
    event.W_vec4D.SetPtEtaPhiM( getattr( event, '%sW_pt'%level ), getattr( event, '%sW_eta'%level ), getattr( event, '%sW_phi'%level ), getattr( event, '%sW_mass'%level     ) )
    event.W_unitVec3D = event.W_vec4D.Vect().Unit()


def makeLeps( event, sample, level ):
    ''' Add important leptons (no full list of leptons is required for now)
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    # Define W leptons
    event.W_l0 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], 0 ) 
    event.W_l1 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], 1 )

    # Add extra vectors
    for p in [ event.W_l0, event.W_l1]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoLepton    as isGoodLepton
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenLepton   as isGoodLepton

    # We may loose some events by cross-cleaning or by thresholds.
    event.foundWl0     = isGoodLepton( event.W_l0 )
    event.foundWl1     = isGoodLepton( event.W_l1 )
    event.found2lep    = event.W_l0['pdgId'] * event.W_l1['pdgId'] > 0

    # choose your selection on leptons
    event.passing_leptons = event.found2lep and event.foundWl0 and event.foundWl1


def makeObservables( event, sample, level):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.bbdPhi = deltaPhi( event.bj0['phi'], event.bj1['phi'] )
    event.bbdR   = deltaR( event.bj0, event.bj1 )

    # double l kinematic
    event.lldPhi = deltaPhi( event.W_l0['phi'], event.W_l1['phi'] )
    event.lldR   = deltaR( event.W_l0, event.W_l1 )

    # signed lepton pt
    event.lep0chargept = event.W_l0['pt'] if event.W_l0['pdgId']>0 else -event.W_l0['pt']
    event.lep1chargept = event.W_l1['pt'] if event.W_l1['pdgId']>0 else -event.W_l1['pt']

    # choose your final selection
    event.passing_checks = event.passing_leptons and event.passing_bjets


def getSequenceList( level, sameFlavor ):
    ''' sequence functions
    '''
    sequence = []

    sequence.append( lambda event, sample: makeJets( event, sample, level ) )
    sequence.append( lambda event, sample: makeMET( event, sample, level ) )
    sequence.append( lambda event, sample: makeW( event, sample, level ) )
    sequence.append( lambda event, sample: makeLeps( event, sample, level ) )
    sequence.append( lambda event, sample: makeObservables( event, sample, level ) )

    return sequence
    

def getPlotList( scaleLumi, level ):

    tag = 'reco' if level == 'reco' else 'genLep'
    preTag = 'reco' if level == 'reco' else 'gen'

    if scaleLumi: y_label = 'norm. diff. xsec'
    else:         y_label = 'diff. x-sec'

    plots = []
    
#    plots.append(Plot( name = 'nbjets',
#      texX = 'Number of b-Jets', texY = 'Number of Events',
#      attribute = lambda event, sample: getattr( event, 'n%sJet'%preTag ) if event.passing_checks else float('nan'),
#      binning=[4,0,4],
#    ))
    
    plots.append(Plot( name = "W_pt",
      texX = 'p_{T}(W) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, "%sW_pt"%level ) if event.passing_checks else float('nan'),
      binning=[20,0,500],
    ))
    
    plots.append(Plot( name = 'W_phi',
      texX = '#phi(W) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, "%sW_phi"%level ) if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    
    plots.append(Plot( name = 'W_eta',
      texX = '#eta(W) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, "%sW_eta"%level ) if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = "W_mass",
      texX = 'm(W) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, "%sW_mass"%level ) if event.passing_checks else float('nan'),
      binning=[20,60,100],
    ))
    
    plots.append(Plot( name = "b0_pt",
      texX = 'p_{T}(b_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.bj0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name = "b0_phi",
      texX = '#phi(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    
    plots.append(Plot( name = "b0_eta",
      texX = '#eta(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = "b1_pt",
      texX = 'p_{T}(b_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.bj1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name = "b1_phi",
      texX = '#phi(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    
    plots.append(Plot( name = "b1_eta",
      texX = '#eta(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = 'deltaPhi_bb',
      texX = '#Delta#phi(bb)', texY = y_label,
      attribute = lambda event, sample: event.bbdPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    
    plots.append(Plot( name = 'deltaR_bb',
      texX = '#DeltaR(bb)', texY = y_label,
      attribute = lambda event, sample: event.bbdR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    
    plots.append(Plot( name = 'deltaPhi_ll',
      texX = '#Delta#phi(ll)', texY = y_label,
      attribute = lambda event, sample: event.lldPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    
    plots.append(Plot( name = 'deltaR_ll',
      texX = '#DeltaR(ll)', texY = y_label,
      attribute = lambda event, sample: event.lldR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    
    plots.append(Plot( name = '2nu_Met_pt',
      texX = 'E_{T}^{miss}(2#nu) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.MET['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name	= '2nu_Met_phi',
      texX = '#phi(E_{T}^{miss}(2#nu))', texY = y_label,
      attribute = lambda event, sample: event.MET['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    
    plots.append(Plot( name = 'l0_pt',
      texX = 'p_{T}(l_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.W_l0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ))
    
    plots.append(Plot( name = 'l0_phi',
      texX = '#phi(l_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.W_l0['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    
    plots.append(Plot( name = 'l0_eta',
      texX = '#eta(l_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.W_l0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = 'l1_pt',
      texX = 'p_{T}(l_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.W_l1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ))
    
    plots.append(Plot( name = 'l1_phi',
      texX = '#phi(l_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.W_l1['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    
    plots.append(Plot( name = 'l1_eta',
      texX = '#eta(l_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.W_l1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = 'l0_pt_charge',
      texX = 'p_{T}(l_{0}) [GeV] signed with lepton charge', texY = y_label,
      attribute = lambda event, sample: event.lep0chargept if event.passing_checks else float('nan'),
      binning=[20,-200,200],
    ))
    
    plots.append(Plot( name = 'l1_pt_charge',
      texX = 'p_{T}(l_{1}) [GeV] signed with lepton charge', texY = y_label,
      attribute = lambda event, sample: event.lep1chargept if event.passing_checks else float('nan'),
      binning=[20,-200,200],
    ))
    
    plots.append(Plot( name = 'njets',
      texX = 'Number of Jets', texY = y_label,
      attribute = lambda event, sample: getattr( event, 'n%sJet'%preTag ) if event.passing_checks else float('nan'),
      binning=[10,0,10],
    ))
    
    plots.append(Plot( name = 'nleps',
      texX = 'Number of Leptons', texY = y_label,
      attribute = lambda event, sample: getattr( event, 'n%sLep'%preTag ) if event.passing_checks else float('nan'),
      binning=[8,0,8],
    ))
    
    plots.append(Plot( name = 'nbjets',
      texX = 'Number of bJets', texY = y_label,
      attribute = lambda event, sample: len(event.bjets) if event.passing_checks else float('nan'),
      binning=[4,0,4],
    ))

    return plots
