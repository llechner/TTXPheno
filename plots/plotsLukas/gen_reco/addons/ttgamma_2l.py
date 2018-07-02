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
        "ref_lumiweight1fb/F",
        "lumiweight1fb/F",
        "genMet_pt/F", "genMet_phi/F",
    
        "ngenJet/I", "genJet[pt/F,eta/F,phi/F]",
        "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I]",
    
        "genPhoton_pt/F", "genPhoton_eta/F", "genPhoton_phi/F", "genPhoton_mass/F",
    
        "genBj0_pt/F", "genBj0_phi/F", "genBj0_eta/F",
        "genBj1_pt/F", "genBj1_phi/F", "genBj1_eta/F",

        "genBjLeadlep_index/I", "genBjLeadhad_index/I",
     ]
     
    # List of variables where genLep is replaced by reco for reco
    read_variables_genLep = [
    ]
    
    if level == 'reco':
        read_variables_gen    = [ variable.replace('gen', 'reco') for variable in read_variables_gen ]
        read_variables_genLep = [ variable.replace('genLep', 'reco') for variable in read_variables_genLep ]

    read_variables = read_variables_gen + read_variables_genLep
    read_variables = list( set( read_variables ) ) # remove double entries
    read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

    return read_variables


def makeJets( event, sample, level ):
    ''' Add a list of filtered jets to the event (full list is required for lepton cross cleaning)
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    # get (second) hardest bjets
    event.bj0 = {'pt':getattr( event, '%sBj0_pt'%preTag ), 'phi':getattr( event, '%sBj0_phi'%preTag ), 'eta':getattr( event, '%sBj0_eta'%preTag )}
    event.bj1 = {'pt':getattr( event, '%sBj1_pt'%preTag ), 'phi':getattr( event, '%sBj1_phi'%preTag ), 'eta':getattr( event, '%sBj1_eta'%preTag )}

    # Define leptonic b-jets
    event.bj0lep = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjLeadlep_index'%preTag ) ) if getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 else nanJet()
    event.bj0had = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjLeadhad_index'%preTag ) ) if getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 else nanJet()

    # Add extra vectors
    for p in [event.bj0, event.bj1, event.bj0lep, event.bj0had]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco':      from TTXPheno.Tools.objectSelection      import isGoodRecoJet     as isGoodJet
    else:                    from TTXPheno.Tools.objectSelection      import isGoodGenJet      as isGoodJet

    # selection checks
    event.Bj0 = isGoodJet( event.bj0 )
#    event.foundBjNonZlep = getattr( event, '%sBjNonZlep_index'%preTag ) >= 0 and isGoodJet( event.bjNonZlep )
#    event.foundBjNonZhad = getattr( event, '%sBjNonZhad_index'%preTag ) >= 0 and isGoodJet( event.bjNonZhad )
    event.foundBj0lep    = getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 and isGoodJet( event.bj0lep )
#    event.foundBj0had    = getattr( event, '%sBjLeadhad_index'%preTag ) >= 0 and isGoodJet( event.bj0had )

    # choose your selection on b-jets
    event.passing_bjets = event.foundBj0 and event.foundBj0lep


def makeMET( event, sample, level ):
    ''' Make a MET vector to facilitate further calculations
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    event.MET = {'pt':getattr(event, '%sMet_pt'%preTag), 'phi':getattr(event, '%sMet_phi'%preTag)}
    addTransverseVector( event.MET )


def makePhoton( event, sample, level ):
    ''' Make a Z vector to facilitate further calculations (either recoZ, genLepZ or genZ)
    '''
    event.gamma_unitVec2D = UnitVectorT2( getattr( event, '%sPhoton_phi'%level     ) )
    event.gamma_vec4D     = ROOT.TLorentzVector()
    event.gamma_vec4D.SetPtEtaPhiM( getattr( event, '%sPhoton_pt'%level ), getattr( event, '%sPhoton_eta'%level ), getattr( event, '%sPhoton_phi'%level ), getattr( event, '%sPhoton_mass'%level     ) )
    event.gamma_unitVec3D = event.gamma_vec4D.Vect().Unit()


def makeLeps( event, sample, level ):
    ''' Add important leptons (no full list of leptons is required for now)
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    # Define hardest leptons
    event.l0 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], 0 )
    event.l1 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], 1 )

    # Add extra vectors
    for p in [ event.l0, event.l1 ]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoLepton  as isGoodLepton
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenLepton   as isGoodLepton

    # We may loose some events by cross-cleaning or by thresholds.
    event.found1lep    = isGoodLepton( event.l0 )
    event.found2lep    = isGoodLepton( event.l1 )
    event.oppositeSign = event.l0['pdgId']*event.l1['pdgId'] > 0.

    # choose your selection on leptons
    event.passing_leptons = event.found1lep and event.found2lep and event.oppositeSign


def makeObservables( event, sample, level):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.bbdPhi = deltaPhi( event.bj0['phi'], event.bj1['phi'] )
    event.bbdR   = deltaR( event.bj0, event.bj1 )

    # double l kinematic
    event.lldPhi = deltaPhi( event.l0['phi'], event.l1['phi'] )
    event.lldR   = deltaR( event.l0, event.l1 )

    # signed lepton pt, Nan if len(event.lepsNotFromZ == 0)
    event.lep0chargept = event.l0['pt'] if event.l0['pdgId']>0 else -event.l0['pt']
    event.lep1chargept = event.l1['pt'] if event.l1['pdgId']>0 else -event.l1['pt']

    # choose your final selection
    event.passing_checks = event.passing_leptons and event.passing_bjets


def getSequenceList( level ):
    ''' sequence functions
    '''
    sequence = []

    sequence.append( lambda event, sample: makeJets( event, sample, level ) )
    sequence.append( lambda event, sample: makeMET( event, sample, level ) )
    sequence.append( lambda event, sample: makePhoton( event, sample, level ) )
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

    plots.append(Plot( name = 'njets',
      texX = 'Number of Jets', texY = 'Number of Events',
      attribute = lambda event, sample: getattr( event, 'n%sJet'%preTag ) if event.passing_checks else float('nan'),
      binning=[10,0,10],
    ))

    plots.append(Plot( name = 'nleps',
      texX = 'Number of Leptons', texY = 'Number of Events',
      attribute = lambda event, sample: getattr( event, 'n%sLep'%preTag ) if event.passing_checks else float('nan'),
      binning=[8,0,8],
    ))

    plots.append(Plot( name = "gamma_pt",
      texX = 'p_{T}(#gamma) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sPhoton_pt'%level ) if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name = "gamma_mass",
      texX = 'm(#gamma) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sPhoton_mass'%level ) if event.passing_checks else float('nan'),
      binning=[20,-5,5],
    ))
    
    plots.append(Plot( name = "b0_pt",
      texX = 'p_{T}(b_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.bj0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name = "b1_pt",
      texX = 'p_{T}(b_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.bj1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name = "b0_eta",
      texX = '#eta(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = "b1_eta",
      texX = '#eta(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = "b0_phi",
      texX = '#phi(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['phi'] if event.passing_checks else float('nan'),
      binning=[20,pi,pi],
    ))
    
    plots.append(Plot( name = "b1_phi",
      texX = '#phi(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['phi'] if event.passing_checks else float('nan'),
      binning=[20,pi,pi],
    ))
    
    plots.append(Plot( name = "l0_pt",
      texX = 'p_{T}(l_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.l0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,300],
    ))
    
    plots.append(Plot( name = "l0_eta",
      texX = '#eta(l_{0})', texY = y_label,
      attribute = lambda event, sample: event.l0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = "l0_phi",
      texX = '#phi(l_{0})', texY = y_label,
      attribute = lambda event, sample: event.l0['phi'] if event.passing_checks else float('nan'),
      binning=[20,pi,pi],
    ))
    
    plots.append(Plot( name = "l1_pt",
      texX = 'p_{T}(l_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.l1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,300],
    ))
    
    plots.append(Plot( name = "l1_eta",
      texX = '#eta(l_{1})', texY = y_label,
      attribute = lambda event, sample: event.l1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = "l1_phi",
      texX = '#phi(l_{1})', texY = y_label,
      attribute = lambda event, sample: event.l1['phi'] if event.passing_checks else float('nan'),
      binning=[20,pi,pi],
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
      texX = 'E_{T}^{miss} [GeV]', texY = y_label,
      attribute = lambda event, sample: event.MET['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name	= '2nu_Met_phi',
      texX = '#phi(E_{T}^{miss})', texY = y_label,
      attribute = lambda event, sample: event.MET['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
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
    
    return plots
