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
    
        "genZ_pt/F", "genZ_eta/F", "genZ_phi/F", "genZ_mass/F", "genZ_cosThetaStar/F",
    
        "genBj0_pt/F", "genBj0_phi/F", "genBj0_eta/F",
        "genBj1_pt/F", "genBj1_phi/F", "genBj1_eta/F",
    
        "genBjLeadlep_index/I", "genBjLeadhad_index/I",
        "genBjNonZlep_index/I", "genBjNonZhad_index/I",
    ]
     
    # List of variables where genLep is replaced by reco for reco
    read_variables_genLep = [
        "genLepZ_pt/F", "genLepZ_eta/F", "genLepZ_phi/F", "genLepZ_mass/F", "genLepZ_cosThetaStar/F",
        "genLepZ_lldPhi/F", "genLepZ_lldR/F",
    
        "genLepZ_l1_index/I", "genLepZ_l2_index/I",
        "genLepNonZ_l1_index/I", "genLepNonZ_l2_index/I",
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
    for p in [event.bj0, event.bj1 ]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoJet     as isGoodJet
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenJet      as isGoodJet

    # selection checks
    event.foundBj0       = isGoodJet( event.bj0 )
#    event.foundBj1       = isGoodJet( event.bj1 )
#    event.foundBjNonZlep = getattr( event, '%sBjNonZlep_index'%preTag ) >= 0 and isGoodJet( event.bjNonZlep )
#    event.foundBjNonZhad = getattr( event, '%sBjNonZhad_index'%preTag ) >= 0 and isGoodJet( event.bjNonZhad )
#    event.foundBj0lep    = getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 and isGoodJet( event.bj0lep )
#    event.foundBj0had    = getattr( event, '%sBjLeadhad_index'%preTag ) >= 0 and isGoodJet( event.bj0had )

    # choose your selection on b-jets
    event.passing_bjets = event.foundBj0


def makeMET( event, sample, level ):
    ''' Make a MET vector to facilitate further calculations
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    event.MET = {'pt':getattr(event, '%sMet_pt'%preTag), 'phi':getattr(event, '%sMet_phi'%preTag)}
    addTransverseVector( event.MET )


def makeZ( event, sample, level ):
    ''' Make a Z vector to facilitate further calculations (either recoZ, genLepZ or genZ)
    '''
    event.Z_unitVec2D = UnitVectorT2( getattr( event, '%sZ_phi'%level     ) )
    event.Z_vec4D     = ROOT.TLorentzVector()
    event.Z_vec4D.SetPtEtaPhiM( getattr( event, '%sZ_pt'%level ), getattr( event, '%sZ_eta'%level ), getattr( event, '%sZ_phi'%level ), getattr( event, '%sZ_mass'%level     ) )
    event.Z_unitVec3D = event.Z_vec4D.Vect().Unit()


def makeLeps( event, sample, level, flavorCheck ):
    ''' Add important leptons (no full list of leptons is required for now)
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    event.leps = getCollection( event, '%sLep'%preTag, ['pt', 'eta', 'phi', 'pdgId', 'motherPdgId'], 'n%sLep'%preTag )
    event.lepsFromZ = list( filter( lambda j: j['motherPdgId'] == 23, event.leps ) )
    event.oneZ = len( event.lepsFromZ ) == 2

    # Define Z leptons
    event.Z_l0 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], getattr( event, '%sZ_l1_index'%tag     ) ) if getattr( event, '%sZ_l1_index'%tag ) >= 0 else nanLepton()
    event.Z_l1 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], getattr( event, '%sZ_l2_index'%tag     ) ) if getattr( event, '%sZ_l2_index'%tag ) >= 0 else nanLepton()

    # Define non-Z leptons
    event.NonZ_l0 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], getattr( event, '%sNonZ_l1_index'%tag     ) ) if getattr( event, '%sNonZ_l1_index'%tag ) >= 0 else nanLepton()
    event.NonZ_l1 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], getattr( event, '%sNonZ_l2_index'%tag     ) ) if getattr( event, '%sNonZ_l2_index'%tag ) >= 0 else nanLepton()

    # Add extra vectors
    for p in [ event.Z_l0, event.Z_l1, event.NonZ_l0, event.NonZ_l1 ]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoLepton  as isGoodLepton
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenLepton   as isGoodLepton

    # We may loose some events by cross-cleaning or by thresholds.
    event.foundZl0     = getattr( event, '%sZ_l1_index'%tag ) >= 0 and isGoodLepton( event.Z_l0 )
    event.foundZl1     = getattr( event, '%sZ_l2_index'%tag ) >= 0 and isGoodLepton( event.Z_l1 )
    event.foundZ       = event.Z_l0['pdgId'] * event.Z_l1['pdgId'] < 0 and abs(event.Z_l0['pdgId']) == abs(event.Z_l1['pdgId'])
    event.found3lep    = getattr( event, '%sNonZ_l1_index'%tag ) >= 0 and isGoodLepton( event.NonZ_l0 )
    event.found4lep    = getattr( event, '%sNonZ_l2_index'%tag ) >= 0 and isGoodLepton( event.NonZ_l1 )
    event.oppositeSign = event.NonZ_l0['pdgId'] * event.NonZ_l1['pdgId'] < 0
    event.sameFlavor   = abs(event.NonZ_l0['pdgId']) == abs(event.NonZ_l1['pdgId'])

    # choose your selection on leptons
    event.passing_leptons = event.found3lep and event.found4lep and event.foundZl0 and event.foundZl1 and event.foundZ and event.oppositeSign
    if flavorCheck == 'same': event.passing_leptons = event.passing_leptons and event.sameFlavor
    elif flavorCheck == 'opposite': event.passing_leptons = event.passing_leptons and not event.sameFlavor


def makeObservables( event, sample, level):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.bbdPhi = deltaPhi( event.bj0['phi'], event.bj1['phi'] )
    event.bbdR   = deltaR( event.bj0, event.bj1 )

    # double l kinematic
    event.nonZ_lldPhi = deltaPhi( event.NonZ_l0['phi'], event.NonZ_l1['phi'] )
    event.nonZ_lldR = deltaR( event.NonZ_l0, event.NonZ_l1 )

    # signed lepton pt, Nan if len(event.lepsNotFromZ == 0)
    event.nonZlep0chargept = event.NonZ_l0['pt'] if event.NonZ_l0['pdgId']>0 else -event.NonZ_l0['pt']
    event.nonZlep1chargept = event.NonZ_l1['pt'] if event.NonZ_l1['pdgId']>0 else -event.NonZ_l1['pt']

    # choose your final selection
    event.passing_checks = event.passing_leptons and event.passing_bjets


def getSequenceList( level, flavorCheck ):
    ''' sequence functions
    '''
    sequence = []

    sequence.append( lambda event, sample: makeJets( event, sample, level ) )
    sequence.append( lambda event, sample: makeMET( event, sample, level ) )
    sequence.append( lambda event, sample: makeZ( event, sample, level ) )
    sequence.append( lambda event, sample: makeLeps( event, sample, level, flavorCheck ) )
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

    plots.append( Plot( name = "Z_pt",
      texX = 'p_{T}(Z) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sZ_pt'%level ) if event.passing_checks else float('nan'),
      binning=[10,0,500],
    ) )
    
    plots.append( Plot( name = 'Z_phi',
      texX = '#phi(Z) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sZ_phi'%level ) if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ) )
    
    plots.append( Plot( name = 'Z_eta',
      texX = '#eta(Z) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sZ_eta'%level ) if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ) )
    
    plots.append( Plot( name = "Z_cosThetaStar",
      texX = 'cos(#theta*)', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sZ_cosThetaStar'%level ) if event.passing_checks else float('nan'),
      binning=[20,-1.2,1.2],
    ) )
    
    plots.append( Plot( name = "Z_mass",
      texX = 'm(ll) [GeV]', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sZ_mass'%level ) if event.passing_checks else float('nan'),
      binning=[20,70,110],
    ) )
    
    plots.append( Plot( name = "b0_pt",
      texX = 'p_{T}(b_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.bj0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ) )
    
    plots.append( Plot( name = "b0_phi",
      texX = '#phi(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ) )
    
    plots.append( Plot( name = "b0_eta",
      texX = '#eta(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ) )
    
    plots.append( Plot( name = "b1_pt",
      texX = 'p_{T}(b_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.bj1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ) )
    
    plots.append( Plot( name = "b1_phi",
      texX = '#phi(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ) )
    
    plots.append( Plot( name = "b1_eta",
      texX = '#eta(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ) )
    
    plots.append( Plot( name = 'deltaPhi_bb',
      texX = '#Delta#phi(bb)', texY = y_label,
      attribute = lambda event, sample: event.bbdPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ) )
    
    plots.append( Plot( name = 'deltaR_bb',
      texX = '#DeltaR(bb)', texY = y_label,
      attribute = lambda event, sample: event.bbdR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ) )
        
    plots.append( Plot( name = "nonZl0_pt",
      texX = 'p_{T}(l^{non-Z}_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,300],
    ) )
    
    plots.append( Plot( name = "nonZl0_phi",
      texX = '#phi(l^{non-Z}_{0})', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l0['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ) )
    
    plots.append( Plot( name = "nonZl0_eta",
      texX = '#eta(l^{non-Z}_{0})', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ) )
    
    plots.append( Plot( name = "nonZl1_pt",
      texX = 'p_{T}(l^{non-Z}_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ) )
    
    plots.append( Plot( name = "nonZl1_phi",
      texX = '#phi(l^{non-Z}_{1})', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l1['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ) )
    
    plots.append( Plot( name = "nonZl1_eta",
      texX = '#eta(l^{non-Z}_{1})', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ) )
        
    plots.append( Plot( name = 'nonZ_deltaPhi_ll',
      texX = '#Delta#phi(ll)^{non-Z}', texY = y_label,
      attribute = lambda event, sample: event.nonZ_lldPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ) )
    
    plots.append( Plot( name = 'nonZ_deltaR_ll',
      texX = '#DeltaR(ll)^{non-Z}', texY = y_label,
      attribute = lambda event, sample: event.nonZ_lldR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ) )
    
    # for gen use genLepZ_lldPhi
    plots.append( Plot( name = 'Z_deltaPhi_ll',
      texX = '#Delta#phi(ll)', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sZ_lldPhi'%tag ) if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ) )
        
    # for gen use genLepZ_lldR
    plots.append( Plot( name = 'Z_deltaR_ll',
      texX = '#DeltaR(ll)', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sZ_lldR'%tag ) if event.passing_checks else float('nan'),
      binning=[20,0,4],
    ) )
        
    plots.append( Plot( name = '2nu_Met_pt',
      texX = 'E_{T}^{miss} [GeV]', texY = y_label,
      attribute = lambda event, sample: event.MET['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ) )
    
    plots.append( Plot( name	= '2nu_Met_phi',
      texX = '#phi(E_{T}^{miss})', texY = y_label,
      attribute = lambda event, sample: event.MET['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ) )
    
    plots.append( Plot( name = 'l0nonZ_pt_charge',
      texX = 'p_{T}(l^{non-Z}) [GeV] signed with lepton charge', texY = y_label,
      attribute = lambda event, sample: event.nonZlep0chargept if event.passing_checks else float('nan'),
      binning=[20,-200,200],
    ) )
    
    plots.append( Plot( name = 'l1nonZ_pt_charge',
      texX = 'p_{T}(l^{non-Z}) [GeV] signed with lepton charge', texY = y_label,
      attribute = lambda event, sample: event.nonZlep1chargept if event.passing_checks else float('nan'),
      binning=[20,-200,200],
    ) )

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

