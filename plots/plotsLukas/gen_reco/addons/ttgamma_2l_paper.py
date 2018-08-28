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
from TTXPheno.Tools.helpers              import deltaPhi, getCollection, deltaR, nanJet, nanLepton, nanPhoton, getObjDict

# Import samples
from TTXPheno.samples.benchmarks         import *

# Import helpers
from plot_helpers                        import *


def getVariableList( level ):
   
    # List of variables where gen is replaced by reco for reco
    read_variables_gen = [
#        "ref_lumiweight1fb/F",
        "lumiweight1fb/F",
#        "nSignalPhotons/I",
        "genMet_pt/F", "genMet_phi/F",
    
        "ngenJet/I", #"genJet[pt/F,eta/F,phi/F]",
        "ngenLep/I", #"genLep[pt/F,eta/F,phi/F,pdgId/I]",
    
        "ngenPhoton/I", #"genPhoton[pt/F,phi/F,eta/F,mass/F,motherPdgId/I]"
#        "genPhoton_pt/F", "genPhoton_eta/F", "genPhoton_phi/F",
    
        "genBj0_pt/F", "genBj0_phi/F", "genBj0_eta/F",
        "genBj1_pt/F", "genBj1_phi/F", "genBj1_eta/F",

        "genBjLeadlep_index/I", "genBjLeadhad_index/I",

        "genZ_mass/F",
     ]
     
    # List of variables where genLep is replaced by reco for reco
    read_variables_genLep = [
    ]
    
    if level == 'reco':
        read_variables_gen    = [ variable.replace('gen', 'reco') for variable in read_variables_gen ]
        read_variables_genLep = [ variable.replace('genLep', 'reco') for variable in read_variables_genLep ]
        read_variables_gen.append("recoJet[pt/F,eta/F,phi/F,bTag/F]")
        read_variables_gen.append("recoPhoton[pt/F,eta/F,phi/F,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F,genIndex/I,minLeptonPt/F,minLeptonDR/F,minJetDR/F]")
        read_variables_gen.append("recoLep[pt/F,eta/F,phi/F,pdgId/I,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F]")
        read_variables_gen.append("genPhoton[motherPdgId/I,relIso04/F]")
    else:
        read_variables_gen.append("genJet[pt/F,eta/F,phi/F,matchBParton/I]")
        read_variables_gen.append("genPhoton[pt/F,phi/F,eta/F,mass/F,motherPdgId/I,relIso04/F,minLeptonDR/F,minJetDR/F]")
        read_variables_gen.append("genLep[pt/F,phi/F,eta/F,pdgId/I,motherPdgId/I]")

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

    # Define leptonic b-jets
    event.bj0lep = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjLeadlep_index'%preTag ) ) if getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 else nanJet()
#    event.bj0had = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjLeadhad_index'%preTag ) ) if getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 else nanJet()

    # Add extra vectors
    for p in [event.bj0, event.bj1, event.bj0lep]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco':      from TTXPheno.Tools.objectSelection      import isGoodRecoJet     as isGoodJet
    else:                    from TTXPheno.Tools.objectSelection      import isGoodGenJet      as isGoodJet

    # selection checks
    event.foundBj0 = isGoodJet( event.bj0 )
#    event.foundBj0lep    = getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 and isGoodJet( event.bj0lep )

    # choose your selection on b-jets
    event.passing_bjets = event.foundBj0


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
    preTag = 'reco' if level == 'reco' else 'gen'

    if level == 'reco':
        photonList = ['pt', 'eta', 'phi', 'isolationVar', 'isolationVarRhoCorr', 'sumPtCharged', 'sumPtNeutral', 'sumPtChargedPU', 'sumPt', 'ehadOverEem', 'genIndex', 'minLeptonDR', 'minLeptonPt', 'minJetDR']
    else:
        photonList = ['pt', 'eta', 'phi', 'mass', 'motherPdgId', 'relIso04', 'minLeptonDR', 'minJetDR'] #'minLeptonPt'

    event.gammas = getCollection( event, '%sPhoton'%preTag, photonList, 'n%sPhoton'%preTag )

    event.gamma0 = event.gammas[0] if len(event.gammas)>0 else nanPhoton()
    event.gamma1 = event.gammas[1] if len(event.gammas)>1 else nanPhoton()
    event.gamma2 = event.gammas[2] if len(event.gammas)>2 else nanPhoton()
#    event.gamma0 = getObjDict( event, '%sPhoton_'%preTag, photonList, 0 )
#    try: event.gamma1 = getObjDict( event, '%sPhoton_'%preTag, photonList, 1 )
#    except: event.gamma1 = nanPhoton()
#    try: event.gamma2 = getObjDict( event, '%sPhoton_'%preTag, photonList, 2 )
#    except: event.gamma2 = nanPhoton()

    if level == 'reco':
        for p in [event.gamma2, event.gamma1, event.gamma0]:
            addIDeltaBeta( p )

    event.gamma_unitVec2D = UnitVectorT2( event.gamma0['phi'] )
    event.gamma_vec4D     = ROOT.TLorentzVector()
    event.gamma_vec4D.SetPtEtaPhiM( event.gamma0['pt'], event.gamma0['eta'], event.gamma0['phi'], 0 )
    event.gamma_unitVec3D = event.gamma_vec4D.Vect().Unit()

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco':      from TTXPheno.Tools.objectSelection      import isGoodRecoPhoton     as isGoodPhoton
    else:                    from TTXPheno.Tools.objectSelection      import isGoodGenPhoton      as isGoodPhoton

    event.foundGamma0    = isGoodPhoton( event.gamma0 )
    event.foundGamma1    = isGoodPhoton( event.gamma1 )
    event.foundGamma2    = isGoodPhoton( event.gamma2 )

    event.passing_photons =  event.foundGamma0

def makeLeps( event, sample, level, flavorCheck ):
    ''' Add important leptons (no full list of leptons is required for now)
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    if level == 'reco':
        leptonList = ['pt', 'eta', 'phi', 'pdgId', 'isolationVar', 'isolationVarRhoCorr', 'sumPtCharged', 'sumPtNeutral', 'sumPtChargedPU', 'sumPt', 'ehadOverEem', 'genIndex']
    else:
        leptonList = ['pt', 'eta', 'phi', 'pdgId', 'motherPdgId']

    event.leps = getCollection( event, '%sLep'%preTag, leptonList, 'n%sLep'%preTag )

    # Define hardest leptons
    event.l0 = event.leps[0]
    event.l1 = event.leps[1]

    # Add extra vectors
    for p in [ event.l0, event.l1 ]:
        addTransverseVector( p )
        addTLorentzVector( p )
        if level == 'reco':
            addIDeltaBeta( p )


    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoLepton  as isGoodLepton
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenLepton   as isGoodLepton

    # We may loose some events by cross-cleaning or by thresholds.
    event.found1lep    = isGoodLepton( event.l0 )
    event.found2lep    = isGoodLepton( event.l1 ) and len(event.leps) == 2
    event.oppositeSign = event.l0['pdgId']*event.l1['pdgId'] < 0.
    event.sameFlavor   = abs(event.l0['pdgId']) == abs(event.l1['pdgId'])

    # choose your selection on leptons
    event.passing_leptons = event.found1lep and event.found2lep and event.oppositeSign
    if flavorCheck == 'same': event.passing_leptons = event.passing_leptons and event.sameFlavor
    elif flavorCheck == 'opposite': event.passing_leptons = event.passing_leptons and not event.sameFlavor


def makeObservables( event, sample, level ):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.bbdPhi = deltaPhi( event.bj0['phi'], event.bj1['phi'] )
    event.bbdR   = deltaR( event.bj0, event.bj1 )

    # lep gamma kinematic
    event.l0gammadPhi = deltaPhi( event.l0['phi'], event.gamma0['phi'] )
    event.l0gammadR   = deltaR( event.l0, event.gamma0 )

    event.l1gammadPhi = deltaPhi( event.l1['phi'], event.gamma0['phi'] )
    event.l1gammadR   = deltaR( event.l1, event.gamma0 )

    # double l kinematic
    event.lldPhi = deltaPhi( event.l0['phi'], event.l1['phi'] )
    event.lldR   = deltaR( event.l0, event.l1 )

    event.minLeptonG0dR   = min([ deltaR( lepton, event.gamma0 ) for lepton in event.leps ]) if len(event.leps)>0 and event.gamma0['pt']==event.gamma0['pt'] else float('nan')
    event.minLeptonG1dR   = min([ deltaR( lepton, event.gamma1 ) for lepton in event.leps ]) if len(event.leps)>0 and event.gamma1['pt']==event.gamma1['pt'] else float('nan')
    event.minLeptonG2dR   = min([ deltaR( lepton, event.gamma2 ) for lepton in event.leps ]) if len(event.leps)>0 and event.gamma2['pt']==event.gamma2['pt'] else float('nan')

    event.minJetG0dR      = min([ deltaR( jet, event.gamma0 ) for jet in event.jets ]) if len(event.jets)>0 and event.gamma0['pt']==event.gamma0['pt'] else float('nan')
    event.minJetG1dR      = min([ deltaR( jet, event.gamma1 ) for jet in event.jets ]) if len(event.jets)>0 and event.gamma1['pt']==event.gamma1['pt'] else float('nan')
    event.minJetG2dR      = min([ deltaR( jet, event.gamma2 ) for jet in event.jets ]) if len(event.jets)>0 and event.gamma2['pt']==event.gamma2['pt'] else float('nan')

    # signed lepton pt, Nan if len(event.lepsNotFromZ == 0)
    event.lep0chargept = event.l0['pt'] if event.l0['pdgId']>0 else -event.l0['pt']
    event.lep1chargept = event.l1['pt'] if event.l1['pdgId']>0 else -event.l1['pt']

    event.mll = (event.l0['vec4D'] + event.l1['vec4D']).M()
    event.mllgamma = (event.l0['vec4D'] + event.l1['vec4D'] + event.gamma_vec4D).M()

    #cut the Z window in same flavor lepton region
    event.offZ = (abs(event.mll - 91.2) > 10 and abs(event.mllgamma - 91.2) > 10) or event.mll is float('nan') or event.mllgamma is float('nan') or abs(event.l0['pdgId']) != abs(event.l1['pdgId'])
    event.mllCut = event.mll > 40

    # choose your final selection
    event.passing_checks = event.passing_leptons and event.passing_bjets and event.passing_photons and event.mllCut and event.offZ


def getSequenceList( level, flavorCheck ):
    ''' sequence functions
    '''
    sequence = []

    sequence.append( lambda event, sample: makeJets( event, sample, level ) )
    sequence.append( lambda event, sample: makeMET( event, sample, level ) )
    sequence.append( lambda event, sample: makePhoton( event, sample, level ) )
    sequence.append( lambda event, sample: makeLeps( event, sample, level, flavorCheck ) )
    sequence.append( lambda event, sample: makeObservables( event, sample, level ) )

    return sequence
    
def getPlotList( scaleLumi, level ):

    tag = 'reco' if level == 'reco' else 'genLep'
    preTag = 'reco' if level == 'reco' else 'gen'

#    if scaleLumi: y_label = 'norm. diff. xsec'
#    else:         y_label = 'diff. x-sec'

#    labelAddon = '#frac{1}{#sigma} ' if scaleLumi else ''
    labelAddon = '1/#sigma ' if scaleLumi else ''
#    unit = ' [GeV^{-1}]' if scaleLumi else ' [fb/GeV]'
    unit = ''

    plots = []
    fisherInfoVariables = []
    
    plots.append(Plot( name = 'l0gammaDPhi',
      texX = '#Delta#phi(l_{0},#gamma_{0})', texY = labelAddon + 'd#sigma/d#Delta#phi(l_{0},#gamma_{0})'+unit,
#      texX = '#Delta#phi(l_{0},#gamma_{0})', texY = labelAddon + '#frac{d#sigma}{d #Delta#phi(l_{0},#gamma_{0})}'+unit,
      attribute = lambda event, sample: event.l0gammadPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'l0gammaDR',
      texX = '#Delta R(l_{0},#gamma_{0})', texY = labelAddon + 'd#sigma/d#Delta R(l_{0},#gamma_{0})'+unit,
#      texX = '#Delta R(l_{0},#gamma_{0})', texY = labelAddon + '#frac{d#sigma}{d #Delta R(l_{0},#gamma_{0})}'+unit,
      attribute = lambda event, sample: event.l0gammadR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'l1gammaDPhi',
      texX = '#Delta#phi(l_{1},#gamma_{0})', texY = labelAddon + 'd#sigma/d#Delta#phi(l_{1},#gamma_{0})'+unit,
#      texX = '#Delta#phi(l_{1},#gamma_{0})', texY = labelAddon + '#frac{d#sigma}{d #Delta#phi(l_{1},#gamma_{0})}'+unit,
      attribute = lambda event, sample: event.l1gammadPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'l1gammaDR',
      texX = '#Delta R(l_{1},#gamma_{0})', texY = labelAddon + 'd#sigma/d#Delta R(l_{1},#gamma_{0})'+unit,
#      texX = '#Delta R(l_{1},#gamma_{0})', texY = labelAddon + '#frac{d#sigma}{d #Delta R(l_{1},#gamma_{0})}'+unit,
      attribute = lambda event, sample: event.l1gammadR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = "gamma_pt10",
      texX = 'p_{T}(#gamma_{0}) [GeV]', texY = labelAddon + 'd#sigma/dp_{T}(#gamma_{0})'+unit,
#      texX = 'p_{T}(#gamma_{0}) [GeV]', texY = labelAddon + '#frac{d#sigma}{d p_{T}(#gamma_{0})}'+unit,
      attribute = lambda event, sample: event.gamma0['pt'] if event.passing_checks else float('nan'),
      binning=[10,0,400],
    ))
    fisherInfoVariables.append('%sPhoton_pt[0]'%preTag)

    plots.append(Plot( name = "gamma_pt20",
      texX = 'p_{T}(#gamma_{0}) [GeV]', texY = labelAddon + 'd#sigma/dp_{T}(#gamma_{0})'+unit,
#      texX = 'p_{T}(#gamma_{0}) [GeV]', texY = labelAddon + '#frac{d#sigma}{d p_{T}(#gamma_{0})}'+unit,
      attribute = lambda event, sample: event.gamma0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    fisherInfoVariables.append('%sPhoton_pt[0]'%preTag)

    plots.append(Plot( name = "gamma_eta",
      texX = '#eta(#gamma_{0})', texY = labelAddon + 'd#sigma/d#eta(#gamma_{0})'+unit,
      attribute = lambda event, sample: event.gamma0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    fisherInfoVariables.append('%sPhoton_eta[0]'%preTag)

    plots.append(Plot( name = "l1_pt",
      texX = 'p_{T}(l_{1}) [GeV]', texY = labelAddon + 'd#sigma/dp_{T}(l_{1})'+unit,
#      texX = 'p_{T}(l_{1}) [GeV]', texY = labelAddon + '#frac{d#sigma}{d p_{T}(l_{1})}'+unit,
      attribute = lambda event, sample: event.l1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,100],
    ))
    fisherInfoVariables.append('%sLep_pt[1]'%preTag)

    plots.append(Plot( name = "l0_pt",
      texX = 'p_{T}(l_{0}) [GeV]', texY = labelAddon + 'd#sigma/dp_{T}(l_{0})'+unit,
#      texX = 'p_{T}(l_{0}) [GeV]', texY = labelAddon + '#frac{d#sigma}{d p_{T}(l_{0})}'+unit,
      attribute = lambda event, sample: event.l0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ))
    fisherInfoVariables.append('%sLep_pt[0]'%preTag)

    plots.append(Plot( name = 'deltaPhi_bb',
      texX = '#Delta#phi(bb)', texY = labelAddon + 'd#sigma/d#Delta#phi(bb)'+unit,
#      texX = '#Delta#phi(bb)', texY = labelAddon + '#frac{d#sigma}{d #Delta#phi(bb)}'+unit,
      attribute = lambda event, sample: event.bbdPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'deltaR_bb',
      texX = '#DeltaR(bb)', texY = labelAddon + 'd#sigma/d#DeltaR(bb)'+unit,
#      texX = '#DeltaR(bb)', texY = labelAddon + '#frac{d#sigma}{d #DeltaR(bb)}'+unit,
      attribute = lambda event, sample: event.bbdR if event.passing_checks else float('nan'),
      binning=[20,0,3.5],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'deltaPhi_ll',
      texX = '#Delta#phi(ll)', texY = labelAddon + 'd#sigma/d#Delta#phi(ll)'+unit,
#      texX = '#Delta#phi(ll)', texY = labelAddon + '#frac{d#sigma}{d #Delta#phi(ll)}'+unit,
      attribute = lambda event, sample: event.lldPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'deltaR_ll',
      texX = '#DeltaR(ll)', texY = labelAddon + 'd#sigma/d#DeltaR(ll)'+unit,
#      texX = '#DeltaR(ll)', texY = labelAddon + '#frac{d#sigma}{d #DeltaR(ll)}'+unit,
      attribute = lambda event, sample: event.lldR if event.passing_checks else float('nan'),
      binning=[20,0,3.5],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'Met_pt',
      texX = 'E_{T}^{miss} [GeV]', texY = labelAddon + 'd#sigma/dE_{T}^{miss}'+unit,
#      texX = 'E_{T}^{miss} [GeV]', texY = labelAddon + '#frac{d#sigma}{d E_{T}^{miss}}'+unit,
      attribute = lambda event, sample: event.MET['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ))
    fisherInfoVariables.append('%sMet_pt'%preTag)

    return plots, fisherInfoVariables

