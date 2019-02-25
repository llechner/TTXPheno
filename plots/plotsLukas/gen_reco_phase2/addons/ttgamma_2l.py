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
        "ngenTop/I", 
    
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
        read_variables_gen.append("genTop[pt/F,eta/F,phi/F,mass/F,pdgId/I]")
        read_variables_gen.append("genPhoton[pt/F,phi/F,eta/F,mass/F,motherPdgId/I,relIso04/F,minLeptonDR/F,minJetDR/F]")
        read_variables_gen.append("genLep[pt/F,phi/F,eta/F,pdgId/I,motherPdgId/I,grandmotherPdgId/I]")

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
        leptonList = ['pt', 'eta', 'phi', 'pdgId', 'motherPdgId', 'grandmotherPdgId']


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


def makeSpinCorrelationObservables( event, sample, level ):

    if level != 'gen': return
    event.tops = getCollection( event, 'genTop', ['pt', 'eta', 'phi', 'pdgId' , 'mass'], 'ngenTop' )
    event.c_nn    = float('nan')
    event.c_rr    = float('nan')
    event.c_kk    = float('nan')
    event.c_rkpkr = float('nan')
    event.c_nrprn = float('nan')
    event.c_nkpkn = float('nan')
    event.c_rkmkr = float('nan')
    event.c_nrmrn = float('nan')
    event.c_nkmkn = float('nan')
    event.b_ntop  = float('nan')
    event.b_ntbar = float('nan')
    event.b_rtop  = float('nan')
    event.b_rtbar = float('nan')
    event.b_ktop  = float('nan')
    event.b_ktbar = float('nan')
    event.b_rstartop  = float('nan')
    event.b_rstartbar = float('nan')
    event.b_kstartop  = float('nan')
    event.b_kstartbar = float('nan')

    if len(event.tops)>=2 \
            and event.l0['pdgId']*event.l1['pdgId']<0\
            and abs(event.l0['motherPdgId'])==abs(event.l1['motherPdgId']) == 24\
            and abs(event.l0['grandmotherPdgId'])==abs(event.l1['grandmotherPdgId']) == 6:

        if event.tops[0]['pdgId'] <0:
            event.top_pos, event.top_neg = event.tops[:2]
        else:
            event.top_neg, event.top_pos = event.tops[:2]
            
        if event.l0['pdgId']<0:
            event.lep_pos, event.lep_neg = event.l0, event.l1
        else:
            event.lep_neg, event.lep_pos = event.l0, event.l1

        # fill top/antitip with event.top_pos, event.top_neg
        top_lab = ROOT.TLorentzVector()
        top_lab.SetPtEtaPhiM(event.top_neg['pt'], event.top_neg['eta'], event.top_neg['phi'], event.top_neg['mass'])
        tbar_lab = ROOT.TLorentzVector()
        tbar_lab.SetPtEtaPhiM(event.top_pos['pt'], event.top_pos['eta'], event.top_pos['phi'], event.top_pos['mass'])
        # create boost to ttbar
        ttbar_lab = ROOT.TLorentzVector(top_lab + tbar_lab)
        boosttottbar = -ttbar_lab.BoostVector()
        # apply ttbar boost to top/tbar
        top_ZMFttbar = ROOT.TLorentzVector(top_lab)
        top_ZMFttbar.Boost(boosttottbar)
        tbar_ZMFttbar = ROOT.TLorentzVector(tbar_lab)
        tbar_ZMFttbar.Boost(boosttottbar)
        # apply ttbar boost to leptons
        lep_ZMFtbar = ROOT.TLorentzVector()
        lep_ZMFtbar.SetPtEtaPhiM(event.lep_neg['pt'], event.lep_neg['eta'], event.lep_neg['phi'], 0. )
        lep_ZMFtbar.Boost(boosttottbar)
        boosttotbar = -tbar_ZMFttbar.BoostVector()
        lep_ZMFtbar.Boost(boosttotbar)
        # anti lepton
        lbar_ZMFtop = ROOT.TLorentzVector()
        lbar_ZMFtop.SetPtEtaPhiM(event.lep_pos['pt'], event.lep_pos['eta'], event.lep_pos['phi'], 0. )
        lbar_ZMFtop.Boost(boosttottbar)
        boosttotop = -top_ZMFttbar.BoostVector()
        lbar_ZMFtop.Boost(boosttotop)
        # caculate correlation obervables
        # cos between leptons
        #event.cosphi = cos(lep_ZMFtbar.Angle(lbar_ZMFtop.Vect()))
        # define bases
        pbone = ROOT.TVector3(0.,0.,1.)
        kbone = ROOT.TVector3(top_ZMFttbar.Vect().Unit())
        # observables based on https://arxiv.org/pdf/1508.05271v2.pdf
        yp = pbone.Dot(kbone)
        rp = sqrt(1-yp*yp)
        signyp = 1.
        if yp<0 :
            signyp = -1.
        rbtwo = ROOT.TVector3(1/rp*(pbone-yp*kbone))
        nbtwo = ROOT.TVector3(1/rp*(pbone.Cross(kbone)))
        lbar_ZMFtop_V3U = ROOT.TVector3(lbar_ZMFtop.Vect().Unit())
        lep_ZMFtbar_V3U = ROOT.TVector3(lep_ZMFtbar.Vect().Unit())
        # calculate spin correlation coefficients and add variables to event 
        event.c_nn    = lbar_ZMFtop_V3U.Dot(signyp*nbtwo) * lep_ZMFtbar_V3U.Dot(-1*signyp*nbtwo)
        event.c_rr    = lbar_ZMFtop_V3U.Dot(signyp*rbtwo) * lep_ZMFtbar_V3U.Dot(-1*signyp*rbtwo)
        event.c_kk    = lbar_ZMFtop_V3U.Dot(kbone) * lep_ZMFtbar_V3U.Dot(-1*kbone)
        event.c_rkpkr = lbar_ZMFtop_V3U.Dot(signyp*rbtwo) * lep_ZMFtbar_V3U.Dot(-1*kbone) + lbar_ZMFtop_V3U.Dot(kbone) * lep_ZMFtbar_V3U.Dot(-1*signyp*rbtwo)
        event.c_nrprn = lbar_ZMFtop_V3U.Dot(signyp*nbtwo) * lep_ZMFtbar_V3U.Dot(-1*signyp*rbtwo) + lbar_ZMFtop_V3U.Dot(signyp*rbtwo) * lep_ZMFtbar_V3U.Dot(-1*signyp*nbtwo)
        event.c_nkpkn = lbar_ZMFtop_V3U.Dot(signyp*nbtwo) * lep_ZMFtbar_V3U.Dot(-1*kbone) + lbar_ZMFtop_V3U.Dot(kbone) * lep_ZMFtbar_V3U.Dot(-1*signyp*nbtwo)
        event.c_rkmkr = lbar_ZMFtop_V3U.Dot(signyp*rbtwo) * lep_ZMFtbar_V3U.Dot(-1*kbone) - lbar_ZMFtop_V3U.Dot(kbone) * lep_ZMFtbar_V3U.Dot(-1*signyp*rbtwo)
        event.c_nrmrn = lbar_ZMFtop_V3U.Dot(signyp*nbtwo) * lep_ZMFtbar_V3U.Dot(-1*signyp*rbtwo) - lbar_ZMFtop_V3U.Dot(signyp*rbtwo) * lep_ZMFtbar_V3U.Dot(-1*signyp*nbtwo)
        event.c_nkmkn = lbar_ZMFtop_V3U.Dot(signyp*nbtwo) * lep_ZMFtbar_V3U.Dot(-1*kbone) - lbar_ZMFtop_V3U.Dot(kbone) * lep_ZMFtbar_V3U.Dot(-1*signyp*nbtwo)
        # calculate polarization observables and add variables to event
        event.b_ntop  = lbar_ZMFtop_V3U.Dot(signyp*nbtwo)
        event.b_ntbar = lep_ZMFtbar_V3U.Dot(-1*signyp*nbtwo)
        event.b_rtop  = lbar_ZMFtop_V3U.Dot(signyp*rbtwo)
        event.b_rtbar = lep_ZMFtbar_V3U.Dot(-1*signyp*rbtwo)
        event.b_ktop  = lbar_ZMFtop_V3U.Dot(kbone)
        event.b_ktbar = lep_ZMFtbar_V3U.Dot(-1*kbone)
        # third base for spin polarization
        signdrap = -1.
        if abs(top_lab.Rapidity())>abs(tbar_lab.Rapidity())  :
            signdrap = 1.
        kbthree = ROOT.TVector3(signdrap * kbone)
        rbthree = ROOT.TVector3(signdrap * signyp * rbtwo)
        event.b_rstartop  = lbar_ZMFtop_V3U.Dot(rbthree)
        event.b_rstartbar = lep_ZMFtbar_V3U.Dot(-1*rbthree)
        event.b_kstartop  = lbar_ZMFtop_V3U.Dot(kbthree)
        event.b_kstartbar = lep_ZMFtbar_V3U.Dot(-1*kbthree)

                
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
    sequence.append( lambda event, sample: makeSpinCorrelationObservables( event, sample, level) )
    sequence.append( lambda event, sample: makeObservables( event, sample, level ) )

    return sequence
    
def getPlotList( scaleLumi, level ):

    tag = 'reco' if level == 'reco' else 'genLep'
    preTag = 'reco' if level == 'reco' else 'gen'

    if scaleLumi: y_label = 'norm. diff. xsec'
    else:         y_label = 'diff. x-sec'

    plots = []
    fisherInfoVariables = []
    
    plots.append(Plot( name = 'l0gammaDPhi',
      texX = '#Delta#Phi(l_{0},#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.l0gammadPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'l0gammaDR',
      texX = '#Delta R(l_{0},#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.l0gammadR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'l1gammaDPhi',
      texX = '#Delta#Phi(l_{1},#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.l1gammadPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'l1gammaDR',
      texX = '#Delta R(l_{1},#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.l1gammadR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = "mll",
      texX = 'm(ll) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.mll if event.passing_checks else float('nan'),
      binning=[50,0,200],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = "mllgamma",
      texX = 'm(ll#gamma) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.mllgamma if event.passing_checks else float('nan'),
      binning=[50,0,200],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = "l0_PdgId",
      texX = 'motherPdgId(l_{0})', texY = y_label,
      attribute = lambda event, sample: event.l0['pdgId'] if event.passing_checks else float('nan'),
      binning=[52,-26,26],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = "l1_PdgId",
      texX = 'motherPdgId(l_{0})', texY = y_label,
      attribute = lambda event, sample: event.l1['pdgId'] if event.passing_checks else float('nan'),
      binning=[52,-26,26],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = "gamma_pt",
      texX = 'p_{T}(#gamma_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.gamma0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    fisherInfoVariables.append('%sPhoton_pt[0]'%preTag)


    plots.append(Plot( name = "gamma_phi",
      texX = '#phi(#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.gamma0['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    fisherInfoVariables.append('%sPhoton_phi[0]'%preTag)


    plots.append(Plot( name = "gamma_eta",
      texX = '#eta(#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.gamma0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    fisherInfoVariables.append('%sPhoton_eta[0]'%preTag)


    plots.append(Plot( name = 'minLeptonDR',
      texX = 'minLeptonDR(#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.gamma0['minLeptonDR'] if event.passing_checks else float('nan'),
      binning=[27,0.3,3],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'minLeptonDR_zoom',
      texX = 'minLeptonDR(#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.gamma0['minLeptonDR'] if event.passing_checks else float('nan'),
      binning=[21,0.3,1],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'minJetDR',
      texX = 'minJetDR(#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.gamma0['minJetDR'] if event.passing_checks else float('nan'),
      binning=[27,0.3,3],
    ))
    fisherInfoVariables.append(None)

    plots.append(Plot( name = 'minJetDR_zoom',
      texX = 'minJetDR(#gamma_{0})', texY = y_label,
      attribute = lambda event, sample: event.gamma0['minJetDR'] if event.passing_checks else float('nan'),
      binning=[21,0.3,1],
    ))
    fisherInfoVariables.append(None)


    plots.append(Plot( name = 'deltaR_lepg0',
      texX = 'min(#DeltaR(lep, #gamma_{0}))', texY = y_label,
      attribute = lambda event, sample: event.minLeptonG0dR if event.passing_checks else float('nan'),
      binning=[20,0.3,5],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_lepg1',
      texX = 'min(#DeltaR(lep, #gamma_{1}))', texY = y_label,
      attribute = lambda event, sample: event.minLeptonG1dR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_lepg2',
      texX = 'min(#DeltaR(lep, #gamma_{2}))', texY = y_label,
      attribute = lambda event, sample: event.minLeptonG2dR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_jetg0',
      texX = 'min(#DeltaR(jet, #gamma_{0}))', texY = y_label,
      attribute = lambda event, sample: event.minJetG0dR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_jetg1',
      texX = 'min(#DeltaR(jet, #gamma_{1}))', texY = y_label,
      attribute = lambda event, sample: event.minJetG1dR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_jetg2',
      texX = 'min(#DeltaR(jet, #gamma_{2}))', texY = y_label,
      attribute = lambda event, sample: event.minJetG2dR if event.passing_checks else float('nan'),
      binning=[20,0.3,3],
    ))
    fisherInfoVariables.append(None)

    
    if level == 'reco':

        plots.append(Plot( name = "gamma_motherPdgId",
          texX = 'motherPdgId(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.genPhoton_motherPdgId[event.recoPhoton_genIndex[0]] if event.recoPhoton_genIndex[0] >= 0 and event.passing_checks else float('nan'),
          binning=[52,-26,26],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma0_isolationVar",
          texX = 'isolationVar(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['isolationVar'] if event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma0_isolationVar_loose",
          texX = 'isolationVar(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['isolationVar'] if event.passing_checks else float('nan'),
          binning=[20,0,0.4],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma0_minLeptonPt",
          texX = 'minLeptonPt(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['minLeptonPt'] if event.passing_checks else float('nan'),
          binning=[20,0,200],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "#gamma0_isolationVarRhoCorr",
#          texX = 'isolationVarRhoCorr(#gamma_{0})', texY = y_label,
#          attribute = lambda event, sample:  event.gamma0['isolationVarRhoCorr'] if event.passing_checks else float('nan'),
#          binning=[20,0,0.05],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma0_sumPtCharged",
          texX = 'sumPtCharged(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample:  event.gamma0['sumPtCharged'] if event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma0_sumPtNeutral",
          texX = 'sumPtNeutral(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample:  event.gamma0['sumPtNeutral'] if event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "gamma0_sumPtChargedPU",
#          texX = 'sumPtChargedPU(#gamma_{0})', texY = y_label,
#          attribute = lambda event, sample:  event.gamma0['sumPtChargedPU'] if event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma0_sumPt",
          texX = 'sumPt(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample:  event.gamma0['sumPt'] if event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

 #       plots.append(Plot( name = "gamma0_ehadOverEem",
 #         texX = 'ehadOverEem(#gamma_{0})', texY = y_label,
 #         attribute = lambda event, sample:  event.gamma0['ehadOverEem'] if event.passing_checks else float('nan'),
 #         binning=[20,0,0.1],
 #       ))
 #       fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma0_genIndex",
          texX = 'genIndex(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample:  event.gamma0['genIndex'] if event.passing_checks else float('nan'),
          binning=[6,-2,4],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma0_IDeltaBeta",
          texX = 'I^{#Delta#beta}_{rel}(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample:  event.gamma0['IDeltaBeta'] if event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_isolationVar_e",
          texX = 'isolationVar(l_{0})', texY = y_label,
          attribute = lambda event, sample: event.l0['isolationVar'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_isolationVar",
          texX = 'isolationVar(l_{0})', texY = y_label,
          attribute = lambda event, sample: event.l0['isolationVar'] if event.passing_checks else float('nan'),
          binning=[20,0,0.5],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_isolationVar",
          texX = 'isolationVar(l_{1})', texY = y_label,
          attribute = lambda event, sample: event.l1['isolationVar'] if event.passing_checks else float('nan'),
          binning=[20,0,0.5],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_isolationVar_e_loose",
          texX = 'isolationVar(l_{0})', texY = y_label,
          attribute = lambda event, sample: event.l0['isolationVar'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,0.4],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l0_isolationVarRhoCorr_e",
#          texX = 'isolationVarRhoCorr(l_{0})', texY = y_label,
#          attribute = lambda event, sample:  event.l0['isolationVarRhoCorr'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_sumPtCharged_e",
          texX = 'sumPtCharged(l_{0})', texY = y_label,
          attribute = lambda event, sample:  event.l0['sumPtCharged'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_sumPtNeutral_e",
          texX = 'sumPtNeutral(l_{0})', texY = y_label,
          attribute = lambda event, sample:  event.l0['sumPtNeutral'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l0_sumPtChargedPU_e",
#          texX = 'sumPtChargedPU(l_{0})', texY = y_label,
#          attribute = lambda event, sample:  event.l0['sumPtChargedPU'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
#          binning=[20,0,2],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_sumPt_e",
          texX = 'sumPt(l_{0})', texY = y_label,
          attribute = lambda event, sample:  event.l0['sumPt'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l0_ehadOverEem_e",
#          texX = 'ehadOverEem(l_{0})', texY = y_label,
#          attribute = lambda event, sample:  event.l0['ehadOverEem'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_IDeltaBeta_e",
          texX = 'I^{#Delta#beta}_{rel}(l_{0})', texY = y_label,
          attribute = lambda event, sample:  event.l0['IDeltaBeta'] if abs( event.l0['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)


        plots.append(Plot( name = "l0_isolationVar_mu",
          texX = 'isolationVar(l_{0})', texY = y_label,
          attribute = lambda event, sample: event.l0['isolationVar'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_isolationVar_mu_loose",
          texX = 'isolationVar(l_{0})', texY = y_label,
          attribute = lambda event, sample: event.l0['isolationVar'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,0.4],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l0_isolationVarRhoCorr_mu",
#          texX = 'isolationVarRhoCorr(l_{0})', texY = y_label,
#          attribute = lambda event, sample:  event.l0['isolationVarRhoCorr'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_sumPtCharged_mu",
          texX = 'sumPtCharged(l_{0})', texY = y_label,
          attribute = lambda event, sample:  event.l0['sumPtCharged'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_sumPtNeutral_mu",
          texX = 'sumPtNeutral(l_{0})', texY = y_label,
          attribute = lambda event, sample:  event.l0['sumPtNeutral'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l0_sumPtChargedPU_mu",
#          texX = 'sumPtChargedPU(l_{0})', texY = y_label,
#          attribute = lambda event, sample:  event.l0['sumPtChargedPU'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
#          binning=[20,0,2],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_sumPt_mu",
          texX = 'sumPt(l_{0})', texY = y_label,
          attribute = lambda event, sample:  event.l0['sumPt'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l0_ehadOverEem_mu",
#          texX = 'ehadOverEem(l_{0})', texY = y_label,
#          attribute = lambda event, sample:  event.l0['ehadOverEem'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_IDeltaBeta_mu",
          texX = 'I^{#Delta#beta}_{rel}(l_{0})', texY = y_label,
          attribute = lambda event, sample:  event.l0['IDeltaBeta'] if abs( event.l0['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)


        plots.append(Plot( name = "l1_isolationVar_e",
          texX = 'isolationVar(l_{1})', texY = y_label,
          attribute = lambda event, sample: event.l1['isolationVar'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_isolationVar_e_loose",
          texX = 'isolationVar(l_{1})', texY = y_label,
          attribute = lambda event, sample: event.l1['isolationVar'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,0.4],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l1_isolationVarRhoCorr_e",
#          texX = 'isolationVarRhoCorr(l_{1})', texY = y_label,
#          attribute = lambda event, sample:  event.l1['isolationVarRhoCorr'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_sumPtCharged_e",
          texX = 'sumPtCharged(l_{1})', texY = y_label,
          attribute = lambda event, sample:  event.l1['sumPtCharged'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_sumPtNeutral_e",
          texX = 'sumPtNeutral(l_{1})', texY = y_label,
          attribute = lambda event, sample:  event.l1['sumPtNeutral'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l1_sumPtChargedPU_e",
#          texX = 'sumPtChargedPU(l_{1})', texY = y_label,
#          attribute = lambda event, sample:  event.l1['sumPtChargedPU'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
#          binning=[20,0,2],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_sumPt_e",
          texX = 'sumPt(l_{1})', texY = y_label,
          attribute = lambda event, sample:  event.l1['sumPt'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l1_ehadOverEem_e",
#          texX = 'ehadOverEem(l_{1})', texY = y_label,
#          attribute = lambda event, sample:  event.l1['ehadOverEem'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_IDeltaBeta_e",
          texX = 'I^{#Delta#beta}_{rel}(l_{1})', texY = y_label,
          attribute = lambda event, sample:  event.l1['IDeltaBeta'] if abs( event.l1['pdgId'] ) == 11 and event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)


        plots.append(Plot( name = "l1_isolationVar_mu",
          texX = 'isolationVar(l_{1})', texY = y_label,
          attribute = lambda event, sample: event.l1['isolationVar'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_isolationVar_mu_loose",
          texX = 'isolationVar(l_{1})', texY = y_label,
          attribute = lambda event, sample: event.l1['isolationVar'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,0.4],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l1_isolationVarRhoCorr_mu",
#          texX = 'isolationVarRhoCorr(l_{1})', texY = y_label,
#          attribute = lambda event, sample:  event.l1['isolationVarRhoCorr'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_sumPtCharged_mu",
          texX = 'sumPtCharged(l_{1})', texY = y_label,
          attribute = lambda event, sample:  event.l1['sumPtCharged'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_sumPtNeutral_mu",
          texX = 'sumPtNeutral(l_{1})', texY = y_label,
          attribute = lambda event, sample:  event.l1['sumPtNeutral'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l1_sumPtChargedPU_mu",
#          texX = 'sumPtChargedPU(l_{1})', texY = y_label,
#          attribute = lambda event, sample:  event.l1['sumPtChargedPU'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
#          binning=[20,0,2],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_sumPt_mu",
          texX = 'sumPt(l_{1})', texY = y_label,
          attribute = lambda event, sample:  event.l1['sumPt'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,50],
        ))
        fisherInfoVariables.append(None)

#        plots.append(Plot( name = "l1_ehadOverEem_mu",
#          texX = 'ehadOverEem(l_{1})', texY = y_label,
#          attribute = lambda event, sample:  event.l1['ehadOverEem'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
#          binning=[20,0,0.1],
#        ))
#        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_IDeltaBeta_mu",
          texX = 'I^{#Delta#beta}_{rel}(l_{1})', texY = y_label,
          attribute = lambda event, sample:  event.l1['IDeltaBeta'] if abs( event.l1['pdgId'] ) == 13 and event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))
        fisherInfoVariables.append(None)

    elif level == 'gen':

        plots.append(Plot( name = "gamma_relIso04_l",
          texX = 'relIso04(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['relIso04'] if abs(event.gamma0['motherPdgId']) in [11,13] and event.passing_checks else float('nan'),
          binning=[20,0,0.5],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_relIso04_q",
          texX = 'relIso04(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['relIso04'] if abs(event.gamma0['motherPdgId']) < 5 and event.passing_checks else float('nan'),
          binning=[20,0,0.5],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_relIso04_g",
          texX = 'relIso04(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['relIso04'] if abs(event.gamma0['motherPdgId']) == 21 and event.passing_checks else float('nan'),
          binning=[20,0,0.5],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_relIso04_all",
          texX = 'relIso04(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['relIso04'] if event.passing_checks else float('nan'),
          binning=[20,0,0.5],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_relIso04_q_zoom",
          texX = 'relIso04(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['relIso04'] if abs(event.gamma0['motherPdgId']) < 5 and event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_relIso04_l_zoom",
          texX = 'relIso04(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['relIso04'] if abs(event.gamma0['motherPdgId']) in [11,13] and event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_relIso04_g_zoom",
          texX = 'relIso04(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['relIso04'] if abs(event.gamma0['motherPdgId']) == 21 and event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_relIso04_all_zoom",
          texX = 'relIso04(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['relIso04'] if event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_motherPdgId",
          texX = 'motherPdgId(#gamma_{0})', texY = y_label,
          attribute = lambda event, sample: event.gamma0['motherPdgId'] if event.passing_checks else float('nan'),
          binning=[52,-26,26],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "gamma_mass",
          texX = 'm(#gamma_{0}) [GeV]', texY = y_label,
          attribute = lambda event, sample: event.gamma0['mass'] if event.passing_checks else float('nan'),
          binning=[20,-1e-4,1e-4],
        ))
        fisherInfoVariables.append(None)


        plots.append(Plot( name = "l0_pt_gamma0FromE",
          texX = 'p_{T}(l_{0}) [GeV] if #gamma_{0} origin = e^{-},e^{+}', texY = y_label,
          attribute = lambda event, sample: event.l0['pt'] if abs(event.gamma0['motherPdgId'])==11 and event.passing_checks else float('nan'),
          binning=[20,0,100],
        ))
        fisherInfoVariables.append(None)
    
    
        plots.append(Plot( name = "l0_pt_gamma0FromMu",
          texX = 'p_{T}(l_{0}) [GeV] if #gamma_{0} origin = #mu^{-},#mu^{+}', texY = y_label,
          attribute = lambda event, sample: event.l0['pt'] if abs(event.gamma0['motherPdgId'])==13 and event.passing_checks else float('nan'),
          binning=[20,0,100],
        ))
        fisherInfoVariables.append(None)
    
    
        plots.append(Plot( name = "l0_pt_gamma0FromW",
          texX = 'p_{T}(l_{0}) [GeV] if #gamma_{0} origin = W^{-},W^{+}', texY = y_label,
          attribute = lambda event, sample: event.l0['pt'] if abs(event.gamma0['motherPdgId'])==24 and event.passing_checks else float('nan'),
          binning=[20,0,100],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l0_motherPdgId",
          texX = 'motherPdgId(l_{0})', texY = y_label,
          attribute = lambda event, sample: event.l0['motherPdgId'] if event.passing_checks else float('nan'),
          binning=[52,-26,26],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "l1_motherPdgId",
          texX = 'motherPdgId(l_{1})', texY = y_label,
          attribute = lambda event, sample: event.l1['motherPdgId'] if event.passing_checks else float('nan'),
          binning=[52,-26,26],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_nn",
          texX = 'C_{nn}', texY = y_label,
          attribute = lambda event, sample: event.c_nn if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_rr",
          texX = 'C_{rr}', texY = y_label,
          attribute = lambda event, sample: event.c_rr if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_kk",
          texX = 'C_{kk}', texY = y_label,
          attribute = lambda event, sample: event.c_kk if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_rkpkr",
          texX = 'C_{rk}+C_{kr}', texY = y_label,
          attribute = lambda event, sample: event.c_rkpkr if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_nrprn",
          texX = 'C_{nr}+C_{rn}', texY = y_label,
          attribute = lambda event, sample: event.c_nrprn if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_nkpkn",
          texX = 'C_{nk}+C_{kn}', texY = y_label,
          attribute = lambda event, sample: event.c_nkpkn if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_rkmkr",
          texX = 'C_{rk}-C_{kr}', texY = y_label,
          attribute = lambda event, sample: event.c_rkmkr if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_nrmrn",
          texX = 'C_{nr}-C_{rn}', texY = y_label,
          attribute = lambda event, sample: event.c_nrmrn if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "c_nkmkn",
          texX = 'C_{nk}-C_{kn}', texY = y_label,
          attribute = lambda event, sample: event.c_nkmkn if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_ntop",
          texX = 'b_{n}(t)', texY = y_label,
          attribute = lambda event, sample: event.b_ntop if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_ntbar",
          texX = 'b_{n}(#overline{t})', texY = y_label,
          attribute = lambda event, sample: event.b_ntbar if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_rtop",
          texX = 'b_{r}(t)', texY = y_label,
          attribute = lambda event, sample: event.b_rtop if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_rtbar",
          texX = 'b_{r}(#overline{t})', texY = y_label,
          attribute = lambda event, sample: event.b_rtbar if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_ktop",
          texX = 'b_{k}(t)', texY = y_label,
          attribute = lambda event, sample: event.b_ktop if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_ktbar",
          texX = 'b_{k}(#overline{t})', texY = y_label,
          attribute = lambda event, sample: event.b_ktbar if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_rstartop",
          texX = 'b_{r^*}(t)', texY = y_label,
          attribute = lambda event, sample: event.b_rstartop if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_rstartbar",
          texX = 'b_{r^*}(#overline{t})', texY = y_label,
          attribute = lambda event, sample: event.b_rstartbar if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_kstartop",
          texX = 'b_{k^*}(t)', texY = y_label,
          attribute = lambda event, sample: event.b_kstartop if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

        plots.append(Plot( name = "b_kstartbar",
          texX = 'b_{k^*}(#overline{t})', texY = y_label,
          attribute = lambda event, sample: event.b_kstartbar if event.passing_checks else float('nan'),
          binning=[20,-1,1],
        ))
        fisherInfoVariables.append(None)

    plots.append(Plot( name = "b0_pt",
      texX = 'p_{T}(b_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.bj0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = "b0_phi",
      texX = '#phi(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = "b0_eta",
      texX = '#eta(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = "b1_pt",
      texX = 'p_{T}(b_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.bj1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = "b1_phi",
      texX = '#phi(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = "b1_eta",
      texX = '#eta(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = "l0_pt",
      texX = 'p_{T}(l_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.l0['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,300],
    ))
    fisherInfoVariables.append('%sLep_pt[0]'%preTag)

    
    plots.append(Plot( name = "l0_pt_zoom",
      texX = 'p_{T}(l_{0}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.l0['pt'] if event.passing_checks else float('nan'),
      binning=[50,15,65],
    ))
    fisherInfoVariables.append('%sLep_pt[0]'%preTag)

    
    plots.append(Plot( name = "l0_phi",
      texX = '#phi(l_{0})', texY = y_label,
      attribute = lambda event, sample: event.l0['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    fisherInfoVariables.append('%sLep_phi[0]'%preTag)

    
    plots.append(Plot( name = "l0_eta",
      texX = '#eta(l_{0})', texY = y_label,
      attribute = lambda event, sample: event.l0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    fisherInfoVariables.append('%sLep_eta[0]'%preTag)


    plots.append(Plot( name = "l1_pt",
      texX = 'p_{T}(l_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.l1['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,300],
    ))
    fisherInfoVariables.append('%sLep_pt[1]'%preTag)

    plots.append(Plot( name = "l1_pt_zoom",
      texX = 'p_{T}(l_{1}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.l1['pt'] if event.passing_checks else float('nan'),
      binning=[50,15,65],
    ))
    fisherInfoVariables.append('%sLep_pt[1]'%preTag)

    
    plots.append(Plot( name = "l1_phi",
      texX = '#phi(l_{1})', texY = y_label,
      attribute = lambda event, sample: event.l1['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    fisherInfoVariables.append('%sLep_phi[1]'%preTag)

    
    plots.append(Plot( name = "l1_eta",
      texX = '#eta(l_{1})', texY = y_label,
      attribute = lambda event, sample: event.l1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    fisherInfoVariables.append('%sLep_eta[1]'%preTag)

    
    plots.append(Plot( name = 'deltaPhi_bb',
      texX = '#Delta#phi(bb)', texY = y_label,
      attribute = lambda event, sample: event.bbdPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_bb',
      texX = '#DeltaR(bb)', texY = y_label,
      attribute = lambda event, sample: event.bbdR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaPhi_ll',
      texX = '#Delta#phi(ll)', texY = y_label,
      attribute = lambda event, sample: event.lldPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_ll',
      texX = '#DeltaR(ll)', texY = y_label,
      attribute = lambda event, sample: event.lldR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    fisherInfoVariables.append(None)


    plots.append(Plot( name = 'deltaR_lepg0',
      texX = 'min(#DeltaR(lep, #gamma_{0}))', texY = y_label,
      attribute = lambda event, sample: event.minLeptonG0dR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_lepg1',
      texX = 'min(#DeltaR(lep, #gamma_{1}))', texY = y_label,
      attribute = lambda event, sample: event.minLeptonG1dR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_lepg2',
      texX = 'min(#DeltaR(lep, #gamma_{2}))', texY = y_label,
      attribute = lambda event, sample: event.minLeptonG2dR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_jetg0',
      texX = 'min(#DeltaR(jet, #gamma_{0}))', texY = y_label,
      attribute = lambda event, sample: event.minJetG0dR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_jetg1',
      texX = 'min(#DeltaR(jet, #gamma_{1}))', texY = y_label,
      attribute = lambda event, sample: event.minJetG1dR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'deltaR_jetg2',
      texX = 'min(#DeltaR(jet, #gamma_{2}))', texY = y_label,
      attribute = lambda event, sample: event.minJetG2dR if event.passing_checks else float('nan'),
      binning=[20,0,6],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = '2nu_Met_pt',
      texX = 'E_{T}^{miss} [GeV]', texY = y_label,
      attribute = lambda event, sample: event.MET['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    fisherInfoVariables.append('%sMet_pt'%preTag)

    
    plots.append(Plot( name	= '2nu_Met_phi',
      texX = '#phi(E_{T}^{miss})', texY = y_label,
      attribute = lambda event, sample: event.MET['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    fisherInfoVariables.append('%sMet_phi'%preTag)

    
    plots.append(Plot( name = 'l0_pt_charge',
      texX = 'p_{T}(l_{0}) [GeV] signed with lepton charge', texY = y_label,
      attribute = lambda event, sample: event.lep0chargept if event.passing_checks else float('nan'),
      binning=[20,-200,200],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'l1_pt_charge',
      texX = 'p_{T}(l_{1}) [GeV] signed with lepton charge', texY = y_label,
      attribute = lambda event, sample: event.lep1chargept if event.passing_checks else float('nan'),
      binning=[20,-200,200],
    ))
    fisherInfoVariables.append(None)

    
    plots.append(Plot( name = 'njets',
      texX = 'Number of Jets', texY = y_label,
      attribute = lambda event, sample: getattr( event, 'n%sJet'%preTag ) if event.passing_checks else float('nan'),
      binning=[10,0,10],
    ))
    fisherInfoVariables.append('n%sJet'%preTag)


    plots.append(Plot( name = 'nleps',
      texX = 'Number of Leptons', texY = y_label,
      attribute = lambda event, sample: getattr( event, 'n%sLep'%preTag ) if event.passing_checks else float('nan'),
      binning=[8,0,8],
    ))
    fisherInfoVariables.append('n%sLep'%preTag)


    plots.append(Plot( name = 'nbjets',
      texX = 'Number of bJets', texY = y_label,
      attribute = lambda event, sample: len(event.bjets) if event.passing_checks else float('nan'),
      binning=[4,0,4],
    ))
    fisherInfoVariables.append(None)


    plots.append(Plot( name = 'ngammas',
      texX = 'Number of photons', texY = y_label,
      attribute = lambda event, sample: len(event.gammas) if event.passing_checks else float('nan'),
      binning=[10,0,10],
    ))
    fisherInfoVariables.append(None)


    return plots, fisherInfoVariables

