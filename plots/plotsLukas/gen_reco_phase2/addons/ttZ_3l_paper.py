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
        "ngenLep/I", #"genLep[pt/F,eta/F,phi/F,pdgId/I]",
    
        "genZ_pt/F", "genZ_eta/F", "genZ_phi/F", "genZ_mass/F", "genZ_cosThetaStar/F",
    
        "genBj0_pt/F", "genBj0_phi/F", "genBj0_eta/F", "genBj0_nNeutrals/I", "genBj0_nCharged/I",
        "genBj1_pt/F", "genBj1_phi/F", "genBj1_eta/F", "genBj0_nNeutrals/I", "genBj0_nCharged/I",
    
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
        read_variables_gen.append("recoJet[pt/F,eta/F,phi/F,bTag_medium/I,nNeutrals/I,nCharged/I]")
        read_variables_gen.append("recoLep[pt/F,eta/F,phi/F,pdgId/I,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F]")
    else:
        read_variables_gen.append("genJet[pt/F,eta/F,phi/F,matchBParton/I]")
        read_variables_gen.append("genLep[pt/F,phi/F,eta/F,pdgId/I]")

    read_variables = read_variables_gen + read_variables_genLep + ['nBTag/I']
    read_variables = list( set( read_variables ) ) # remove double entries
#    read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

    return read_variables


def makeJets( event, sample, level ):
    ''' Add a list of filtered jets to the event (full list is required for lepton cross cleaning)
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoJet     as isGoodJet
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenJet      as isGoodJet

    # load jets
    btag = 'bTag_medium' if level == 'reco' else 'matchBParton'
    event.jets = getCollection( event, '%sJet'%preTag, ['pt', 'eta', 'phi', btag, 'nNeutrals', 'nCharged' ], 'n%sJet'%preTag )
#    event.jets = list( filter( lambda j: isGoodJet(j), event.jets ) )

    # Define leptonic b-jets
    event.bj0lep = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi', 'nNeutrals', 'nCharged'], getattr( event, '%sBjLeadlep_index'%preTag ) ) if getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 else nanJet()
    event.bj0had = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi', 'nNeutrals', 'nCharged'], getattr( event, '%sBjLeadhad_index'%preTag ) ) if getattr( event, '%sBjLeadhad_index'%preTag ) >= 0 else nanJet()

    # Define non-Z b-jets
    event.bjNonZlep = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi', 'nNeutrals', 'nCharged'], getattr( event, '%sBjNonZlep_index'%preTag ) ) if getattr( event, '%sBjNonZlep_index'%preTag ) >= 0 else nanJet()
    event.bjNonZhad = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi', 'nNeutrals', 'nCharged'], getattr( event, '%sBjNonZhad_index'%preTag ) ) if getattr( event, '%sBjNonZhad_index'%preTag ) >= 0 else nanJet()

    # get (second) hardest bjets
    event.bj0 = {'pt':getattr( event, '%sBj0_pt'%preTag ), 'phi':getattr( event, '%sBj0_phi'%preTag ), 'eta':getattr( event, '%sBj0_eta'%preTag ), 'nNeutrals':getattr( event, '%sBj0_nNeutrals'%preTag ), 'nCharged':getattr( event, '%sBj0_nCharged'%preTag )}
    event.bj1 = {'pt':getattr( event, '%sBj1_pt'%preTag ), 'phi':getattr( event, '%sBj1_phi'%preTag ), 'eta':getattr( event, '%sBj1_eta'%preTag ),  'nNeutrals':getattr( event, '%sBj0_nNeutrals'%preTag ), 'nCharged':getattr( event, '%sBj0_nCharged'%preTag )}

    # Add extra vectors
    for p in [event.bj0, event.bj1, event.bj0lep, event.bj0had, event.bjNonZlep, event.bjNonZhad]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # selection checks
#    event.foundBj0       = isGoodJet( event.bj0 )

    # choose your selection on b-jets
#    event.passing_bjets = event.foundBj0
#    event.passing_jets  = len( list( filter( lambda j: abs(j['eta'])>2.4, event.jets ) ) ) >= 1


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


def makeLeps( event, sample, level, leptonFlavor ):
    ''' Add important leptons (no full list of leptons is required for now)
    '''
    preTag = 'reco' if level == 'reco' else 'gen'
    tag    = 'reco' if level == 'reco' else 'genLep'

    if level == 'reco':
        leptonList = ['pt', 'eta', 'phi', 'pdgId', 'isolationVar', 'isolationVarRhoCorr', 'sumPtCharged', 'sumPtNeutral', 'sumPtChargedPU', 'sumPt', 'ehadOverEem', 'genIndex']
    else:
        leptonList = ['pt', 'eta', 'phi', 'pdgId']

    # Define Z leptons
    event.Z_l0 = getObjDict( event, '%sLep_'%preTag, leptonList, getattr( event, '%sZ_l1_index'%tag     ) ) if getattr( event, '%sZ_l1_index'%tag ) >= 0 else nanLepton()
    event.Z_l1 = getObjDict( event, '%sLep_'%preTag, leptonList, getattr( event, '%sZ_l2_index'%tag     ) ) if getattr( event, '%sZ_l2_index'%tag ) >= 0 else nanLepton()

    # Define non-Z leptons
    event.NonZ_l0 = getObjDict( event, '%sLep_'%preTag, leptonList, getattr( event, '%sNonZ_l1_index'%tag     ) ) if getattr( event, '%sNonZ_l1_index'%tag ) >= 0 else nanLepton()
    event.NonZ_l1 = getObjDict( event, '%sLep_'%preTag, leptonList, getattr( event, '%sNonZ_l2_index'%tag     ) ) if getattr( event, '%sNonZ_l2_index'%tag ) >= 0 else nanLepton()

    # Add extra vectors
    for p in [ event.Z_l0, event.Z_l1, event.NonZ_l0, event.NonZ_l1 ]:
        addTransverseVector( p )
        addTLorentzVector( p )
        if level == 'reco':
            addIDeltaBeta( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoLepton    as isGoodLepton
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenLepton   as isGoodLepton

    # We may loose some events by cross-cleaning or by thresholds.
    event.foundZl0     = getattr( event, '%sZ_l1_index'%tag ) >= 0 #and isGoodLepton( event.Z_l0 )
    event.foundZl1     = getattr( event, '%sZ_l2_index'%tag ) >= 0 #and isGoodLepton( event.Z_l1 )
    event.foundZ       = event.Z_l0['pdgId'] * event.Z_l1['pdgId'] < 0 and abs(event.Z_l0['pdgId']) == abs(event.Z_l1['pdgId'])
    event.found3lep    = getattr( event, '%sNonZ_l1_index'%tag ) >= 0 #and isGoodLepton( event.NonZ_l0 )

    # choose your selection on leptons
    event.passing_leptons = event.found3lep and event.foundZl0 and event.foundZl1 and event.foundZ

    lepFlavors = [ int(abs(event.NonZ_l0['pdgId'])), int(abs(event.Z_l0['pdgId'])), int(abs(event.Z_l1['pdgId'])) ]

    if leptonFlavor == 'e': event.passing_leptons = event.passing_leptons and abs(event.NonZ_l0['pdgId'])==11
    elif leptonFlavor == 'mu': event.passing_leptons = event.passing_leptons and abs(event.NonZ_l0['pdgId'])==13
    elif leptonFlavor == 'mumumu': event.passing_leptons = event.passing_leptons and lepFlavors.count(13)==3
    elif leptonFlavor == 'mumue': event.passing_leptons = event.passing_leptons and lepFlavors.count(13)==2 and lepFlavors.count(11)==1
    elif leptonFlavor == 'muee': event.passing_leptons = event.passing_leptons and lepFlavors.count(13)==1 and lepFlavors.count(11)==2
    elif leptonFlavor == 'eee': event.passing_leptons = event.passing_leptons and lepFlavors.count(11)==3


def makeObservables( event, sample, level):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.bbdPhi = deltaPhi( event.bj0['phi'], event.bj1['phi'] )
    event.bbdR   = deltaR( event.bj0, event.bj1 )

    # Make leptonic W, Nan if any component is Nan
    event.Wlep_vec2D = event.MET['vec2D'] + event.NonZ_l0['vec2D']
    # Lp, Nan if any component is Nan
    event.Wlep_Lp = ( event.Wlep_vec2D*event.NonZ_l0['vec2D'] ) / ( event.Wlep_vec2D*event.Wlep_vec2D )
    # classic MT, Nan if any component is Nan
    event.Wlep_MT = sqrt( MTSquared( event.MET, event.NonZ_l0     ) )

    # blep+lep subsystem, Nan if len(event.lepsNotFromZ == 0)
    event.bleplep_vec2D = event.bjNonZlep['vec2D'] + event.NonZ_l0['vec2D']
    event.bleplep_vec4D = event.bjNonZlep['vec4D'] + event.NonZ_l0['vec4D']

    # transverse mass of top, Nan if len(event.lepsNotFromZ == 0)
    event.t_MT = sqrt( MTSquared( event.NonZ_l0, event.MET ) + MSquared( event.NonZ_l0, event.bjNonZlep ) + MTSquared( event.MET, event.bjNonZlep     ) )
    event.t_vec2D = event.NonZ_l0['vec2D'] + event.MET['vec2D'] + event.bjNonZlep['vec2D']

    # signed lepton pt, Nan if len(event.lepsNotFromZ == 0)
    event.getnonZlepchargept = event.NonZ_l0['pt'] if event.NonZ_l0['pdgId']>0 else -event.NonZ_l0['pt']

    # choose your final selection
    event.passing_checks = event.passing_leptons #and event.recoZ_lldPhi<1 and event.recoZ_lldR<1.5 #and event.passing_bjets


def getSequenceList( level, leptonFlavor ):
    ''' sequence functions
    '''
    sequence = []

    sequence.append( lambda event, sample: makeJets( event, sample, level ) )
    sequence.append( lambda event, sample: makeMET( event, sample, level ) )
    sequence.append( lambda event, sample: makeZ( event, sample, level ) )
    sequence.append( lambda event, sample: makeLeps( event, sample, level, leptonFlavor ) )
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

    y_label = 'Number of Events'

    plots = []
    fisherInfoVariables = []
    

    plots.append( Plot( name = "nonZl0_isolationVar",
      texX = '#isolationVar(l_{0}^{non-Z})', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l0['isolationVar'] if event.passing_checks else float('nan'),
      binning=[50,0,0.2],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Zl0_isolationVar",
      texX = '#isolationVar(l_{0}^{Z})', texY = y_label,
      attribute = lambda event, sample: event.Z_l0['isolationVar'] if event.passing_checks else float('nan'),
      binning=[50,0,0.2],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Zl1_isolationVar",
      texX = '#isolationVar(l_{1}^{Z})', texY = y_label,
      attribute = lambda event, sample: event.Z_l1['isolationVar'] if event.passing_checks else float('nan'),
      binning=[50,0,0.2],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "nonZl0_sumPtCharged",
      texX = '#sumPtCharged(l_{0}^{non-Z})', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l0['sumPtCharged'] if event.passing_checks else float('nan'),
      binning=[50,0,10],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Zl0_sumPtCharged",
      texX = '#sumPtCharged(l_{0}^{Z})', texY = y_label,
      attribute = lambda event, sample: event.Z_l0['sumPtCharged'] if event.passing_checks else float('nan'),
      binning=[50,0,10],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Zl1_sumPtCharged",
      texX = '#sumPtCharged(l_{1}^{Z})', texY = y_label,
      attribute = lambda event, sample: event.Z_l1['sumPtCharged'] if event.passing_checks else float('nan'),
      binning=[50,0,10],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "nonZl0_sumPtNeutral",
      texX = '#sumPtNeutral(l_{0}^{non-Z})', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l0['sumPtNeutral'] if event.passing_checks else float('nan'),
      binning=[50,0,10],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Zl0_sumPtNeutral",
      texX = '#sumPtNeutral(l_{0}^{Z})', texY = y_label,
      attribute = lambda event, sample: event.Z_l0['sumPtNeutral'] if event.passing_checks else float('nan'),
      binning=[50,0,10],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Zl1_sumPtNeutral",
      texX = '#sumPtNeutral(l_{1}^{Z})', texY = y_label,
      attribute = lambda event, sample: event.Z_l1['sumPtNeutral'] if event.passing_checks else float('nan'),
      binning=[50,0,10],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "nonZl0_eta",
      texX = '#eta(l_{0}^{non-Z})', texY = y_label,
      attribute = lambda event, sample: event.NonZ_l0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Zl0_eta",
      texX = '#eta(l_{0}^{Z})', texY = y_label,
      attribute = lambda event, sample: event.Z_l0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Zl1_eta",
      texX = '#eta(l_{1}^{Z})', texY = y_label,
      attribute = lambda event, sample: event.Z_l1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "l0_pt20",
      texX = 'p_{T}(l_{0}) [GeV]', texY = y_label + ' / 20 GeV',
      attribute = lambda event, sample: getattr( event, '%sLep_pt[0]'%level ) if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "l1_pt20",
      texX = 'p_{T}(l_{1}) [GeV]', texY = y_label + ' / 20 GeV',
      attribute = lambda event, sample: getattr( event, '%sLep_pt[1]'%level ) if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "l2_pt20",
      texX = 'p_{T}(l_{2}) [GeV]', texY = y_label + ' / 20 GeV',
      attribute = lambda event, sample: getattr( event, '%sLep_pt[2]'%level ) if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "l0_pt10",
      texX = 'p_{T}(l_{0}) [GeV]', texY = y_label + ' / 10 GeV',
      attribute = lambda event, sample: getattr( event, '%sLep_pt[0]'%level ) if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "l1_pt10",
      texX = 'p_{T}(l_{1}) [GeV]', texY = y_label + ' / 10 GeV',
      attribute = lambda event, sample: getattr( event, '%sLep_pt[1]'%level ) if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "l2_pt10",
      texX = 'p_{T}(l_{2}) [GeV]', texY = y_label + ' / 10 GeV',
      attribute = lambda event, sample: getattr( event, '%sLep_pt[2]'%level ) if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "jet0_pt",
      texX = 'p_{T}(jet_{0}) [GeV]', texY = y_label + ' / 20 GeV',
      attribute = lambda event, sample: event.jets[0]['pt'] if event.passing_checks and len(event.jets)>0 else float('nan'),
      binning=[20,0,400],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "jet1_pt",
      texX = 'p_{T}(jet_{1}) [GeV]', texY = y_label + ' / 20 GeV',
      attribute = lambda event, sample: event.jets[1]['pt'] if event.passing_checks and len(event.jets)>1 else float('nan'),
      binning=[20,0,400],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "jet0_eta",
      texX = '#eta(j_{0})', texY = y_label,
      attribute = lambda event, sample: event.jets[0]['eta'] if event.passing_checks and len(event.jets)>0 else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "jet1_eta",
      texX = '#eta(j_{1})', texY = y_label,
      attribute = lambda event, sample: event.jets[1]['eta'] if event.passing_checks and len(event.jets)>1 else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "jet2_eta",
      texX = '#eta(j_{2})', texY = y_label,
      attribute = lambda event, sample: event.jets[2]['eta'] if event.passing_checks and len(event.jets)>2 else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)


    plots.append( Plot( name = "b0_eta",
      texX = '#eta(b_{0})', texY = y_label,
      attribute = lambda event, sample: event.bj0['eta'] if event.passing_checks else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "b1_eta",
      texX = '#eta(b_{1})', texY = y_label,
      attribute = lambda event, sample: event.bj1['eta'] if event.passing_checks else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "Z_eta",
      texX = '#eta(Z)', texY = y_label,
      attribute = lambda event, sample: getattr( event, '%sZ_eta'%level ) if event.passing_checks else float('nan'),
      binning=[20,-4,4],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "jet0_nCharged",
      texX = 'N_{charged}(jet_{0})', texY = y_label,
      attribute = lambda event, sample: event.jets[0]['nCharged'] if event.passing_checks and len(event.jets)>0 else float('nan'),
      binning=[20,0,20],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "jet1_nCharged",
      texX = 'N_{charged}(jet_{1})', texY = y_label,
      attribute = lambda event, sample: event.jets[1]['nCharged'] if event.passing_checks and len(event.jets)>1 else float('nan'),
      binning=[20,0,20],
    ) )
    fisherInfoVariables.append(None)


    plots.append( Plot( name = "jet0_nNeutrals",
      texX = 'N_{neutral}(jet_{0})', texY = y_label,
      attribute = lambda event, sample: event.jets[0]['nNeutrals'] if event.passing_checks and len(event.jets)>0 else float('nan'),
      binning=[20,0,20],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = "jet1_nNeutrals",
      texX = 'N_{neutral}(jet_{1})', texY = y_label,
      attribute = lambda event, sample: event.jets[1]['nNeutrals'] if event.passing_checks and len(event.jets)>1 else float('nan'),
      binning=[20,0,20],
    ) )
    fisherInfoVariables.append(None)


    plots.append( Plot( name = "Z_pt20",
      texX = 'p_{T}(ll) [GeV]', texY = y_label + ' / 20 GeV',
#      texX = 'p_{T}(ll) [GeV]', texY = labelAddon + 'd#sigma/dp_{T}(ll)'+unit,
#      texX = 'p_{T}(ll) [GeV]', texY = labelAddon + '#frac{d#sigma}{d p_{T}(ll)}'+unit,
      attribute = lambda event, sample: getattr( event, '%sZ_pt'%level ) if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ) )
    fisherInfoVariables.append('%sZ_pt'%level)

    plots.append( Plot( name = "Z_pt10",
      texX = 'p_{T}(ll) [GeV]', texY = y_label + ' / 40 GeV',
#      texX = 'p_{T}(ll) [GeV]', texY = labelAddon + 'd#sigma/dp_{T}(ll)'+unit,
#      texX = 'p_{T}(ll) [GeV]', texY = labelAddon + '#frac{d#sigma}{d p_{T}(ll)}'+unit,
      attribute = lambda event, sample: getattr( event, '%sZ_pt'%level ) if event.passing_checks else float('nan'),
      binning=[10,0,400],
    ) )
    fisherInfoVariables.append('%sZ_pt'%level)

    # for gen use genLepZ_lldPhi
    plots.append( Plot( name = 'Z_deltaPhi_ll',
      texX = '#Delta#phi(ll)', texY = y_label,
#      texX = '#Delta#phi(ll)', texY = labelAddon + 'd#sigma/d#Delta#phi(ll)'+unit,
#      texX = '#Delta#phi(ll)', texY = labelAddon + '#frac{d#sigma}{d #Delta#phi(ll)}'+unit,
      attribute = lambda event, sample: getattr( event, '%sZ_lldPhi'%tag ) if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ) )
    fisherInfoVariables.append('%sZ_lldPhi'%tag)

    # for gen use genLepZ_lldR
    plots.append( Plot( name = 'Z_deltaR_ll',
      texX = '#DeltaR(ll)', texY = y_label,
#      texX = '#DeltaR(ll)', texY = labelAddon + 'd#sigma/d#DeltaR(ll)'+unit,
#      texX = '#DeltaR(ll)', texY = labelAddon + '#frac{d#sigma}{d #DeltaR(ll)}'+unit,
      attribute = lambda event, sample: getattr( event, '%sZ_lldR'%tag ) if event.passing_checks else float('nan'),
      binning=[20,0,3.5],
    ) )
    fisherInfoVariables.append('%sZ_lldR'%tag)
        
    plots.append( Plot( name = "Z_cosThetaStar10",
      texX = 'cos(#theta_{Z}*)', texY = y_label,
#      texX = 'cos(#theta*)', texY = labelAddon + 'd#sigma/dcos(#theta*)'+unit,
#      texX = 'cos(#theta*)', texY = labelAddon + '#frac{d#sigma}{d cos(#theta*)}'+unit,
      attribute = lambda event, sample: getattr( event, '%sZ_cosThetaStar'%level ) if event.passing_checks else float('nan'),
      binning=[10,-1.,1.],
    ) )
    fisherInfoVariables.append('%sZ_cosThetaStar'%level)

    plots.append( Plot( name = "Z_cosThetaStar20",
      texX = 'cos(#theta_{Z}*)', texY = y_label,
#      texX = 'cos(#theta*)', texY = labelAddon + 'd#sigma/dcos(#theta*)'+unit,
#      texX = 'cos(#theta*)', texY = labelAddon + '#frac{d#sigma}{d cos(#theta*)}'+unit,
      attribute = lambda event, sample: getattr( event, '%sZ_cosThetaStar'%level ) if event.passing_checks else float('nan'),
      binning=[20,-1.,1.],
    ) )
    fisherInfoVariables.append('%sZ_cosThetaStar'%level)

    plots.append( Plot( name = "Z_cosThetaStar5",
      texX = 'cos(#theta*_{Z})', texY = y_label,
#      texX = 'cos(#theta*)', texY = labelAddon + 'd#sigma/dcos(#theta*)'+unit,
#      texX = 'cos(#theta*)', texY = labelAddon + '#frac{d#sigma}{d cos(#theta*)}'+unit,
      attribute = lambda event, sample: getattr( event, '%sZ_cosThetaStar'%level ) if event.passing_checks else float('nan'),
      binning=[5,-1.,1.],
    ) )
    fisherInfoVariables.append('%sZ_cosThetaStar'%level)
    
    plots.append( Plot( name = 'deltaPhi_bb',
      texX = '#Delta#phi(bb)', texY = y_label,
#      texX = '#Delta#phi(bb)', texY = labelAddon + 'd#sigma/d#Delta#phi(bb)'+unit,
#      texX = '#Delta#phi(bb)', texY = labelAddon + '#frac{d#sigma}{d #Delta#phi(bb)}'+unit,
      attribute = lambda event, sample: event.bbdPhi if event.passing_checks else float('nan'),
      binning=[20,0,pi],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = 'deltaR_bb',
      texX = '#DeltaR(bb)', texY = y_label,
#      texX = '#DeltaR(bb)', texY = labelAddon + 'd#sigma/d#DeltaR(bb)'+unit,
#      texX = '#DeltaR(bb)', texY = labelAddon + '#frac{d#sigma}{d #DeltaR(bb)}'+unit,
      attribute = lambda event, sample: event.bbdR if event.passing_checks else float('nan'),
      binning=[20,0,3.5],
    ) )
    fisherInfoVariables.append(None)

    plots.append( Plot( name = 'Met_pt',
      texX = 'E_{T}^{miss} [GeV]', texY = y_label + ' / 10 GeV',
#      texX = 'E_{T}^{miss} [GeV]', texY = labelAddon + 'd#sigma/dE_{T}^{miss}'+unit,
#      texX = 'E_{T}^{miss} [GeV]', texY = labelAddon + '#frac{d#sigma}{d E_{T}^{miss}}'+unit,
      attribute = lambda event, sample: event.MET['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,200],
    ) )
    fisherInfoVariables.append('%sMet_pt'%preTag)

    plots.append(Plot( name = 'njets',
      texX = 'Number of Jets', texY = y_label,
      attribute = lambda event, sample: len(event.jets) if event.passing_checks else float('nan'),
      binning=[10,0,10],
    ) )
    fisherInfoVariables.append('n%sJet'%preTag)


    plots.append(Plot( name = 'nleps',
      texX = 'Number of Leptons', texY = y_label,
      attribute = lambda event, sample: getattr( event, 'n%sLep'%preTag ) if event.passing_checks else float('nan'),
      binning=[8,0,8],
    ) )
    fisherInfoVariables.append('n%sLep'%preTag)


    plots.append(Plot( name = 'nbjets',
      texX = 'Number of bJets', texY = y_label,
      attribute = lambda event, sample: getattr( event, 'nBTag' ) if event.passing_checks else float('nan'),
      binning=[4,0,4],
    ))
    fisherInfoVariables.append(None)

    if len(plots) != len(fisherInfoVariables):
        raise ValueError('Wrong plot list in ttZ_3l! plots and fisherInfoVariables lists must be same size!')

    return plots, fisherInfoVariables

