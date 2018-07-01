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
    event.bj0lep = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjLeadlep_index'%preTag ) ) if getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 else NanJet()
    event.bj0had = getObjDict( event, '%sJet_'%preTag, ['pt', 'eta', 'phi'], getattr( event, '%sBjLeadhad_index'%preTag ) ) if getattr( event, '%sBjLeadlep_index'%preTag ) >= 0 else NanJet()

    # Add extra vectors
    for p in [event.bj0, event.bj1, event.bj0lep, event.bj0had]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco':      from TTXPheno.Tools.objectSelection      import isGoodRecoJet     as isGoodJet
    else:                    from TTXPheno.Tools.objectSelection      import isGoodGenJet      as isGoodJet

    # selection checks
    event.foundBj0 = isGoodJet( event.bj0 )
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
#    event.l1 = getObjDict( event, '%sLep_'%preTag, ['pt', 'eta', 'phi', 'pdgId'], 1 )

    # Add extra vectors
    for p in [ event.l0 ]:
        addTransverseVector( p )
        addTLorentzVector( p )

    # Import additional functions/classes specified for the level of reconstruction
    if level == 'reco': from TTXPheno.Tools.objectSelection      import isGoodRecoLepton  as isGoodLepton
    else:               from TTXPheno.Tools.objectSelection      import isGoodGenLepton   as isGoodLepton

    # We may loose some events by cross-cleaning or by thresholds.
    event.found1lep    = isGoodLepton( event.l0 )

    # choose your selection on leptons
    event.passing_leptons = event.found1lep


def makeObservables( event, sample, level):
    ''' Compute all relevant observables
    '''
    # double b kinematic
    event.bbdPhi = deltaPhi( event.bj0['phi'], event.bj1['phi'] )
    event.bbdR   = deltaR( event.bj0, event.bj1 )

    # Make leptonic W, Nan if any component is Nan
    event.Wlep_vec2D = event.MET['vec2D'] + event.l0['vec2D']
    # Lp, Nan if any component is Nan
    event.Wlep_Lp = ( event.Wlep_vec2D*event.l0['vec2D'] ) / ( event.Wlep_vec2D*event.Wlep_vec2D )
    # classic MT, Nan if any component is Nan
    event.Wlep_MT = sqrt( MTSquared( event.MET, event.l0 ) )

    # blep+lep subsystem, Nan if len(event.lepsNotFromZ == 0)
    event.bleplep_vec2D = event.bj0lep['vec2D'] + event.l0['vec2D']
    event.bleplep_vec4D = event.bj0lep['vec4D'] + event.l0['vec4D']

    # transverse mass of top, Nan if len(event.lepsNotFromZ == 0)
    event.t_MT = sqrt( MTSquared( event.l0, event.MET ) + MSquared( event.l0, event.bj0lep ) + MTSquared( event.MET, event.bj0lep ) )
    event.t_vec2D = event.l0['vec2D'] + event.MET['vec2D'] + event.bj0lep['vec2D']

    # signed lepton pt, Nan if len(event.lepsNotFromZ == 0)
    event.lep0chargept = event.l0['pt'] if event.l0['pdgId']>0 else -event.l0['pt']
#    event.lep1chargept = event.l1['pt'] if event.l1['pdgId']>0 else -event.l1['pt']

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
    
    plots.append(Plot( name = 'Met_pt',
      texX = 'E_{T}^{miss} [GeV]', texY = y_label,
      attribute = lambda event, sample: event.MET['pt'] if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name	= 'Met_phi',
      texX = '#phi(E_{T}^{miss})', texY = y_label,
      attribute = lambda event, sample: event.MET['phi'] if event.passing_checks else float('nan'),
      binning=[20,-pi,pi],
    ))
    
    plots.append(Plot( name = 'W_pt',
      texX = 'p_{T}(W_{lep}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.Wlep_vec2D.Mod() if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name = 'W_Lp',
      texX = 'L_{p} from W_{lep}', texY = y_label,
      attribute = lambda event, sample: event.Wlep_Lp if event.passing_checks else float('nan'),
      binning=[20,-3,3],
    ))
    
    plots.append(Plot( name = 'bleplep_dot_ngamma_2D',
      texX = 'p_{T}(b_{lep} + l) [GeV] #upoint n(#gamma) (2D)', texY = y_label,
      attribute = lambda event, sample: event.bleplep_vec2D*event.gamma_unitVec2D if event.passing_checks else float('nan'),
      binning=[20,-400,400],
    ))
    
    plots.append(Plot( name = 'bleplep_dot_ngamma_3D',
      texX = 'p_{T}(b_{lep} + l) [GeV] #upoint n(#gamma) (3D)', texY = y_label,
      attribute = lambda event, sample: event.bleplep_vec4D.Vect()*event.gamma_unitVec3D if event.passing_checks else float('nan'),
      binning=[20,-400,400],
    ))
    
    plots.append(Plot( name = 'top_dot_ngamma',
      texX = 'p_{T}(t_{lep}) [GeV] #upoint n(#gamma)', texY = y_label,
      attribute = lambda event, sample: event.t_vec2D*event.gamma_unitVec2D if event.passing_checks else float('nan'),
      binning=[20,-400,400],
    ))
    
    plots.append(Plot( name = 'top_lep_pt',
      texX = 'p_{T}(t_{lep}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.t_vec2D.Mod() if event.passing_checks else float('nan'),
      binning=[20,0,400],
    ))
    
    plots.append(Plot( name = 'l0_pt_charge',
      texX = 'p_{T}(l_{0}) [GeV] signed with lepton charge', texY = y_label,
      attribute = lambda event, sample: event.lep0chargept if event.passing_checks else float('nan'),
      binning=[20,-200,200],
    ))
    
    plots.append(Plot( name = 'mT_W',
      texX = 'm_{T}(W_{lep}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.Wlep_MT if event.passing_checks else float('nan'),
      binning=[20,0,150],
    ))
    
    plots.append(Plot( name = 'mT_t',
      texX = 'm_{T}(t_{lep}) [GeV]', texY = y_label,
      attribute = lambda event, sample: event.t_MT if event.passing_checks else float('nan'),
      binning=[20,0,300],
    ))
    
    return plots
