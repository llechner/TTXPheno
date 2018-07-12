# Standard imports and batch mode
import ROOT
ROOT.gROOT.SetBatch(True)
from math                                import cos, sin, sinh, cosh

# RootTools
from RootTools.core.standard             import *

# Helpers
def addTransverseVector( p_dict ):
    ''' add a transverse vector for further calculations
    '''
    p_dict['vec2D'] = ROOT.TVector2( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']) )

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vec4D'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), 0 )

def UnitVectorT2( phi ):
    ''' 2D Unit Vector
    '''
    return ROOT.TVector2( cos(phi), sin(phi) )

def MTSquared( p1, p2 ):
    ''' compute MT from 2 particles
    '''
    return 2*p1['pt']*p2['pt']*( 1-cos(p1['phi']-p2['phi']) )

def MSquared( p1, p2 ):
    ''' compute MassSquared from 2 particles
    '''
    return 2*p1['pt']*p2['pt']*( cosh(p1['eta']-p2['eta'])-cos(p1['phi']-p2['phi']) )

def addIDeltaBeta( p_dict ):
    ''' add I_rel^{DeltaBeta}
    '''
    p_dict['IDeltaBeta'] = ( p_dict['sumPtCharged'] + p_dict['sumPtNeutral'] - 0.5* p_dict['sumPtChargedPU'] ) / ( p_dict['pt'] )

def getLeptonIsolationPlotList( particleList, y_label ):
    ''' return plotlist for Isolation plots
    '''
    plots = []

    for p in particleList:

        plots.append(Plot( name = "%s_isolationVar_%s"%( p['particleString'], p['eventString'] ),
          texX = 'isolationVar(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample: getattr( event, p['eventString'] )['isolationVar'] if abs( getattr( event, p['eventString'] )['pdgId'] ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.25],
        ))

        plots.append(Plot( name = "%s_isolationVarRhoCorr_%s"%( p['particleString'], p['eventString'] ),
          texX = 'isolationVarRhoCorr(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['isolationVarRhoCorr'] if abs( getattr( event, p['eventString'] )['pdgId'] ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.25],
        ))

        plots.append(Plot( name = "%s_sumPtCharged_%s"%( p['particleString'], p['eventString'] ),
          texX = 'sumPtCharged(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['sumPtCharged'] if abs( getattr( event, p['eventString'] )['pdgId'] ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,15],
        ))

        plots.append(Plot( name = "%s_sumPtNeutral_%s"%( p['particleString'], p['eventString'] ),
          texX = 'sumPtNeutral(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['sumPtNeutral'] if abs( getattr( event, p['eventString'] )['pdgId'] ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,15],
        ))

        plots.append(Plot( name = "%s_sumPtChargedPU_%s"%( p['particleString'], p['eventString'] ),
          texX = 'sumPtChargedPU(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['sumPtChargedPU'] if abs( getattr( event, p['eventString'] )['pdgId'] ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))

        plots.append(Plot( name = "%s_sumPt_%s"%( p['particleString'], p['eventString'] ),
          texX = 'sumPt(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['sumPt'] if abs( getattr( event, p['eventString'] )['pdgId'] ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[15,0,25],
        ))

        plots.append(Plot( name = "%s_ehadOverEem_%s"%( p['particleString'], p['eventString'] ),
          texX = 'ehadOverEem(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['ehadOverEem'] if abs( getattr( event, p['eventString'] )['pdgId'] ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))

        plots.append(Plot( name = "%s_IDeltaBeta_%s"%( p['particleString'], p['eventString'] ),
          texX = 'I^{#Delta#beta}_{rel}(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['IDeltaBeta'] if abs( getattr( event, p['eventString'] )['pdgId'] ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.3],
        ))

    return plots


def getPhotonIsolationPlotList( particleList, y_label ):
    ''' return plotlist for Isolation plots for photons
    '''
    plots = []

    for p in particleList:

        plots.append(Plot( name = "%s_isolationVar_%s"%( p['particleString'], p['eventString'] ),
          texX = 'isolationVar(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample: getattr( event, p['eventString'] )['isolationVar'] if event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))

        plots.append(Plot( name = "%s_isolationVarRhoCorr_%s"%( p['particleString'], p['eventString'] ),
          texX = 'isolationVarRhoCorr(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['isolationVarRhoCorr'] if event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))

        plots.append(Plot( name = "%s_sumPtCharged_%s"%( p['particleString'], p['eventString'] ),
          texX = 'sumPtCharged(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['sumPtCharged'] if event.passing_checks else float('nan'),
          binning=[20,0,10],
        ))

        plots.append(Plot( name = "%s_sumPtNeutral_%s"%( p['particleString'], p['eventString'] ),
          texX = 'sumPtNeutral(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['sumPtNeutral'] if event.passing_checks else float('nan'),
          binning=[20,0,10],
        ))

        plots.append(Plot( name = "%s_sumPtChargedPU_%s"%( p['particleString'], p['eventString'] ),
          texX = 'sumPtChargedPU(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['sumPtChargedPU'] if event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))

        plots.append(Plot( name = "%s_sumPt_%s"%( p['particleString'], p['eventString'] ),
          texX = 'sumPt(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['sumPt'] if event.passing_checks else float('nan'),
          binning=[15,0,20],
        ))

        plots.append(Plot( name = "%s_ehadOverEem_%s"%( p['particleString'], p['eventString'] ),
          texX = 'ehadOverEem(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['ehadOverEem'] if event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))

        plots.append(Plot( name = "%s_genIndex_%s"%( p['particleString'], p['eventString'] ),
          texX = 'genIndex(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['genIndex'] if event.passing_checks else float('nan'),
          binning=[6,-2,4],
        ))

        plots.append(Plot( name = "%s_IDeltaBeta_%s"%( p['particleString'], p['eventString'] ),
          texX = 'I^{#Delta#beta}_{rel}(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] )['IDeltaBeta'] if event.passing_checks else float('nan'),
          binning=[20,0,0.15],
        ))

    return plots

