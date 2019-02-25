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

def getLeptonIsolationPlotList( particleList, y_label, zoom=False ):
    ''' return plotlist for Isolation plots
    '''
    plots = []

    for p in particleList:

        plots.append(Plot( name = "%s_isolationVar_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'isolationVar(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample: getattr( event, p['eventString'] + '_isolationVar[%i]'%p['index'] ) if abs( getattr( event, p['eventString'] + '_pdgId[%i]'%p['index'] ) ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.1] if zoom else [20,0,0.25],
        ))

        plots.append(Plot( name = "%s_isolationVarRhoCorr_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'isolationVarRhoCorr(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_isolationVarRhoCorr[%i]'%p['index'] ) if abs( getattr( event, p['eventString'] + '_pdgId[%i]'%p['index'] ) ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.1] if zoom else [20,0,0.25],
        ))

        plots.append(Plot( name = "%s_sumPtCharged_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'sumPtCharged(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_sumPtCharged[%i]'%p['index'] ) if abs( getattr( event, p['eventString'] + '_pdgId[%i]'%p['index'] ) ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,2] if zoom else [20,0,15],
        ))

        plots.append(Plot( name = "%s_sumPtNeutral_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'sumPtNeutral(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_sumPtNeutral[%i]'%p['index'] ) if abs( getattr( event, p['eventString'] + '_pdgId[%i]'%p['index'] ) ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,2] if zoom else [20,0,15],
        ))

        plots.append(Plot( name = "%s_sumPtChargedPU_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'sumPtChargedPU(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_sumPtChargedPU[%i]'%p['index'] ) if abs( getattr( event, p['eventString'] + '_pdgId[%i]'%p['index'] ) ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,2] if zoom else [20,0,0.1],
        ))

        plots.append(Plot( name = "%s_sumPt_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'sumPt(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_sumPt[%i]'%p['index'] ) if abs( getattr( event, p['eventString'] + '_pdgId[%i]'%p['index'] ) ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,5] if zoom else [15,0,25],
        ))

        plots.append(Plot( name = "%s_ehadOverEem_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'ehadOverEem(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_ehadOverEem[%i]'%p['index'] ) if abs( getattr( event, p['eventString'] + '_pdgId[%i]'%p['index'] ) ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))

        plots.append(Plot( name = "%s_IDeltaBeta_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'I^{#Delta#beta}_{rel}(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_IDeltaBeta[%i]'%p['index'] ) if abs( getattr( event, p['eventString'] + '_pdgId[%i]'%p['index'] ) ) == p['pdg'] and event.passing_checks else float('nan'),
          binning=[20,0,0.05] if zoom else [20,0,0.3],
        ))

    return plots


def getPhotonIsolationPlotList( particleList, y_label, zoom=False ):
    ''' return plotlist for Isolation plots for photons
    '''
    plots = []

    for p in particleList:

        plots.append(Plot( name = "%s_isolationVar_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'isolationVar(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample: getattr( event, p['eventString'] + '_isolationVar[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[20,0,0.05] if zoom else [20,0,0.15],
        ))

        plots.append(Plot( name = "%s_isolationVarRhoCorr_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'isolationVarRhoCorr(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_isolationVarRhoCorr[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[20,0,0.05] if zoom else [20,0,0.15],
        ))

        plots.append(Plot( name = "%s_sumPtCharged_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'sumPtCharged(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_sumPtCharged[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[20,0,2] if zoom else [20,0,10],
        ))

        plots.append(Plot( name = "%s_sumPtNeutral_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'sumPtNeutral(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_sumPtNeutral[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[20,0,2] if zoom else [20,0,10],
        ))

        plots.append(Plot( name = "%s_sumPtChargedPU_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'sumPtChargedPU(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_sumPtChargedPU[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))

        plots.append(Plot( name = "%s_sumPt_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'sumPt(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_sumPt[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[20,0,5] if zoom else [15,0,20],
        ))

        plots.append(Plot( name = "%s_ehadOverEem_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'ehadOverEem(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_ehadOverEem[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[20,0,0.1],
        ))

        plots.append(Plot( name = "%s_genIndex_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'genIndex(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_genIndex[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[6,-2,4],
        ))

        plots.append(Plot( name = "%s_IDeltaBeta_%s%s"%( p['particleString'], p['eventString'], '_zoom' if zoom else ''),
          texX = 'I^{#Delta#beta}_{rel}(%s)'%p['eventString'], texY = y_label,
          attribute = lambda event, sample:  getattr( event, p['eventString'] + '_IDeltaBeta[%i]'%p['index'] ) if event.passing_checks else float('nan'),
          binning=[20,0,0.05] if zoom else [20,0,0.15],
        ))

    return plots

