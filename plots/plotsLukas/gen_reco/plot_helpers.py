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

