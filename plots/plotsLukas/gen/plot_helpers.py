# Standard imports and batch mode
import ROOT
ROOT.gROOT.SetBatch(True)
from math                                import cos, sin, sinh, cosh

# Helpers
def addTransverseVector( p_dict ):
    ''' add a transverse vector for further calculations
    '''
    p_dict['vec2D'] = ROOT.TVector2( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']) )

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vec4D'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), 0 )

def NanJet():
    ''' return a dict in Jet format filled with Nan
    '''
    return {'index':float('nan'), 'pt':float('nan'), 'phi':float('nan'), 'eta':float('nan'), 'matchBParton':float('nan'), 'vec2D':ROOT.TVector2( float('nan'), float('nan') ), 'vec4D':ROOT.TLorentzVector( float('nan'), float('nan'), float('nan'), float('nan') )}

def NanLepton():
    ''' return a dict in Lepton format filled with Nan
    '''
    return {'index':float('nan'), 'pt':float('nan'), 'phi':float('nan'), 'pdgId':float('nan'), 'eta':float('nan'), 'motherPdgId':float('nan'), 'vec2D':ROOT.TVector2( float('nan'), float('nan') ), 'vec4D':ROOT.TLorentzVector( float('nan'), float('nan'), float('nan'), float('nan') )}

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
