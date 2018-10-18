''' ObjectSelections
'''


def genFraction( fwlite_jet, pdgId ):
    return sum( [ fwlite_jet.getGenConstituent(i).pt() for i in range(fwlite_jet.numberOfSourceCandidatePtrs()) if abs(fwlite_jet.getGenConstituent(i).pdgId()) == pdgId ], 0 )/fwlite_jet.pt()

def genJetId( fwlite_jet ):
    return genFraction(fwlite_jet, 13)<0.8 and genFraction(fwlite_jet, 11)<0.8

def isGoodGenJet( j ):
    ''' jet object selection
    '''
    return j['pt'] > 30 and abs( j['eta'] ) < 4.0 #2.4

def isGoodGenPhoton( j ):
    ''' jet object selection
    '''
    return j['pt'] > 15 and abs( j['eta'] ) < 2.1

def isGoodGenLepton( l ):
    ''' lepton object selection
    '''
    return l['pt'] > 10 and abs( l['eta'] ) < 3.0 and abs( int(l['pdgId']) ) in [11,13] #eta < 2.5

def isGoodRecoLepton( l ):
    return l['pt'] > 10 and abs( l['eta'] ) < 3.0 and abs( int(l['pdgId']) ) in [11,13] #eta < 2.5

def isGoodRecoMuon( l ):
    return abs( l['pdgId'] ) == 13 and abs( l['eta'] ) < 3.0 and l['pt'] > 10 #eta < 2.5

def isGoodRecoElectron( l ):
    return abs( l['pdgId'] ) == 11 and abs( l['eta'] ) < 3.0 and l['pt'] > 10 #eta < 2.5

def isGoodRecoJet( j, pt_var = 'pt'):
    return  abs( j['eta'] ) < 4.0 and j[pt_var] > 30 and j['nCharged']>1 and j['nNeutrals']>0 #eta < 2.4

def isGoodRecoPhoton( g ):
    return  abs( g['eta'] ) < 2.1 and g['pt'] > 15

