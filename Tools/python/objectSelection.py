''' ObjectSelections
'''

def isGoodGenJet( j ):
    ''' jet object selection
    '''
    return j['pt'] > 30 and abs( j['eta'] ) < 2.4

def isGoodGenLepton( l ):
    ''' lepton object selection
    '''
    return l['pt'] > 10 and abs( l['eta'] ) < 2.5 and abs( int(l['pdgId']) ) in [11,13]
