''' Upgrade JEC
'''
import os
from TTXPheno.Tools.helpers import getObjFromFile

class UpgradeJECUncertainty:

    def __init__( self ):
       self.light_jet_JEC = getObjFromFile( os.path.expandvars( "$CMSSW_BASE/src/TTXPheno/Tools/data/HL_YR_JEC.root"), "TOTAL_DIJET_AntiKt4EMTopo_YR2018" )
       self.b_jet_JEC     = getObjFromFile( os.path.expandvars( "$CMSSW_BASE/src/TTXPheno/Tools/data/HL_YR_JEC.root"), "TOTAL_BJES_AntiKt4EMTopo_YR2018" )

    def applyJECInfo(self, jet, flavor ):
        if abs(flavor) == 5:
            h = self.b_jet_JEC
        else:
            h = self.light_jet_JEC
        fac = h.GetBinContent( h.FindBin( jet['pt'] ) )
        jet['pt_JEC_up'] = jet['pt']*(1+fac)
        jet['pt_JEC_down'] = jet['pt']/(1+fac) 
