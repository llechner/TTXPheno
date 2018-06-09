''' Class for reading Delphes files.
    Based on Suchitas class
'''

# Standard imports
import ROOT

# Suchis reader
ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/TTXPheno/Tools/scripts/DelphesRecoClass.C")

class DelphesReader:


    def __init__( self, filename ):
        ''' Return a Suchita-reader
        '''
        self.reader =  ROOT.DelphesRecoClass(filename)

    # Read a vector collection from the Delphes reader
    def read_collection( self, collection, variables ):
        nColl   = getattr( self.reader, collection+"_size" )
        buffers = {var_old: getattr( self.reader, collection+'_'+var_old) for var_old, var_new in variables}
        return [{var_new:buffers[var_old][i] for var_old, var_new in variables} for i in range(nColl)]

    def muons( self ):
        res = self.read_collection( 'Muon', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'), 
                ('Charge', 'charge'), ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),  
                ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt') 
            ])
        for r in res:
            r['pdgId'] = -13*r['charge']
        return res

    def electrons( self ):
        res = self.read_collection( 'Electron', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'), 
                ('Charge', 'charge'), ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),  
                ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt'),
                ('EhadOverEem','ehadOverEem')
            ])
        for r in res:
            r['pdgId'] = -11*r['charge']
        return res

    def jets( self ):
        return self.read_collection( 'Jet', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('BTag', 'bTag'), ( 'bTagPhys', 'bTagPhys')
            ])

    def photons( self ):
        return self.read_collection( 'Jet', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
            ])

    def met( self ):
        return self.read_collection( 'MissingET', [('MET', 'pt'), ('Phi', 'phi')] )

    def getEntry( self, entry ):
        # Set the Suchi-Reader 
        self.reader.GetEntry( entry )
