''' Class to interpret weight info pkl file
'''

# General
import pickle
import scipy.special
import itertools

from operator import mul

# Logger
import logging
logger = logging.getLogger(__name__)

class WeightInfo:
    def __init__( self, filename ):
        self.data      = pickle.load(file(filename))

        self.variables = self.data.keys()[0].split('_')[::2]
        self.nvar      = len(self.variables)

        # Sort wrt to position in ntuple
        self.id = self.data.keys()
        self.id.sort(key=lambda w: self.data[w])
        self.nid = len(self.id)

        logger.debug( "Found %i variables: %s. Found %i weights." %(self.nvar, ",".join( self.variables ), self.nid) )

    def set_order( self, order):
#        gp_order = get_pkl_order(...) #get order of gridpack from pkl file
#        if order > gp_order:
#            raise ValueError( "Polynomial order is greater than in the gridpack (order %i)" % gp_order )
        self.order = order

    @staticmethod
    def get_ndof( nvar, order ):
        return sum( [ int(scipy.special.binom(nvar + o - 1, o)) for o in xrange(order+1) ] )

    def weight_string(self):
        substrings = []
        counter = 0
        for o in xrange(self.order+1):
            for comb in itertools.combinations_with_replacement( self.variables, o ):
                substrings.append(  "*".join( ["p_C[%i]"%counter] + [ "rw_%s"%v for v in  comb] )  )
                counter += 1

        return "+".join( substrings )

    def arg_weight_string(self, **kwargs):
        if len(kwargs)==0: return 'p_C[0]'
        unused_args = set(kwargs.keys()) - set(self.variables)
        if len(unused_args) > 0:
            raise ValueError( "Variable %s not in the gridpack" % ' && '.join(unused_args) )
        substrings = []
        counter = -1
        for o in xrange(self.order+1):
            for comb in itertools.combinations_with_replacement( self.variables, o ):
                counter += 1
                if False in [v in kwargs for v in comb]: continue
                substrings.append( "p_C[%i]*%s" %(counter, str(reduce(mul,[kwargs.get(v) for v in comb],1)).rstrip('0').rstrip('.')) )
        return "+".join( substrings )

    @staticmethod
    def differentiate( comb, var ):
        ''' Differentiate a product
        '''
        prefac = comb.count( var ) 
        if prefac==0:
            diff_comb = tuple()
        else:
            diff_comb = list(comb)
            diff_comb.remove( var )
        return prefac, tuple(diff_comb)

    def diff_weight_string(self, var):
        if var not in self.variables:
            raise ValueError( "Variable %s not in list of variables %r" % (var, self.variables) )
        substrings = []
        counter = 0
        for o in xrange(self.order+1):
            for comb in itertools.combinations_with_replacement( self.variables, o ):
                prefac, diff_comb = WeightInfo.differentiate( comb, var)
                if prefac!=0:
                    substrings.append(  "*".join( ["%i*p_C[%i]"%(prefac, counter) if prefac!=1 else "p_C[%i]"% counter] + [ "rw_%s"%v for v in diff_comb] )  )
                counter += 1

        return "+".join( substrings )

    def FisherParametrization( self, var1, var2):
        if var1==var2:
            return "(%s)**2/(%s)"%( self.diff_weight_string( var1 ), self.weight_string() )
        else:
            return "(%s)*(%s)/(%s)"%( self.diff_weight_string( var1 ), self.diff_weight_string( var2 ), self.weight_string( ) )


if __name__ == "__main__":

    #w = WeightInfo("/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_currentplane_highStat_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl")
    #w.FisherParametrization(2, 'cpt', 'cpt')
    #c = ROOT.TChain( "Events" )
    #c.Add("/afs/hephy.at/data/rschoefbeck02/TopEFT/skims/gen/v2/fwlite_ttZ_ll_LO_currentplane_highStat_scan/fwlite_ttZ_ll_LO_currentplane_highStat_scan_0.root" )
    #fisher_string = ":".join( [ w.FisherParametrization(2, 'cpt', 'cpt'),  w.FisherParametrization(2, 'cpt', 'cpQM'),  w.FisherParametrization(2, 'cpQM', 'cpQM') ] )
    #c.Scan(fisher_string)
    import ROOT
    c = ROOT.TChain("Events")
    c.Add("/afs/hephy.at/data/rschoefbeck02/TopEFT/skims/gen/v2_small/fwlite_ttZ_ll_LO_highStat_scan/fwlite_ttZ_ll_LO_highStat_scan.root")
    w = WeightInfo("/afs/hephy.at/data/rschoefbeck02/TopEFT/results/gridpacks/ttZ0j_rwgt_patch_625_slc6_amd64_gcc630_CMSSW_9_3_0_tarball.pkl")
    w.set_order( 3 )
#    fisher_string = ":".join( [ w.FisherParametrization( 'cpt', 'cpt'),  w.FisherParametrization( 'cpt', 'cpQM'),  w.FisherParametrization('cpQM', 'cpQM') ] )

    print(w.arg_weight_string(ctZI=2, cpt=5, ctZ=.4, cmk=3))
#     print(w.arg_weight_string())
