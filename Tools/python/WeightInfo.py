''' Class to interpret weight info pkl file
'''

# General
import pickle
import scipy.special
import itertools

from operator import mul

import numpy as np

# Logger
import logging
logger = logging.getLogger(__name__)

class WeightInfo:
    def __init__( self, filename ):
        data = pickle.load(file(filename))

        if 'rw_dict' in data.keys(): self.data = data['rw_dict']
        else: self.data = data

        if 'order' in data.keys(): self.pkl_order = data['order']['order']
        else: self.pkl_order = None

        if 'ref_point' in data.keys(): self.ref_point = data['ref_point']
        else: self.ref_point = None

        # store all variables (Wilson coefficients)
        self.variables = self.data.keys()[0].split('_')[::2]
        self.nvar      = len(self.variables)

        # compute reference point coordinates
        self.ref_point_coordinates = { var: float( self.ref_point[var] ) if (self.ref_point is not None and var in self.ref_point.keys()) else 0 for var in self.variables }

        # Sort wrt to position in ntuple
        self.id = self.data.keys()
        self.id.sort(key=lambda w: self.data[w])
        self.nid = len(self.id)

        logger.debug( "Found %i variables: %s. Found %i weights." %(self.nvar, ",".join( self.variables ), self.nid) )

    def set_order( self, order):
        if self.pkl_order == None:
            print( "WARNING: Could not find the polynomial order of the gridpack!")
        elif order > self.pkl_order:
            raise ValueError( "Polynomial order is greater than in the gridpack (order %i)" % self.pkl_order )
        self.order = order

    @staticmethod
    def get_ndof( nvar, order ):
        return sum( [ int(scipy.special.binom(nvar + o - 1, o)) for o in xrange(order+1) ] )

    # compute combinations on demand
    @property
    def combinations( self ):
        if hasattr( self, "_combinations"):
            return self._combinations
        else:
            self._combinations = []
            for o in xrange(self.order+1):
                self._combinations.extend( list(itertools.combinations_with_replacement( self.variables, o )) )
            return self._combinations

    def weight_string(self):
        ''' get the full reweight string
        '''
        substrings = []
        for i_comb, comb in enumerate(self.combinations):
            substrings.append(  "*".join( ["p_C[%i]"%i_comb] + [ "(rw_%s-%s)" %( v, str(float(self.ref_point[v])).rstrip('0') ) if self.ref_point is not None and v in self.ref_point.keys() else "rw_%s"%(v) for v in  comb] )  )

        return "+".join( substrings )

    def complement_args(self, args ):
        ''' prepare the args; add the ref_point ones and check that there is no inconsistency
        '''

        # Remove zeros
        for x,y in args.iteritems():
            if y==0.: del args[x] 

        # add WC that are in the ref point but not in args
        if self.ref_point is not None:
            for item in self.ref_point.keys():
                if item not in args.keys(): args[item] = 0

        # check if WC in args that are not in the gridpack
        unused_args = set(args.keys()) - set(self.variables)
        if len(unused_args) > 0:
            raise ValueError( "Variable %s not in the gridpack! Please use only the following variables: %s" % (' && '.join(unused_args), ', '.join(self.variables)) )

    def get_weight_string(self, **kwargs):
        '''make a root draw string that evaluates the weight in terms of the p_C coefficient vector using the kwargs as WC
        '''
        # add the arguments from the ref-point 
        self.complement_args( kwargs )

        substrings = []
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            substrings.append( "p_C[%i]*%s" %(i_comb, str(float(reduce(mul, [ (float(kwargs[v]) - self.ref_point_coordinates[v]) for v in comb],1))).rstrip('0') ) )
        return "+".join( substrings )

    def get_weight_func(self, **kwargs):
        '''construct a lambda function that evaluates the weight in terms of the event.p_C coefficient vector using the kwargs as WC
        '''
        # add the arguments from the ref-point 
        self.complement_args( kwargs )

        terms = []
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            # store [ ncoeff, factor ]
            terms.append( [ i_comb, float(reduce(mul,[ ( float(kwargs[v]) - self.ref_point_coordinates[v]) for v in comb],1)) ] )

        return lambda event, sample: sum( event.p_C[term[0]]*term[1] for term in terms )

    def get_weight_yield(self, coeffLists, **kwargs):
        '''compute yield from a list of coefficients (in the usual order of p_C) using the kwargs as WC
        '''

        # add the arguments from the ref-point 
        self.complement_args( kwargs )

        result = 0 
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            for coeffList in coeffLists:
                if coeffList[i_comb]==0: continue
                result += coeffList[i_comb]*float(reduce(mul,[ ( float(kwargs[v]) - self.ref_point_coordinates[v]) for v in comb],1))

        return result

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
        ''' return string of the differentiated full weight string
        '''

        if var not in self.variables:
            raise ValueError( "Variable %s not in list of variables %r" % (var, self.variables) )

        substrings = []
        for i_comb, comb in enumerate(self.combinations):
            prefac, diff_comb = WeightInfo.differentiate( comb, var)
            if prefac!=0:
                substrings.append(  "*".join( ["%i*p_C[%i]"%(prefac, i_comb) if prefac!=1 else "p_C[%i]"% i_comb] + [ "(rw_%s-%s)"%(v, str(float(self.ref_point[v])).rstrip('0')) if self.ref_point is not None and v in self.ref_point.keys() else "rw_%s"%(v) for v in  diff_comb] ) )

        return "+".join( substrings )

    def get_diff_weight_string(self, var, **kwargs):
        '''make a root draw string that evaluates the diff weight 
           in terms of the p_C coefficient vector using the kwargs as WC
        '''

        # add the arguments from the ref-point 
        self.complement_args( kwargs )

        substrings = []
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            prefac, diff_comb = WeightInfo.differentiate( comb, var)
            if prefac == 0: continue
            substrings.append( "p_C[%i]*%s" %(i_comb, str(float(reduce(mul,[ ( float(kwargs[v]) - float(self.ref_point[v]) ) if self.ref_point is not None and v in self.ref_point.keys() else float(kwargs[v]) for v in diff_comb] + [prefac],1))).rstrip('0') ) )

        return "+".join( substrings )


    def get_diff_weight_func(self, var, **kwargs):
        '''construct a lambda function that evaluates the diff weight in terms of the event.p_C coefficient vector using the kwargs as WC
        '''

        # add the arguments from the ref-point 
        self.complement_args( kwargs )

        terms = []
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            prefac, diff_comb = WeightInfo.differentiate( comb, var)
            if prefac == 0: continue
            # store [ ncoeff, factor ]
            terms.append( [ i_comb, float(reduce(mul,[ ( float(kwargs[v]) - float(self.ref_point[v]) ) if self.ref_point is not None and v in self.ref_point.keys() else float(kwargs[v]) for v in diff_comb] + [prefac],1)) ] )

        return lambda event, sample: sum( event.p_C[term[0]]*term[1] for term in terms )


    def get_diff_weight_yield(self, var, coeffLists, **kwargs):
        '''compute diff yield from a list of coefficients (in the usual order of p_C) using the kwargs as WC
        '''

        # add the arguments from the ref-point 
        self.complement_args( kwargs )

        result = 0 
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            prefac, diff_comb = WeightInfo.differentiate( comb, var)
            if prefac == 0: continue
            for coeffList in coeffLists:
                if coeffList[i_comb]==0: continue
                result += coeffList[i_comb]*float(reduce(mul,[ ( float(kwargs[v]) - float(self.ref_point[v]) ) if self.ref_point is not None and v in self.ref_point.keys() else float(kwargs[v]) for v in diff_comb] + [prefac],1))

        return result

    def fisherParametrization( self, var1, var2):
        ''' return string of the fisher information matrix entry ij
        '''

        if var1==var2:
            return "(%s)**2/(%s)"%( self.diff_weight_string( var1 ), self.weight_string() )
        else:
            return "(%s)*(%s)/(%s)"%( self.diff_weight_string( var1 ), self.diff_weight_string( var2 ), self.weight_string( ) )

    def get_fisherParametrization_entry( self, var1, var2, coeffLists, **kwargs):
        ''' return the value of the fisher information matrix entry ij
        '''

        if var1==var2:
            return self.get_diff_weight_yield( var1, coeffLists, **kwargs )**2 / self.get_weight_yield( coeffLists, **kwargs )
        else:
            return self.get_diff_weight_yield( var1, coeffLists, **kwargs ) * self.get_diff_weight_yield( var2, coeffLists, **kwargs ) / self.get_weight_yield( coeffLists, **kwargs )

    def get_fisherInformation_matrix( self, coeffLists, **kwargs):
        ''' return the full fisher information matrix
        '''

        # initialize FI matrix
        fi_matrix = np.full((self.nvar, self.nvar), float('nan'))

        for i in range( self.nvar ):
            for j in range( self.nvar ):
                if i>j: fi_matrix[i,j] = fi_matrix[j,i]
                fi_matrix[i,j] = self.get_fisherParametrization_entry( self.variables[i], self.variables[j], coeffLists, **kwargs )

        return fi_matrix


    def matrix_to_string( self, matrix ):
        ''' return the matrix in a terminal visualization string (print)
        '''
        return '\n'.join( ['\t'.join( (map('{:.5f}'.format, item) + [self.variables[i-1]] if i!=0 else item) ) for i, item in enumerate([self.variables] + matrix.tolist()) ] )                

# Make a list from the bin contents from a histogram that resulted from a 'Draw' of p_C 
def histo_to_list( histo ):
    return [histo.GetBinContent(i) for i in range(1,histo.GetNbinsX()+1)]

if __name__ == "__main__":

    # RootTools
    from RootTools.core.standard import *
    # Logger
    import TTXPheno.Tools.logger as logger
    import RootTools.core.logger as logger_rt
    logger    = logger.get_logger(   'INFO', logFile = None)
    logger_rt = logger_rt.get_logger('INFO', logFile = None)
    # TTXPheno
    from TTXPheno.Tools.cutInterpreter import cutInterpreter
    from TTXPheno.samples.benchmarks import * 

    # Sample
    sample = fwlite_ttZ_ll_LO_order2_15weights_ref
    sample.reduceFiles( to = -1 )
    w = WeightInfo(sample.reweight_pkl)
    w.set_order( 2 )

    selection_string = cutInterpreter.cutString('lepSel3-onZ-njet3p-nbjet1p')

    # Make a coeff histo from a sample
    def getCoeffListFromDraw( sample, selectionString, weightString = None):
        histo = sample.get1DHistoFromDraw( 
            "Iteration$", 
            [ len(w.combinations), 0, len(w.combinations) ], 
            selectionString = selectionString, 
            weightString = 'p_C*(%s)'%weightString if weightString is not None else 'p_C' )
        return histo_to_list( histo )

    # Make a coeff histo from a sample
    def getCoeffPlotFromDraw( sample, variableString, binning, selectionString, weightString = None):
        # 2D Plot, Iteration$ is on x
        histo = sample.get2DHistoFromDraw( 
            "Iteration$:%s"%variableString, 
            [ len(w.combinations), 0, len(w.combinations) ] + binning, 
            selectionString = selectionString, 
            weightString = 'p_C*(%s)'%weightString if weightString is not None else 'p_C' )

        return [ histo_to_list(histo.ProjectionX("%i_px"%i, i+1, i+1)) for i in range( histo.GetNbinsY() ) ]
#        return histo_to_list( histo )

    # Fisher information in ptZ histo
    coeff_Z_pt = getCoeffPlotFromDraw( sample, 'Z_pt', [ 20, 0, 500 ], selection_string, weightString='150*lumiweight1fb')
    # Fisher information in x-sec
    coeff_tot = getCoeffListFromDraw( sample, selection_string, weightString='150*lumiweight1fb')

#
    print w.matrix_to_string(w.get_fisherInformation_matrix([coeff_tot]))
    print w.matrix_to_string(w.get_fisherInformation_matrix(coeff_Z_pt))
