''' Class to interpret weight info pkl file
'''

# General
import pickle
import scipy.special
import scipy.linalg
import itertools

# TTXPheno
import TTXPheno.Tools.helpers as helpers

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
        else: 
            self.ref_point = None
            logger.warning( "No reference point found in pkl file!" )

        # store all variables (Wilson coefficients)
        self.variables = self.data.keys()[0].split('_')[::2]
        self.nvar      = len(self.variables)

        # compute reference point coordinates
        self.ref_point_coordinates = { var: float( self.ref_point[var] ) if ( self.ref_point is not None and var in self.ref_point.keys() ) else 0 for var in self.variables }

        # Sort wrt to position in ntuple
        self.id = self.data.keys()
        self.id.sort(key=lambda w: self.data[w])
        self.nid = len(self.id)

        logger.debug( "Found %i variables: %s. Found %i weights." %(self.nvar, ",".join( self.variables ), self.nid) )

    def set_order( self, order ):
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

    def weight_string_WC( self ):
        ''' get the full reweight string
        '''
        substrings = []
        for i_comb, comb in enumerate( self.combinations ):
            subsubstrings = [ "p_C[%i]" %i_comb ]
            for v in comb:
                if self.ref_point_coordinates[v] == 0:
                    subsubstrings.append( 'rw_%s' %v ) 
                else:
                    subsubstrings.append( "(rw_%s-%s)" %( v, str(self.ref_point_coordinates[v]).rstrip('0')) ) 

            substrings.append(  "*".join( subsubstrings )  )

        return "+".join( substrings )

    def set_default_args( self, args ):
        ''' prepare the args; add the ref_point ones and check that there is no inconsistency
        '''

        # append reference point
        for var in self.variables:
            if var not in args.keys():
                args[var]=0.

        # check if WC in args that are not in the gridpack
        unused_args = set(args.keys()) - set(self.variables)
        if len(unused_args) > 0:
            raise ValueError( "Variable %s not in the gridpack! Please use only the following variables: %s" % (' && '.join(unused_args), ', '.join(self.variables)) )

    def get_weight_string( self, **kwargs ):
        '''make a root draw string that evaluates the weight in terms of the p_C coefficient vector using the kwargs as WC
        '''
        # add the arguments from the ref-point 
        self.set_default_args( kwargs )

        substrings = []
        for i_comb, comb in enumerate( self.combinations ):
            if False in [v in kwargs for v in comb]: continue
            # remove 0 entries
            fac = float( reduce( mul, [ (float(kwargs[v]) - self.ref_point_coordinates[v]) for v in comb ], 1 ) )
            if abs(fac) == 0.: continue
            substrings.append( "p_C[%i]*%s" %(i_comb, str(fac).rstrip('0') ) )
        return "+".join( substrings )

    @staticmethod
    @helpers.memoized
    def differentiate( comb, var ):
        ''' Differentiate a polynomial wrt to a variable represented by a combination of terms.
            Returns prefactor new combination.
            d\dv_i (v_i^n * X) -> n v_i^(n-1) * X 
        '''

        if type(var)==type(""):

            prefac = comb.count( var ) 

            if prefac==0:
                diff_comb = tuple()
            else:
                diff_comb = list( comb )
                diff_comb.remove( var )

            return prefac, tuple( diff_comb )

        elif type(var)==type(()) or type(var)==type([]):
            if len(var)==0: 
                return 1, comb
            elif len(var)==1:
                return WeightInfo.differentiate( comb, var[0] )
            else:
                prefac0, comb_diff  = WeightInfo.differentiate( comb, var[0] )
                prefac1, comb_diff2 = WeightInfo.differentiate( comb_diff, var[1:] ) 
            return prefac0*prefac1, comb_diff2
            
    # String methods
    def diff_weight_string_WC(self, var):
        ''' return string of the full weight string, differentiated wrt to var as a function of all WC
        '''

        if var not in self.variables:
            raise ValueError( "Variable %s not in list of variables %r" % (var, self.variables) )

        substrings = []
        for i_comb, comb in enumerate( self.combinations ):
            prefac, diff_comb = WeightInfo.differentiate( comb, var )
            if prefac != 0:
                subsubstrings = [ "%i*p_C[%i]" %(prefac, i_comb) if prefac != 1 else "p_C[%i]" %i_comb ]
                for v in diff_comb:
                    if self.ref_point_coordinates[v] == 0:
                        subsubstrings.append( 'rw_%s'%v )
                    else:
                        subsubstrings.append(  "(rw_%s-%s)"%(v, str(float(self.ref_point[v])).rstrip('0')) )
                substrings.append( "*".join( subsubstrings ) ) 
        
        return "+".join( substrings )


    def fisher_parametrization_string_WC( self, var1, var2 ):
        ''' return a string for the fisher information vor variables var1, vars as a function of the weight coefficients and all WC 
        '''

        if var1 == var2:
            return "(%s)**2/(%s)"%( self.diff_weight_string_WC( var1 ), self.weight_string_WC() )
        else:
            return "(%s)*(%s)/(%s)"%( self.diff_weight_string_WC( var1 ), self.diff_weight_string_WC( var2 ), self.weight_string_WC() )

    def diff_weight_string( self, var, **kwargs ):
        '''make a root draw string that evaluates the diff weight 
           in terms of the p_C coefficient vector using the kwargs as WC
        '''

        if var not in self.variables:
            raise ValueError( "Variable %s not in list of variables %r" % (var, self.variables) )

        # add the arguments from the ref-point 
        self.set_default_args( kwargs )

        substrings = []
        for i_comb, comb in enumerate( self.combinations ):
            if False in [v in kwargs for v in comb]: continue
            prefac, diff_comb = WeightInfo.differentiate( comb, var )
            if prefac == 0: continue
            fac = prefac
            for v in diff_comb:
                fac *= kwargs[v] - self.ref_point_coordinates[v]
                if fac == 0.: break
            if fac == 0.: continue
            elif fac == 1:
                substrings.append( "+p_C[%i]"%i_comb  )
            else:
                substrings.append( ("%+f"%fac).rstrip('0')+"*p_C[%i]"%i_comb  )

        return "".join( substrings ).lstrip('+')

    def get_weight_func(self, **kwargs):
        '''construct a lambda function that evaluates the weight in terms of the event.p_C coefficient vector using the kwargs as WC
        '''

        # add the arguments from the ref-point 
        self.set_default_args( kwargs )

        terms = []
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            # remove 0 entries
            fac = float( reduce( mul, [ (float(kwargs[v]) - self.ref_point_coordinates[v]) for v in comb ], 1 ) )
            if abs(fac) == 0.: continue
            # store [ ncoeff, factor ]
            terms.append( [ i_comb, fac ] )

        return lambda event, sample: sum( event.p_C[term[0]]*term[1] for term in terms )

    def get_diff_weight_func(self, var, **kwargs):
        '''construct a lambda function that evaluates the diff weight in terms of the event.p_C coefficient vector using the kwargs as WC
        '''

        if var not in self.variables:
            raise ValueError( "Variable %s not in gridpack: %r" % ( var, self.variables ) ) 
        
        # add the arguments from the ref-point 
        self.set_default_args( kwargs )

        terms = []
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            prefac, diff_comb = WeightInfo.differentiate( comb, var )
            if prefac == 0: continue
            # store [ ncoeff, factor ]
            fac = prefac
            for v in diff_comb:
                fac *= kwargs[v] - self.ref_point_coordinates[v]
                if fac == 0.: break
            if fac == 0.: continue 
            terms.append( [ i_comb, fac ] )

        return lambda event, sample: sum( event.p_C[term[0]]*term[1] for term in terms )

    def get_total_weight_yield( self, coeffLists, **kwargs ):
        '''compute yield from a list of coefficients (in the usual order of p_C) using the kwargs as WC
        '''

        # add the arguments from the ref-point 
        self.set_default_args( kwargs )

        # combine lists of weights to one list of weights ( sum_events(w0 + wi*Ci + wij*Ci*Cj) = sum(w0) + sum(wi)*Ci + sum(wij)*Ci*Cj )
        coeffList = [sum(i) for i in zip(*coeffLists)]

        return self.get_weight_yield( coeffList, **kwargs )


    def get_weight_yield( self, coeffList, **kwargs ):
        '''compute yield from a list of coefficients (in the usual order of p_C) using the kwargs as WC
        '''
        # check if coeffList is filled with 0
        if all([ v == 0 for v in coeffList ]): return 0.

        # add the arguments from the ref-point 
        self.set_default_args( kwargs )

        result = 0. 
        for i_comb, comb in enumerate(self.combinations):
            if False in [ v in kwargs for v in comb ]: continue
            if coeffList[i_comb]==0: continue
            # remove 0 entries
            fac = float( reduce( mul, [ (float(kwargs[v]) - self.ref_point_coordinates[v]) for v in comb ], 1 ) )
            if abs(fac) == 0.: continue
            result += coeffList[i_comb] * fac

        return result

    def get_diff_weight_yield( self, vars, coeffList, **kwargs ):
        '''compute diff yield from a list of coefficients (in the usual order of p_C) using the kwargs as WC
        '''

        if type(vars)==type(""): vars = (vars,)

        for var in vars:
            if var not in self.variables:
                raise ValueError( "Variable %s not in gridpack: %r" % ( var, self.variables ) ) 

        # check if coeffList is filled with 0
        if all([ v == 0 for v in coeffList ]): return 0.

        # add the arguments from the ref-point 
        self.set_default_args( kwargs )

        result = 0. 
        for i_comb, comb in enumerate(self.combinations):
            if False in [v in kwargs for v in comb]: continue
            prefac, diff_comb = WeightInfo.differentiate( comb, vars)
            # skip entries which are zero
            if prefac == 0: continue
            if coeffList[i_comb] == 0: continue
            fac = prefac
            for v in diff_comb:
                fac *= kwargs[v] - self.ref_point_coordinates[v]
                # skip entries which are zero
                if fac == 0.: break
            if fac == 0.: continue 
            result += coeffList[i_comb]*fac

        return result

    def get_fisherInformation_matrix( self, coeffList, variables = None, **kwargs ):
        ''' return the fisher information matrix for a single event (coefflist)
        '''

        # check if coeffList is filled with 0
        if all([ v == 0 for v in coeffList ]): return variables, np.zeros( ( len(variables), len(variables) ) )

        # If no argument given, provide all
        if variables is None: variables = self.variables

        # calculate derivatives for all variables
        diff_weight_yield = { var:self.get_diff_weight_yield( var, coeffList, **kwargs ) for var in variables }

        # initialize FI matrix with 1/weight (same for all entries)
        weight_yield = self.get_weight_yield( coeffList, **kwargs ) 
        fi_matrix = np.full( ( len(variables), len(variables) ), 1. / weight_yield if weight_yield != 0 else 0)

        for i, var_i in enumerate(variables):
            for j, var_j in enumerate(variables):
                if fi_matrix[i,j] == 0: continue
                if i<=j: 
                    fi_matrix[i,j] *= diff_weight_yield[var_i] * diff_weight_yield[var_j]
                else:
                    fi_matrix[i,j] = fi_matrix[j,i]

        return variables, fi_matrix

    def get_total_fisherInformation_matrix( self, coeffLists, variables = None, **kwargs ):
        ''' return the full fisher information matrix, sum the FI matrices over all coefflists
        '''

        fi_matrix = np.sum( [ self.get_fisherInformation_matrix( coeffList, variables, **kwargs )[1] for coeffList in coeffLists if not all([ v == 0 for v in coeffList ]) ], 0 )

        return variables, fi_matrix


    def matrix_to_string( self, variables, matrix ):
        ''' return the matrix in a terminal visualization string (print)
        '''

        if variables is None: variables = self.variables

        res = [ ' '.join( map( "{:>9}".format, variables ) ) ]
        for i_line, line in enumerate(matrix.tolist()):
            res.append( ' '.join( map('{:+.2E}'.format, line) + [variables[i_line]] ) )

        return '\n'.join( res ) 

#    def get_diff_fisherInformation_matrix( self, coeffList, variable, variables = None, **kwargs ):
#        ''' return the differentiated fisher information matrix wrt 'variable' for a single event (coefflist)
#        '''
#
#        # If no argument given, provide all
#        if variables is None: variables = self.variables
#
#        ## calculate derivatives for all variables
#        #diff_weight_yield = { var:self.get_diff_weight_yield( var, coeffList, **kwargs ) for var in variables }
#
#        ## initialize FI matrix with 1/weight (same for all entries)
#        #weight_yield = self.get_weight_yield( coeffList, **kwargs ) 
#        #fi_matrix = np.full( ( len(variables), len(variables) ), 1. / weight_yield if weight_yield != 0 else 0)
#
#        #for i, var_i in enumerate(variables):
#        #    for j, var_j in enumerate(variables):
#        #        if fi_matrix[i,j] == 0: continue
#        #        if i<=j: 
#        #            fi_matrix[i,j] *= diff_weight_yield[var_i] * diff_weight_yield[var_j]
#        #        else:
#        #            fi_matrix[i,j] = fi_matrix[j,i]
#
#        #return variables, fi_matrix
#
#
#    def get_total_diff_fisherInformation_matrix( self, coeffLists, variable, variables = None, **kwargs ):
#        ''' return the full differentiated fisher information matrix wrt to 'variable', sum the FI matrices over all coefflists
#        '''
#
#        fi_matrix = np.sum( [ self.get_diff_fisherInformation_matrix( coeffList, variable, variables, **kwargs )[1] for coeffList in coeffLists ], 0 )
#
#        return variables, fi_matrix


    def get_christoffels( self, coeffLists, variables = None): 
        ''' Compute christoffel symbols Gamma^i_jk for coefflist in 
            subspace spanned by variables at the point specified by kwargs

            Gamma^i_jk = 0.5*g^il Sum(1/lambda (dl lambda)(dj lambda)(dk lambda) + 2./lambda^2 (dl lambda)(dj dk  lambda ) )
        '''

        # Restrict to subspace
        _variables = self.variables if variables is None else variables

        # Define a function that accepts an index and a position
        def christoffel_symbols( index, position ):
            ''' Compute christoffel i at position in parameter space'''
            # Metric and Metric-inverse in subspace
            ## Make kwargs dict from position
            kwargs_        = {_variables[i_p]:p for i_p,p in enumerate(position)} 
            metric         = self.get_total_fisherInformation_matrix( coeffLists, variables = _variables, **kwargs_ ) [1]
            metric_inverse = scipy.linalg.inv( metric ) 

            # 3D zeros
            christoffel = np.zeros( (len(_variables), len(_variables) ) )

            for coeffList in coeffLists:
                weight_yield       = self.get_weight_yield( coeffList, **kwargs_ )
                if weight_yield == 0.: continue
                #print "weight_yield", weight_yield
                diff_weight_yield  =  { i_var:self.get_diff_weight_yield( var, coeffList, **kwargs_ ) for i_var, var in enumerate(_variables) }
                diff2_weight_yield = { (i_var_1, i_var_2):self.get_diff_weight_yield( (var_1, var_2), coeffList, **kwargs_ ) for i_var_1, var_1 in enumerate(_variables) for i_var_2, var_2 in enumerate(_variables) }
                for l in xrange(len(_variables)):
                    gil = metric_inverse[index][l]
                    #print "i,l,gil",i,l,gil
                    if gil==0.: continue
                    #print index, gil, dg[index] 
                    for j in range(len(_variables)):
                        for k in range(len(_variables)):
                            d_christoffel_jk = gil*( 0.5/weight_yield*diff_weight_yield[l]*diff_weight_yield[j]*diff_weight_yield[k] + 1./weight_yield**2*diff_weight_yield[l]*diff2_weight_yield[(j,k)] )
                            if j==k:
                                christoffel[j][k] += d_christoffel_jk 
                            elif j>k:
                                christoffel[j][k] += d_christoffel_jk 
                                christoffel[k][j] += d_christoffel_jk

            return christoffel
        return christoffel_symbols 

# Make a list from the bin contents from a histogram that resulted from a 'Draw' of p_C 
def histo_to_list( histo ):
    return [ histo.GetBinContent(i) for i in range( 1, histo.GetNbinsX() + 1 ) ]

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
    # Debug 1 event
    sample.reduceFiles( to = 1 )
    w = WeightInfo(sample.reweight_pkl)
    w.set_order( 2 )

    selection_string = cutInterpreter.cutString('lepSel3-onZ-njet3p-nbjet1p-Zpt0')
    #selection_string = "evt==955001&&run==1&&lumi==9551"

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
    coeff_Z_pt = getCoeffPlotFromDraw( sample, 'Z_pt', [ 20, 0, 500 ], selection_string, weightString='150*ref_lumiweight1fb')
    # Fisher information in x-sec
    #coeff_tot = getCoeffListFromDraw( sample, selection_string, weightString='150*lumiweight1fb')
    coeff_tot = getCoeffListFromDraw( sample, selection_string, weightString='150*ref_lumiweight1fb' )

#    w.get_weight_yield(coeff_Z_pt, ctG=10)


    #print(len(coeff_Z_pt))
    #print w.matrix_to_string(*w.get_total_fisherInformation_matrix(coeff_Z_pt, variables))
    #print np.linalg.eigvals(w.get_total_fisherInformation_matrix(coeff_Z_pt, variables)[1])

    from TTXPheno.Tools.Geodesic import Geodesic

    variables = ('cpQM', 'cpt')

    def phase_space_dict(  point ):
        return { var:val for var,val in zip( sum( [ [v, v+'_dot'] for v in variables ], [] ), point ) }

    christoffel_symbols = w.get_christoffels(  coeff_Z_pt,  variables = variables) 

    initial_point      = (0, 0.)
    initial_derivative = (.002, .002 )

    # How far we want to go in the parameter q
    q_max = 15.
    nq    = 500

    # Initialize Geodesic
    geodesic = Geodesic( initial_point, initial_derivative, christoffel_symbols )

    # Define q values & solve
    q_values = np.linspace(0, q_max, nq+1)
    solution = map( phase_space_dict, geodesic.solve(q_values) )
    # Add the parameter value
    for i_q_value, q_value in enumerate( q_values ):
        solution[i_q_value]['q'] = q_value

    import ROOT
    import array

    gr = ROOT.TGraph(nq, array.array('d', [ y[variables[0]] for y in solution ]),  array.array('d', [ y[variables[1]] for y in solution ]) )
    c1 = ROOT.TCanvas()
    gr.Draw("AC*")
    c1.Print("/afs/hephy.at/user/r/rschoefbeck/www/etc/info_geodesic.png")
