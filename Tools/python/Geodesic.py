''' Solve the GR geodesic initial value problem
'''

# standard imports
from math import *
import numpy as np

class Geodesic:

    def __init__( self, christoffel_symbols ):

        # dimension
        self.dim                    = len(christoffel_symbols)
        # Initialize start point
        self.q0                     = 0.

        # Expect a tuple of functions of a point returning Gamma^i_{jk} where i is the tuple index
        # christoffel_symbols[0](0,1) should evaluate to Gamma^0_{0,1}
        self.christoffel_symbols    = christoffel_symbols

        # The geodesic equation written in first order is:
        # That's 2*dim equations which we write as dx1/dq = p1, dp1/dq = ..., dx2/dq = p2, dp2/dq = ..., ...
        # d/dq x^i = p^i            -> even places 
        # d/dq p^i = - Christoffel^i_(jk) p^j p^k -> odd places

        # Compute r.h.s of 1st order system
        def rhs( phase_space_point, q): 
            # Initialize with zeros
            result = np.zeros( 2*self.dim )
            # Specify coordinate derivatives
            for i in range(self.dim):
                # d/dq x^i = p^i -> the even places are for the coordinates 
                result[2*i]   = phase_space_point[2*i + 1] 
                # d/dq p^i = - Christoffel^i_(jk) p^j p^k -> the odd places are for the derivatives of the coordinates )
                christoffel_i = self.christoffel_symbols[i](phase_space_point[::2])
                result[2*i+1] = - sum( phase_space_point[2*j + 1]*phase_space_point[2*k + 1]*christoffel_i[j][k] for j in range(self.dim) for k in range(self.dim )) 

            return result

        self.rhs = rhs


    def solve_ivp( self, initial_point, initial_derivative, q_values):
        ''' Solve the initial value problem '''
        if len(initial_point)!=len(self.christoffel_symbols) or len(initial_derivative)!=len(self.christoffel_symbols):
            raise RuntimeError( "Inconsistent dimensions: initial_point %i initial_derivative %i christoffel_symbol %i" \
                    %(len(initial_point), len(initial_derivative), len(christoffel_symbols)) 
                ) 
        # Expect a tuple of floats
        initial_conditions     = sum( [ [ initial_point[i], initial_derivative[i] ] for i in range(self.dim)], [])

        from scipy.integrate import odeint
        return odeint(self.rhs, initial_conditions, q_values )
        #.set_integrator('vode', method='bdf') 
        #integrator.set_initial_value(self.initial_conditions, q0).set_f_params(2.0).set_jac_params(2.0) 

    #def solve_bvp( self, initial_point, end_point, q_values):
    #    ''' Solve the boundary value problem '''
    #    def bc( ya, yb ):
    #        # Return the positions as boundary conditions
    #        return np.array([ya[2*i]-initial_point[i] for i in range(self.dim)] + [yb[2*i]-end_point[i] for i in range(self.dim)])
    #    
    #    from scipy.integrate import solve_bvp
    #    return solve_bvp( lambda x,y: self.rhs(y,x), bc, q_values )        

def phase_space_dict(  point ):
    return { var:val for var,val in zip( sum( [ [v, v+'_dot'] for v in variables ], [] ), point ) }

if __name__ == '__main__':
    
    # Schwartschild radius
    r_s = 1.
    variables = ( 't', 'r', 'phi' )

    # Schwartzschild from Metric
    ## https://en.wikipedia.org/wiki/Schwarzschild_geodesics#Geodesic_equation
    metric_tensor = lambda position: [\
        [ (1.-r_s/position[1]), 0, 0],
        [ 0, -1./(1.-r_s/position[1]), 0],
        [ 0, 0, -position[1]**2 ],
        #[ 0, 0, 0, -position[1]**2*sin(position[2])**2 ],
    ]
    metric_tensor_inverse = lambda position: [\
        [ 1./(1.-r_s/position[1]), 0, 0],
        [ 0, -(1.-r_s/position[1]), 0],
        [ 0, 0, -1./position[1]**2 ],
        #[ 0, 0, 0, -1./(position[1]**2*sin(position[2])**2) ],
    ]
    
    # Gamma^i_jk = 1/2*g^il( - d_l g_jk + d_j g_lk + d_k g_jl ) 

    # Deltas used for differentiating
    diff = 0.00001
    delta_diff = np.array( [ diff]*len(variables) ) 

    def christoffel(index):
        delta_diff_vec = [ [0]*len(variables) for i in range(len(variables)) ]
        for i in range(len(variables)):
            delta_diff_vec[i][i] = delta_diff[i]
        def __christoffel( position ):
            # differentiate metric
            dg = [(np.array( metric_tensor(position + delta_diff_vec[l])) - np.array(metric_tensor(position)))/delta_diff[l] for l in range(len(variables)) ] 
            result = [ [ 0 for i in range(len(variables))] for j in range(len(variables)) ]
            for l in xrange(len(variables)):
                gil = metric_tensor_inverse(position)[index][l]
                if gil==0.: continue
                #print index, gil, dg[index] 
                for j in range(len(variables)):
                    for k in range(len(variables)):
                        #print index, l, j, k, 0.5*gil, -dg[l][j][k],  0.5*gil*(-dg[l][j][k]), 0.5*gil*dg[j][k][l], 0.5*gil*dg[k][j][l]
                        result[j][k] += 0.5*gil*( -dg[l][j][k] + dg[j][k][l] + dg[k][j][l] )
            return result
        return __christoffel


    christoffel_symboles_differentiated = [ christoffel(i) for i in range(len(variables)) ]

    # Initialize Geodesic
    geodesic_differentiated = Geodesic( christoffel_symboles_differentiated )

    ## Schwartzschild from Christoffel
    ## https://en.wikipedia.org/wiki/Schwarzschild_geodesics#Geodesic_equation
    christoffel_symboles_exact = [\

        lambda position: [
            [ 0, r_s/(2.*position[1]**2*(1.-r_s/position[1])), 0],
            [ r_s/(2.*position[1]**2*(1.-r_s/position[1])), 0, 0],
            [ 0, 0, 0]
        ],
        lambda position: [
            [ r_s/(2*position[1]**3)*(position[1]-r_s), 0, 0],
            [ 0, r_s/(2*position[1]*(position[1]-r_s)), 0],
            [ 0, 0, -(position[1]-r_s)]
        ],
        lambda position: [
            [ 0, 0, 0],
            [ 0, 0, 1./position[1]],
            [ 0, 1./position[1], 0]
        ],
    ] 

    # Initialize Geodesic
    geodesic_exact          = Geodesic( christoffel_symboles_exact )

    ## Boundary value problem
    ##                     t  r  phi
    #initial_point      = (0, 50., 0)
    #end_point          = (1, 30, pi/2.)
    ## How far we want to go in the parameter q
    #q_max = 10000
    #nq    = 50000
    ## Define q values & solve
    #q_values = np.linspace(0, q_max, nq+1)

    #y_initial_values

    #solution_exact = map( phase_space_dict, geodesic_exact.solve_bvp(initial_point, end_point, q_values) )

    # Initial value problem
    #                     t  r  phi
    initial_point      = (0, 50., 0)
    initial_derivative = (1., -.01, 3./initial_point[1]**2 )

    # How far we want to go in the parameter q
    q_max = 10000
    nq    = 50000
    # Define q values & solve
    q_values = np.linspace(0, q_max, nq+1)
    solution_exact = map( phase_space_dict, geodesic_exact.solve_ivp(initial_point, initial_derivative, q_values) )
    # Add the parameter value
    for i_q_value, q_value in enumerate( q_values ):
        solution_exact[i_q_value]['q'] = q_value
    solution_diff = map( phase_space_dict, geodesic_differentiated.solve_ivp(initial_point, initial_derivative, q_values) )
    # Add the parameter value
    for i_q_value, q_value in enumerate( q_values ):
        solution_diff[i_q_value]['q'] = q_value

    # Make a plot    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Conserved quantities
    # Angular momentum : r^2 dphi/dq
    #angular_momentum = [y['r']**2*y['phi_dot'] for y in solution]
    # Energy: dt/dq * w, w=1-r_s/r
    #energy           = [y['t_dot']*(1. - r_s/y['r']) for y in solution]
    #print angular_momentum   
    #print energy
 
    ax = plt.subplot(111, projection='polar')
#    ax.plot([y['phi'] for y in solution_exact ], [y['r'] for y in solution_exact ])
#    ax.plot([y['phi'] for y in solution_diff ], [y['r'] for y in solution_diff ])
    ax.set_rmax(60)
    ax.set_rticks([5, 10, 20, 30, 40, 50])  # radial ticks
    ax.set_rlabel_position(-22.5)      # get radial labels away from plotted line
    ax.grid(True)

    ax.set_title("A Schwartzschild geodesic (rS = 1)", va='bottom')
    plt.savefig('/afs/hephy.at/user/r/rschoefbeck/www/etc/geodesic_bvp.png')
