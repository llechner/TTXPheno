''' Solve the GR geodesic initial value problem
'''

# standard imports
from math import *
import numpy as np
# scipy 
from scipy.integrate import odeint

class Geodesic:

    def __init__( self,  initial_point, initial_derivative, christoffel_symbols):

        if len(variables)!=len(initial_point) or len(initial_derivative)!=len(christoffel_symbols):
            raise RuntimeError( "Inconsistent dimensions: initial_point %i initial_derivative %i christoffel_symbol %i" \
                    %(len(initial_point), len(initial_derivative), len(christoffel_symbols)) 
                ) 

        # dimension
        self.dim                    = len(initial_point)
        # Expect a tuple of floats
        self.initial_point          = initial_point

        # Initialize start point
        self.q0                     = 0.
        self.initial_conditions     = sum( [ [ initial_point[i], initial_derivative[i] ] for i in range(self.dim)], [])

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


    def solve( self, q_values):

        return odeint(self.rhs, self.initial_conditions, q_values )
        #.set_integrator('vode', method='bdf') 
        #integrator.set_initial_value(self.initial_conditions, q0).set_f_params(2.0).set_jac_params(2.0) 

if __name__ == '__main__':

    # https://en.wikipedia.org/wiki/Schwarzschild_geodesics#Geodesic_equation

    # Schwartzschield
    variables = ( 't', 'r', 'phi' )
    
    # Schwartschild radius
    r_s = 1.
    #                     t  r  phi
    initial_point      = (0, 50, 0)
    initial_derivative = (1., -.1, 0.0027 )

    christoffel_symboles = [\

        lambda position: [
            [ 0, r_s/(2.*position[1]**2*(1.-r_s/position[1])), 0],
            [ r_s/(2.*position[1]**2*(1.-r_s/position[1])), 0, 0],
            [ 0, 0, 0]
        ],
        lambda position: [
            [ r_s/position[1]**3*(position[1]-r_s), 0, 0],
            [ 0, r_s/(2*position[1]*(position[1]-r_s)), 0],
            [ 0, 0, -(position[1]-r_s)]
        ],
        lambda position: [
            [ 0, 0, 0],
            [ 0, 0, 1./position[1]],
            [ 0, 1./position[1], 0]
        ],
    ] 

    def phase_space_dict(  point ):
        return { var:val for var,val in zip( sum( [ [v, v+'_dot'] for v in variables ], [] ), point ) }
    
    # Initialize Geodesic
    geodesic = Geodesic( initial_point, initial_derivative, christoffel_symboles )

    # How far we want to go in the parameter q
    q_max = 80000
    nq    = 1000 

    # Define q values & solve
    q_values = np.linspace(0, q_max, nq+1)
    solution = map( phase_space_dict, geodesic.solve(q_values) )
    # Add the parameter value
    for i_q_value, q_value in enumerate( q_values ):
        solution[i_q_value]['q'] = q_value

    # Make a plot    
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Conserved quantities
    # Angular momentum : r^2 dphi/dq
    angular_momentum = [y['r']**2*y['phi_dot'] for y in solution]
    # Energy: dt/dq * w, w=1-r_s/r
    energy           = [y['t_dot']*(1. - r_s/y['r']) for y in solution]
    #print angular_momentum   
    #print energy
 

    ax = plt.subplot(111, projection='polar')
    ax.plot([y['phi'] for y in solution], [y['r'] for y in solution])
    ax.set_rmax(160)
    ax.set_rticks([20, 50, 100, 150])  # radial ticks
    ax.set_rlabel_position(-22.5)      # get radial labels away from plotted line
    ax.grid(True)

    ax.set_title("A Schwartzschild geodesic (rS = 1)", va='bottom')
    plt.savefig('/afs/hephy.at/user/r/rschoefbeck/www/etc/geodesic5.png')
