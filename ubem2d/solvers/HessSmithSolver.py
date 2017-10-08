from collections import namedtuple
import numpy as np
import numpy.linalg as nla
from ubem2d.panel.BodyInfluence import source_influence_matrices_body
from ubem2d.panel.BodyInfluence import vortex_influence_matrices_body

__all__ = ['solve_hess_smith_body']

def solve_hess_smith_body(uinf, body, At=None, An=None, Bt=None, Bn=None):
    '''
    Solve for the steady flow past an airfoil using the Hess-Smith method.
    The Kutta condition is that there be no pressure difference across the
    trailing edge of the airfoil, which implies (via Bernoulli) that no
    vorticity is being shed from the trailing edge.
    '''
    # Compute influence matrices as needed
    if ((At is not None and An is None) or (At is None and An is not None)):
        raise ValueError('Must specify zero or two influence matrices')
    if ((Bt is not None and Bn is None) or (Bt is None and Bn is not None)):
        raise ValueError('Must specify zero or two influence matrices')
    tx,ty,nx,ny = body.tx, body.ty, body.nx, body.ny
    if (At is None):
        (At,An) = source_influence_matrices_body(body)
    if (Bt is None):
        (Bt,Bn) = vortex_influence_matrices_body(body)

    # Compute the right-hand side of Hess-Smith system.  The last row encodes
    # the Kutta condition; the remaining rows enforce flow tangency.
    n = body.nedge
    rhs = np.zeros(n+1)
    rhs[0:n] = -(uinf[0]*nx + uinf[1]*ny)
    rhs[n] = -(uinf[0]*(tx[0]+tx[-1]) + uinf[1]*(ty[0]+ty[-1]))

    # Compute the Hess-Smith matrix.  The Kutta condition is that no vorticity
    # is shed aft of the airfoil.  Since the flow is tangent to the airfoil's
    # boundary, this amounts to the requirement that the flow velocity at the
    # upper trailing-edge surface and lower trailing-edge surface be equal.
    A = np.zeros((n+1,n+1))
    A[0:n,0:n] = An                         # Flow tangency (source terms)
    A[0:n,n] = np.sum(Bn,1)                 # Flow tangency (vorticity terms)
    A[n,0:n] = At[0,:] + At[-1,:]           # Kutta condition (source terms)
    A[n,n] = Bt[0,:].sum() + Bt[-1,:].sum() # Kutta condition (vorticity terms)

    # Solve Hess-Smith system
    soln = nla.solve(A,rhs)
    sigma = soln[0:-1]  # source strengths along each panel
    gamma = soln[-1]    # common value of vortex strength along all panels

    # Compute flow velocity at panel midpoints
    qt = np.dot(At,sigma) + gamma*np.sum(Bt,1) + uinf[0]*tx + uinf[1]*ty
    qn = np.dot(An,sigma) + gamma*np.sum(Bn,1) + uinf[0]*nx + uinf[1]*ny

    # Compute pressure distribution via the steady Bernoulli equation
    cp = 1. - (qt/nla.norm(uinf))**2

    # Return the pressure distribution, vorticity, and source strenghts
    return namedtuple('soln','sigma,gamma,cp,qt,qn,At,An,Bt,Bn')(sigma,gamma,
        cp,qt,qn,At,An,Bt,Bn)
