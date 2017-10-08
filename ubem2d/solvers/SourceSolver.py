from collections import namedtuple
import numpy as np
import numpy.linalg as nla
from ubem2d.panel.PanelInfluence import source_influence_matrices

__all__ = ['solve_source_body']

def solve_source_body(Uinf, body, At = None, An = None):
    '''
    Solve for the steady flow past a body using the source panel method.

    Along the panels of each body are distributed sources of constant 
    strength per unit length along each panel but varying from panel to panel.
    The strengths are solved by the requirement that the net flow, due to
    the background flow and all the source distributions, is tangent to each
    panel at the panel's midpoint.  This produces a linear system of equations,
    with one equation for each panel.

    This method returns a named tuple with the following fields:
    sigma: the source strengths for each panel
    cp: the pressure coefficient at each panel (computed via Bernoulli)
    qt: the net tangential flow speed at panel midpoints
    qn: the net normal flow speed at panel midpoints
    At: the tangential influence matrix
    An: the normal influence matrix
    '''
    # Check if influence matrices are provided; create them if not
    if ((At is None and An is not None) or (At is not None and An is None)):
        raise ValueError('Must specify zero or two influence matrices')
    if (An is None or At is None):
        (At,An) = source_influence_matrices(body.x[:-1],body.y[:-1],
            body.tx, body.ty, body.nx, body.ny, body.edge)

    rhs = -(Uinf[0]*body.nx + Uinf[1]*body.ny)
    sigma = nla.solve(An,rhs)
    qt = At.dot(sigma) + Uinf[0]*body.tx + Uinf[1]*body.ty
    qn = An.dot(sigma) + Uinf[0]*body.nx + Uinf[1]*body.ny
    cp = 1. - (qt/nla.norm(Uinf))**2

    return namedtuple('soln',['sigma','cp','qt','qn','At','An'])(sigma,cp,
        qt,qn,At,An)

if __name__ == '__main__':
    from ubem2d.geometry.CircularCylinder import CircularCylinder
    body = CircularCylinder(50)
    Uinf = (1,0)
    soln = solve_source_body(Uinf,body)
    print(nla.norm(soln.qn))

    theta = np.arctan2(body.ymid,body.xmid)
    exact_theta = np.linspace(-np.pi,np.pi,100)
    exact_cp = 1.-4.*np.sin(exact_theta)**2

    import matplotlib.pyplot as plt
    I = theta.argsort()
    plt.plot(theta[I], soln.cp[I], 'o')
    plt.plot(exact_theta, exact_cp)
    plt.xlabel('Angle (radians)')
    plt.ylabel('Pressure coefficient')
    plt.legend(['Panel method','Theoretical solution'])
    plt.title('Pressure distribution for uniform flow past a cylinder')
    plt.grid(True)
    plt.show()
