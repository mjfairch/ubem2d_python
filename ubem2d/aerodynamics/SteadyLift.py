import math
import numpy as np
import numpy.linalg as nla
from .ForceAndMoment import airfoil_cdclcm
from ubem2d.panel.BodyInfluence import source_influence_matrices_body
from ubem2d.panel.BodyInfluence import vortex_influence_matrices_body
from ubem2d.geometry.Orientation import Orientation
from ubem2d.solvers.HessSmithSystem import HessSmithSystem
from ubem2d.aerodynamics.ForceAndMoment import airfoil_cdclcm

__all__ = ['steady_lift_model']

def steady_lift_model(foil, aoa0, aoa1, n):
    '''
    Fit a linear model CL = CL0 + m*(aoa) for constants CL0 and m, which
    represents the lift coefficient CL as a function of angle of attack
    (in degrees) for the given airfoil.  The lift coefficient is computed
    at n equally spaced angles in the interval [aoa0,aoa1], and least-squares
    fit is used to determine the constants CL0 and m.
    '''
    (At,An) = source_influence_matrices_body(foil)
    (Bt,Bn) = vortex_influence_matrices_body(foil)
    angles = np.linspace(aoa0,aoa1,n)
    A = np.zeros((n,2))
    A[:,0] = 1.
    A[:,1] = np.linspace(aoa0,aoa1,n)
    b = np.zeros((n,1))
    for i in range(n):
        aoa_rad = A[i,1]*math.pi/180
        if (foil.pitch_up is Orientation.CW):
            # Airfoil points leftward and flow is to the right
            uinf = (math.cos(aoa_rad), math.sin(aoa_rad))
        else:
            # Airfoil points rightward and flow is to the left
            uinf = (-math.cos(aoa_rad), math.sin(aoa_rad))
        sys = HessSmithSystem(foil)
        soln = sys.solve(uinf)
        cp = sys.pressure_self(uinf, soln)[0]
        CD,CL,CM = airfoil_cdclcm(uinf, foil, cp)
        b[i] = CL
    x = nla.lstsq(A,b,rcond=None)[0][:,0]
    return (x[0],x[1],A[:,1],b)
