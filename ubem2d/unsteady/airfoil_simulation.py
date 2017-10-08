import numpy as np
from ubem2d.aerodynamics.ForceAndMoment import airfoil_cdclcm
from ubem2d.aerodynamics.ForceAndMoment import drag_lift_vectors
from ubem2d.solvers.BasuHancockSolver import BasuHancockSolver

__all__ = ['airfoil_sensor', 'airfoil_stepper']

def airfoil_sensor(uinf, foil, pp, cp, gamma, dalp, dy, dt):
    '''
    Return the aerodynamic coefficients, energies, and bound circulation
    associated with the most recent solution step.

    CT:     Thrust coefficient
    CL:     Lift coefficient
    CM:     Moment coefficient (about specified axis)
    Ein:    The energy consumed by the step
    Eout:   The energy produced by the step
    bcirc:  The bound circulation at the end of the step
    '''
    CD, CL, CM = airfoil_cdclcm(uinf, foil, cp, pp)
    CT = -CD
    D, L = drag_lift_vectors(uinf)
    CFx = -D[0]*CT + L[0]*CL
    CFy = -D[1]*CT + L[1]*CL
    Ein = -(CFy*dy + CM*dalp)
    Eout = -(CFx*np.linalg.norm(uinf)*dt)
    bound_circ = gamma*foil.perimeter
    return CT, CL, CM, Ein, Eout, bound_circ

def airfoil_stepper(foil, motion, pp=0., wake=None, uinf=(1,0)):
    '''
    A generator which returns the kinematic data (time, pitch, heave) and
    the post-solution data given:

    foil:   Airfoil object in its initial orientation
    motion: Generator which returns tuple: time, (pitch,heave)
    pp:     Pitch position (0=LE, .5=midchord, 1=TE)
    uinf:   Background flow
    '''
    if (wake is None):
        wake = ubem.PointVortexWake(eps=1.e-6)
    solver = BasuHancockSolver(foil, wake)
    t0, alp0, y0 = 0, 0, 0
    for t, (alp, y) in motion:
        dt, dalp, dy = t-t0, (alp-alp0)*np.pi/180, (y-y0)*foil.chord
        foil.pitch(dalp, pp)
        foil.heave(dy)
        sig, gam, cp = solver.step(dt, uinf)[0:3]
        CT, CL, CM, Ein, Eout, bcirc = airfoil_sensor(uinf, foil, pp, cp, gam,
            dalp, dy, dt)
        t0, alp0, y0 = t, alp, y
        yield ((t, alp, y), (CT, CL, CM, Ein, Eout, bcirc))
