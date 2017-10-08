import math
from ubem2d.util.coroutines import arithmetic_stepper

__all__ = ['time_step', 'time_stepper']

def time_step(res, tau, Tfast, T = None):
    '''
    Return the largest time step dt satisfying the following constraints:
    1. There are at least res (resolution) time steps at the fast time scale,
       given by Tfast
    2. If the motion is period with period T, the time step dt will divide T

    See the associated documenation on computation of an appropriate time step.
    '''
    steps_per_cycle = None
    if (T is None):
        dt = min(tau, Tfast)/res
    else:
        steps_per_cycle = math.ceil(res*T/min(tau, Tfast))
        dt = T/steps_per_cycle
    return (dt, steps_per_cycle)

def time_stepper(dt, stepmax=None, tmax=None, t=0.):
    if (tmax is None and stepmax is None):
        return arithmetic_stepper(dt, t)
    elif (tmax is not None and stepmax is not None):
        return arithmetic_stepper(dt, t, lambda step, t: (step >= stepmax or 
            t >= tmax))
    elif (stepmax is None and tmax is not None):
        return arithmetic_stepper(dt, t, lambda step, t: t >= tmax)
    else:
        return arithmetic_stepper(dt, t, lambda step, t: step >= stepmax)
