import numpy as np

__all__ = ['exponential_ramp', 'finite_ramp', 'smooth_ramp']

def exponential_ramp(x):
    '''
    Return the exponential ramp exp(x) defined by

    exp(x<=0) := 0
    exp(x>0)  := exp(-1/x)

    The function is smooth, increasing, and satisfies lim_{x->infty} = 1
    '''
    y = np.zeros(x.shape)
    ip = np.where(x > 0)
    y[ip] = np.exp(-1./x[ip])
    return y

def finite_ramp(x):
    '''
    Return a smooth, increasing ramp which is zero for x <= 0 and 1 for x >= 1.
    The ramp is 'finite' in the sense that it achieves its minimum and maximum
    values on the compact set [0,1].
    '''
    return exponential_ramp(x)/(exponential_ramp(x)+exponential_ramp(1.-x))

def smooth_ramp(x, y0=0., y1=1.):
    '''
    Return a smooth ramp function rising from y0 at x[0] to y1 at x[-1].
    The resulting function has zero derivative at x[0] and x[-1].
    The input x should be increasing.
    '''
    u = (x-x[0])/(x[-1]-x[0])  # map interval [x[0],x[-1]] linearly into [0,1]
    dy = y1-y0                 # ramp amplitude
    return y0 + dy*finite_ramp(u)
