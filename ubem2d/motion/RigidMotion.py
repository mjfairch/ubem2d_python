import numpy as np
from .SE2 import SE2
from ubem2d.util.arrayify import arrayify

__all__ = ['RigidMotion']

class RigidMotion():
    '''
    This class represents a parameterized path in the Lie group SE(2) of
    rigid motions of the Euclidean plane.
    '''
    def __init__(self, ftheta = None, fx = None, fy = None):
        '''
        ftheta: callable which returns theta, given the time
        fx: callable which returns x, given the time
        fy: callable which returns y, given the time
        '''
        if (ftheta is None):
            ftheta = lambda t: np.zeros(arrayify(t).shape)
        if (fx is None):
            fx = lambda t: np.zeros(arrayify(t).shape)
        if (fy is None):
            fy = lambda t: np.zeros(arrayify(t).shape)

        if (not callable(ftheta) or not callable(fy) or not callable(fy)):
            raise ValueError('Arguments must be callable')

        self._ftheta = ftheta
        self._fx = fx
        self._fy = fy
    
    @property
    def ftheta(self):
        return self._ftheta
    
    @property
    def fx(self):
        return self._fx
    
    @property
    def fy(self):
        return self._fy
    
    def __call__(self, t):
        scalar_input = False
        if (np.isscalar(t)):
            scalar_input = True
            t = np.array([t],dtype=np.float64)
        if (type(t) in [list,tuple]):
            t = np.array(t)
        if (t.ndim != 1):
            raise ValueError('Time parameter must be a 1-dimensional array')
        if (len(t) == 0):
            raise ValueError('Must provide at least one time value')

        n = len(t)
        th = self._ftheta(t)
        x = self._fx(t)
        y = self._fy(t)
        g = n*[None]
        for i in range(n):
            g[i] = SE2(th[i],x[i],y[i])

        if (scalar_input):
            return g[0]
        else:
            return g
