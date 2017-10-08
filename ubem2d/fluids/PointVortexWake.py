import numpy as np
from ubem2d.Errors import SizeMismatchError
from ubem2d.util.arrayify import arrayify

__all__ = ['PointVortexWake']

class PointVortexWake():
    def __init__(self, gam = None, x = None, y = None, eps = 1.e-6):
        if (gam is not None and (x is None or y is None) or
            x is not None and (gam is None or y is None) or
            y is not None and (gam is None or x is None)):
            raise ValueError('Must give all or none of gam, x, y')
        if (gam is not None):
            if (gam.shape != x.shape or x.shape != y.shape):
                raise SizeMismatchError()
            self._gam, self._x, self._y = arrayify(gam, x, y)
        else:
            self._gam, self._x, self._y = [np.zeros(0) for i in range(3)]
        self._eps = eps
    
    @property
    def gam(self):
        return self._gam

    @property
    def x(self):
        return self._x
    
    @property
    def y(self):
        return self._y
    
    @property
    def circulation(self):
        return np.sum(self._gam)
    
    def __len__(self):
        return len(self._gam)
    
    def append(self, gam, x, y):
        self._gam = np.append(self._gam, float(gam))
        self._x = np.append(self._x, float(x))
        self._y = np.append(self._y, float(y))
    
    def velocity(self, x, y):
        '''
        Compute the velocity induced by the wake at the points (x,y), where
        x and y are scalars or numpy arrays of arbitrary positive dimension
        and equal shape.
        '''
        x,y = arrayify(x,y)
        if (x.shape != y.shape):
            raise SizeMismatchError()
        n = len(self._gam)
        u = np.zeros(x.shape)
        v = np.zeros(y.shape)
        for i in range(n):
            dx = x - self._x[i]
            dy = y - self._y[i]
            d2 = dx*dx + dy*dy + self._eps*self._eps  # >= self._eps**2
            u -= self._gam[i]*dy/(2*np.pi*d2)
            v += self._gam[i]*dx/(2*np.pi*d2)
        return (u,v)
    
    def self_velocity(self):
        return self.velocity(self._x, self._y)
    
    def self_advect(self,dt):
        (u,v) = self.self_velocity()
        self._x += u*dt
        self._y += v*dt
    
    def advect(self, vx, vy, dt):
        self._x += vx*dt
        self._y += vy*dt
    
    def vortex_cores(self):
        # Eliminate vortices of zero strength
        in0 = np.where(self.gam != 0)
        mu = self.gam[in0]
        xx = self.x[in0]
        yy = self.y[in0]
        # Find indices where sign of vorticity changes
        isc = np.where(np.diff(np.sign(mu)) != 0)[0]
        # i0/i1 = start/end ranges for regions with constant sign of vorticity
        i0 = np.hstack([[0], 1+isc])
        i1 = np.hstack([1+isc, [len(mu)]])
        # Compute core strengths and core x,y locations, by analogy with
        # the center-of-mass formula with vortex strength in place of mass.
        n = len(i0)
        cs = [np.sum(mu[i0[j]:i1[j]]) for j in range(n)]
        cx = [np.sum(mu[i0[j]:i1[j]]*xx[i0[j]:i1[j]])/cs[j] 
            for j in range(n)]
        cy = [np.sum(mu[i0[j]:i1[j]]*yy[i0[j]:i1[j]])/cs[j] 
            for j in range(n)]
        return cs, cx, cy
