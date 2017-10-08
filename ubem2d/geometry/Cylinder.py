import numpy as np
from .Body import Body
from .Orientation import Orientation

__all__ = ['Cylinder']

class Cylinder(Body):
    def __init__(self, r, npanels=50, x0=0, y0=0):
        '''
        Cylinder of arbitrary cross section.  The user provides a lambda
        function, r, which returns the radius as a function of the
        angle theta, i.e., r=r(theta).  For example, the lambda function

        r = lambda t: 1

        would return a circle of radius 1.
        '''
        npanels = int(npanels)
        if (npanels < 3):
            raise ValueError('Number of panels must be at least 3')
        self._x0 = float(x0)
        self._y0 = float(y0)
        theta = np.linspace(0,2*np.pi,npanels+1)
        R = r(theta)
        x = self._x0 + R*np.cos(theta)
        y = self._y0 + R*np.sin(theta)

        super().__init__(x,y,Orientation.CW)

    @property
    def x0(self):
        return self._x0
    
    @property
    def y0(self):
        return self._y0
    
    @property
    def center(self):
        return (self._x0,self._y0)
