import numpy as np
from ubem2d.Errors import SizeMismatchError
from ubem2d.motion.SE2 import SE2
from ubem2d.util.arrayify import arrayify

__all__ = ['Scatter']

class Scatter():
    def __init__(self, x, y):
        '''
        This class represents a collection of arbitrary points in the 
        Euclidean plane, R^2.  Coordinates are stored in floating-point
        format.
        '''
        self.set_corners(x,y)
    
    @property
    def x(self):
        return self._x
    
    @property
    def y(self):
        return self._y
    
    @property
    def diameter(self):
        '''
        The diameter of a collection of points is the maximum value of the 
        distance between any two points.  Computing this is O(n^2), where n 
        is the number of corners, so it is computed only when needed.
        '''
        if (self._diameter < 0):
            d2 = 0.  # square of the diameter
            n = len(self)
            for i in range(n):
                for j in range(i+1,n):
                    dx = self._x[j]-self._x[i]
                    dy = self._y[j]-self._y[i]
                    dij2 = dx*dx + dy*dy
                    if (dij2 > d2):
                        d2 = dij2
            self._diameter = np.sqrt(d2)
        return self._diameter
    
    @property
    def centroid(self):
        return (np.mean(self._x),np.mean(self._y))
    
    def glide(self, g, x0 = None, y0 = None):
        '''
        Apply the rigid Euclidean motion g to all points of this body, 
        with rotations taken about the axis (x0,y0).  If the axis is
        not specified, it is taken to be the centroid.
        '''
        diam = self._diameter # Store diameter: unaffected by rigid motions
        if ((x0 is None and y0 is not None) or (x0 is not None and y0 is None)):
            raise ValueError('Must specify both axis coordinates, or neither')
        if (x0 is None):
            (x0,y0) = self.centroid
        (xx,yy) = g.map_point(self._x,self._y,x0,y0)
        self.set_corners(xx,yy)
        self._diameter = diam # Restore diameter
        return self
    
    def scale(self, scale):
        '''
        Scales the cloud of points about its centroid by the given amount.
        The centroid remains unchanged, and the diameter scales by |scale|
        '''
        diam = self._diameter
        (x0,y0) = self.centroid
        xx = (self._x-x0)*scale + x0
        yy = (self._y-y0)*scale + y0
        self.set_corners(xx,yy)
        self._diameter = np.abs(scale)*diam
        return self
    
    def translate(self, x0, y0):
        return self.glide(SE2(0,x0,y0))
    
    def rotate(self, theta, x0 = None, y0 = None):
        return self.glide(SE2(theta,0,0),x0,y0)
    
    def center(self):
        '''
        Centers the body at the origin, i.e. translates the body so that its
        centroid coincides with the origin.
        '''
        (x0,y0) = self.centroid
        return self.translate(-x0,-y0)
    
    def set_corners(self,x,y):
        x,y = arrayify(x,y)
        if (len(x) != len(y)):
            raise SizeMismatchError()
        self._x = x
        self._y = y
        # A negative diameter signifies it must be recomputed
        self._diameter = -1.
        return self
    
    def __len__(self):
        return len(self._x)

if __name__ == '__main__':
    x = [0,1]
    y = [0,0]
    sc = Scatter(x,y)
    sc.glide(SE2(np.pi/2,0,0),0,0)
    print('x',sc.x)
    print('y',sc.y)
    print('Diameter', sc.diameter)
    print('Centroid', sc.centroid)
