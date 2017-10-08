import math
import numpy as np

__all__ = ['SE2']

class SE2():
    '''
    This class represents an element of the Euclidean Lie group SE(2) of
    orientation-preserving isometries of the Euclidean plane.  It follows
    from the Mazur-Ulam theorem that such an isometry represents a
    rotation about the origin followed by a translation, and hence local
    coordinates for SE(2) are given by the rotation angle theta, and the
    horizontal and vertical translations x and y.
    '''
    def __init__(self, theta = 0., x = 0., y = 0.):
        self.theta = float(theta)
        self.x = float(x)
        self.y = float(y)
    
    @classmethod
    def id(cls):
        return SE2(0,0,0)
    
    @classmethod
    def update(cls,g0,g1):
        return g1*(g0.inv())

    @property
    def theta(self):
        return self._theta
    
    @theta.setter
    def theta(self,value):
        self._theta = value
        self._cth = math.cos(value)
        self._sth = math.sin(value)
    
    @property
    def cth(self):
        return self._cth
    
    @property
    def sth(self):
        return self._sth
    
    def inv(self):
        '''
        Return the inverse of this element of SE(2).
        '''
        return SE2(-self.theta, 
            -self.x*self.cth - self.y*self.sth,
            self.x*self.sth - self.y*self.cth)
    
    def mat(self):
        '''
        Return the matrix representation of this element of SE(2).
        '''
        return np.array([
            [self._cth,-self._sth,self.x],
            [self._sth,self._cth,self.y],
            [0.,0,1.]])
    
    def r3(self):
        '''
        Return a 3-tuple of the local coordinates for this element of SE(2).
        '''
        return (self._theta, self.x, self.y)
    
    def map_vector(self,vx,vy):
        return (self.cth*vx - self.sth*vy, self.sth*vx + self.cth*vy)

    def map_point(self, x, y, x0 = 0., y0 = 0.):
        '''
        Apply the transformation represented by this element of SE(2) to the
        given points (x,y), allowing for rotations about an arbitrary axis
        given in (x0,y0).
        '''
        if (x0 == 0 and y0 == 0):
            return (self.cth*x - self.sth*y + self.x,
                self.sth*x + self.cth*y + self.y)
        else:
            return (self.cth*(x-x0) - self.sth*(y-y0) + self.x + x0,
                self.sth*(x-x0) + self.cth*(y-y0) + self.y + y0)
    
    def __eq__(self, g):
        return (self._theta == g._theta and self.x == g.x and
            self.y == g.y)

    def __repr__(self):
        return '({},{},{})'.format(self.theta,self.x,self.y)
    
    def __mul__(self, g):
        '''
        Multiplication of two elements of SE(2) corresponds to composition of
        two successive rigid motions.  In particular, given two elements
        g,h of SE(2), the product g*h represents the composition h followed
        by g.  Note that translations commute with translations, and 
        rotations commute with rotations, but translations and rotations
        do not commute with one another.  For instance, if:

        rot = SE2(pi/4,0,0),
        tr = SE2(0,1,0),

        then:

        rot*tr == SE2(pi/4,sqrt(2),sqrt(2)),
        tr*rot == SE2(pi/4,1,0)
        '''
        return SE2(self._theta + g._theta,
            g.x*self._cth - g.y*self._sth + self.x,
            g.x*self._sth + g.y*self._cth + self.y)
