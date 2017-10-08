import os
import math
import numpy as np
from ubem2d.geometry.Body import Body
from ubem2d.geometry.Orientation import Orientation
from ubem2d.motion.SE2 import SE2
from ubem2d.motion.RigidMotion import RigidMotion
from ubem2d.math.turning_angle import turning_angle
from ubem2d.util.arrayify import arrayify

__all__ = ['Airfoil']

class Airfoil(Body):
    '''
    This class represents a symmetric or cambered airfoil.
    '''

    def __init__(self, x, y, le, pitch_up = None):
        '''
        x,y: Coordinates of panel corners; vectors must be of equal length.
        le:  Index of the leading-edge corner (trailing-edge corner must have
             index zero according to assumption 1 above.)
        pitch_up: Sense of the pitch-up orientation for a flow moving over 
             the leading edge and heading towards the trailing edge.

        The corners x,y must be ordered from the trailing edge to the leading
        edge and back again, although whether they are ordered clockwise or
        counterclockwise is immaterial.

        If pitch_up is not specified, it is automatically determined by 
        assuming that the airfoil is oriented horizontally and that the
        camber, if present, is concave down.
        '''
        ta = turning_angle(x,y,True)
        tvec_to_nvec = Orientation.CW
        if (ta < 0):
            tvec_to_nvec = Orientation.CCW
            raise NotImplementedError('Clockwise airfoil ordering not yet ' \
                'supported; it breaks unsteady solver.')
        if (le == 0 or le >= len(x)):
            raise ValueError('Bad leading-edge index')
        self._le = le              # index of leading-edge corner
        if (pitch_up is None):
            if (x[le] < x[0]):     # leading edge points leftward
                pitch_up = Orientation.CW
            elif (x[le] > x[0]):   # leading edge points rightward
                pitch_up = Orientation.CCW
            else:
                raise ValueError('Unable to determine pitch-up orientation')
        self._pitch_up = pitch_up  # pitch-up orientation (CCW/CW)
        super().__init__(x,y,tvec_to_nvec)
    
    @classmethod
    def load_from_file(cls, path):
        '''
        Read airfoil coordinates from a text file.  The first line of the file
        should contain two integers, the first encoding the number of corners,
        and the second giving the leading-edge index.  The remaining rows
        contain the x,y coordinates of the corners.
        '''
        if (not os.path.exists(path)):
            raise ValueError('Bad path: {}'.format(path))
        with open(path,'r') as fid:
            (n,le) = [int(x) for x in fid.readline().split()]
            x, y = np.zeros(n), np.zeros(n)
            for i in range(n):
                x[i], y[i] = [float(x) for x in fid.readline().split()]
            fid.close()
        return Airfoil(x,y,le)
    
    def export_to_file(self, path):
        '''
        Save an airfoil to a text file.  The file format is the same as
        described in load_from_file().
        '''
        with open(path,'w') as fid:
            n = len(self)
            fid.write('{}\t{}\n'.format(n,self._le))
            for i in range(n):
                fid.write('{}\t{}\n'.format(self._x[i],self._y[i]))
            fid.close()
    
    @property
    def le(self):
        return self._le
    
    @property
    def pitch_up(self):
        return self._pitch_up
    
    @property
    def leading_edge(self):
        return (self.x[self.le],self.y[self.le])
    
    @property
    def trailing_edge(self):
        return ([self.x[0],self.y[0]])
    
    @property
    def chord(self):
        (xle,yle) = self.leading_edge
        (xte,yte) = self.trailing_edge
        return math.sqrt((xte-xle)**2+(yte-yle)**2)
    
    def chord_point(self, a):
        '''
        Returns the point LE + a*(TE-LE) on the chord line, where LE/TE
        signify the leading/trailing edges, respectively.  So, for example,
        chord_point(0) is the leading edge, chord_point(.25) is the quarter-,
        chord, chord_point(.5) is the midchord, chord_point(1) is the trailing
        edge, etc.
        '''
        xle,yle = self.leading_edge
        xte,yte = self.trailing_edge
        return (xle+a*(xte-xle), yle+a*(yte-yle))

    def pitch(self, alpha, a = 0.):
        '''
        Pitch by alpha radians about an axis located at self.chord_point(a).
        '''
        (xp,yp) = self.chord_point(a)
        if (self._pitch_up == Orientation.CCW):
            return self.glide(SE2(alpha),xp,yp) 
        else:
            return self.glide(SE2(-alpha),xp,yp) 

    def surge(self, dx):
        '''
        Translate by dx in the x direction.
        '''
        return self.glide(SE2(0,dx,0))

    def heave(self, dy):
        '''
        Translate by dy in the y direction.
        '''
        return self.glide(SE2(0,0,dy))

    def glide(self, g, x0 = None, y0 = None):
        super().glide(g,x0,y0)    
        return self

    def scale(self, scale):
        super().scale(scale)
        return self
    
    def translate(self, x0, y0):
        super().translate(x0,y0)
        return self
    
    def rotate(self, theta, x0 = None, y0 = None):
        super().rotate(theta,x0,y0)
        return self
    
    def center(self):
        super().center()
        return self
    
    def set_corners(self,x,y):
        super().set_corners(x,y)
        return self
