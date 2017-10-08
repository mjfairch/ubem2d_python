import numpy as np
from ubem2d.geometry.BrokenLine import BrokenLine
from ubem2d.geometry.Orientation import Orientation
from ubem2d.math.turning_angle import turning_angle

__all__ = ['Body']

class Body(BrokenLine):
    '''
    The Body class extends a BrokenLine by adding normal vectors along each
    line segment.  The normal vectors are obtained by rotating the unit tangent
    vectors either clockwise or counterclockwise.  If the sense of rotation is
    not specified, the code attempts to determine it automatically by computing
    the turning angle of the underlying broken line.
    '''
    def __init__(self, x, y, rot = None):
        '''
        rot: sense in which tangent vectors are rotated into normal vectors.
        '''
        if (rot is None):
            ta = turning_angle(x, y, True)
            if (ta > 0):
                rot = Orientation.CW
            elif (ta < 0):
                rot = Orientation.CCW
            else:
                raise ValueError('Cannot determine normal vector orientation')
        if (rot is not Orientation.CCW and rot is not Orientation.CW):
            raise ValueError('Invalid orientation')
        self._rot = rot
        super().__init__(x,y)
        self.update_body()
    
    @property
    def nx(self):
        return self._nx
    
    @property
    def ny(self):
        return self._ny
    
    @property
    def rot(self):
        return self._rot
    
    @property
    def beta(self):
        return self._beta
    
    def update_body(self):
        if (self._rot == Orientation.CCW):
            self._nx = -self._ty
            self._ny = self._tx
        else:
            self._nx = self._ty
            self._ny = -self._tx
        self._beta = np.arctan2(self._ny,self._nx)
        return self
    
    def glide(self, g, x0 = None, y0 = None):
        super().glide(g,x0,y0)    
        return self.update_body()

    def scale(self, scale):
        super().scale(scale)
        return self.update_body()
    
    def translate(self, x0, y0):
        return super().translate(x0,y0)
    
    def rotate(self, theta, x0 = None, y0 = None):
        super().rotate(theta,x0,y0)
        return self.update_body()
    
    def center(self):
        return super().center()
    
    def set_corners(self,x,y):
        super().set_corners(x,y)
        return self.update_body()
