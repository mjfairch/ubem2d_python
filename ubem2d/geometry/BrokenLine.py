import numpy as np
from .Scatter import Scatter

__all__ = ['BrokenLine']

class BrokenLine(Scatter):
    '''
    The BrokenLine class represents a sequence of straight-line segments
    connected end to end.  It adds unit tangent vectors to the underlying
    Scatter object.  If there are N points, then there are N-1 line
    segments, and the ith tangent vector is the unit vector joining the
    ith point to the (i+1)th point.
    '''
    def __init__(self, x, y):
        super().__init__(x,y)
        self.update_broken_line()

    @classmethod
    def join(cls,x1,y1,x2,y2,n):
        return BrokenLine(np.linspace(x1,x2,n),np.linspace(y1,y2,n))

    @property
    def tx(self):
        return self._tx

    @property
    def ty(self):
        return self._ty

    @property
    def xmid(self):
        return self._xmid

    @property
    def ymid(self):
        return self._ymid

    @property
    def edge(self):
        return self._edge

    @property
    def theta(self):
        return self._theta

    @property
    def nedge(self):
        return len(self._x) - 1

    @property
    def perimeter(self):
        return self._perimeter

    @property
    def closed(self):
        tol = 10*np.finfo(float).eps
        return (np.abs(self._x[0] - self._x[-1]) < tol and
            np.abs(self._y[0] - self._y[-1]) < tol)

    @property
    def centroid(self):
        '''
        Override implementation in Scatter, because for a closed body the
        final and initial points are the same and should not be counted twice.
        '''
        if (self.closed):
            return (np.mean(self._x[:-1]), np.mean(self._y[:-1]))
        else:
            return (np.mean(self._x), np.mean(self._y))

    def __len__(self):
        return self.nedge

    def update_broken_line(self):
        dx = self._x[1:] - self._x[:-1]
        dy = self._y[1:] - self._y[:-1]
        self._edge = np.sqrt(dx**2 + dy**2)
        self._tx = dx/self._edge
        self._ty = dy/self._edge
        self._theta = np.arctan2(self._ty,self._tx)
        self._xmid = self._x[:-1] + .5*dx
        self._ymid = self._y[:-1] + .5*dy
        self._perimeter = np.sum(self._edge)
        return self

    def glide(self, g, x0 = None, y0 = None):
        super().glide(g,x0,y0)
        return self.update_broken_line()

    def scale(self, scale):
        super().scale(scale)
        return self.update_broken_line()

    def translate(self, x0, y0):
        super().translate(x0,y0)
        return self.update_broken_line()

    def rotate(self, theta, x0 = None, y0 = None):
        super().rotate(theta,x0,y0)
        return self.update_broken_line()

    def center(self):
        super().center()
        return self.update_broken_line()

    def set_corners(self,x,y):
        super().set_corners(x,y)
        return self.update_broken_line()

    def pop_first(self):
        self._x = self._x[1:]
        self._y = self._y[1:]
        return self.update_broken_line()

    def pop_last(self):
        self._y = self._x[:-1]
        self._y = self._y[:-1]
        return self.update_broken_line()

    def __add__(self, bl):
        '''
        Return a new BrokenLine obtained by adjoining the points in bl,
        also assuemd to be a BrokenLine, to those in this BrokenLine.
        '''
        return BrokenLine(np.concatenate([self._x,bl._x]),
            np.concatenate([self._y,bl._y]))

if __name__ == '__main__':
    bl = BrokenLine.join(0,0,1,1,5) + BrokenLine.join(1,1,2,0,3).pop_first()
    print('x',bl.x)
    print('y',bl.y)
    print('xmid',bl.xmid)
    print('ymid',bl.ymid)
    print('tx',bl.tx)
    print('ty',bl.ty)
    print('edge',bl.edge)
    print('theta',bl.theta)
