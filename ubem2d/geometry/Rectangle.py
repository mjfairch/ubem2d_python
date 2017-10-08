import numpy as np
from .Body import Body
from .Orientation import Orientation

__all__ = ['Rectangle']

class Rectangle(Body):
    def __init__(self, width=1, height=1, nx=10, ny=10, x0=0, y0=0, ctr=True):
        '''
        width: length of top/bottom edges (i.e. edges parallel to x-axis)
        height: length of left/right edges (i.e. edges parallel to y-axis)
        nx: number of panels on the top/bottom edges
        ny: number of panels on the left/right edges
        x0: x-coordinate of rectangle's center or lower-left corner
        y0: y-coordinate of rectangle's center or lower-left corner
        ctr: if True, then (x0,y0) gives the coordinates of the
             center of the rectangle.  Otherwise, it gives the coordinates of
             the lower-left corner.
        '''
        width = float(width)
        height = float(height)
        if (width < 0 or height < 0):
            raise ValueError('Width and height must be positive')
        nx = int(nx)
        ny = int(ny)
        if (nx < 1 or ny < 1):
            raise ValueError('Number of edge panels must be at least 1')
        x0 = float(x0)
        y0 = float(y0)

        # Begin with lower-left corner at (x0,y0) and proceed counterclockwise
        xx = np.linspace(0,width,nx+1)   # coordinates along top/bottom edges
        yy = np.linspace(0,height,ny+1)  # coordinates along left/right edges
        x = np.concatenate([xx, ny*[xx[-1]], np.flipud(xx[0:-1]), ny*[xx[0]]])
        y = np.concatenate([(nx+1)*[yy[0]], yy[1:], nx*[yy[-1]], 
            np.flipud(yy[:-1])])
        if (ctr):
            x -= width/2
            y -= height/2
        x += x0
        y += y0

        super().__init__(x,y,Orientation.CW)

        self._width = width
        self._height = height
    
    @property
    def width(self):
        return self._width
    
    @property
    def height(self):
        return self._height
    
    @property
    def center(self):
        return self.centroid

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    width = float(input('Enter width: '))
    height = float(input('Enter height: '))
    nx = int(input('Enter number of panels on top/bottom edges: '))
    ny = int(input('Enter number of panels on left/right edges: '))
    rect = Rectangle(width, height, nx, ny)
    n = len(rect.edge)
    for i in range(n):
        print('{}\t{}'.format(rect.x[i],rect.y[i]))
    print('\nWidth:', rect.width)
    print('Height:', rect.width)
    print('Center:', rect.center)
    print('Closed:', rect.closed)
    print('Centroid:', rect.centroid)
    print('Perimeter:', rect.perimeter)
    print('Diameter:', rect.diameter)

    plt.plot(rect.x,rect.y,'.-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Rectangle with {} panels'.format(rect.nedge))
    plt.axis('equal')
    plt.grid(True)
    plt.show()
