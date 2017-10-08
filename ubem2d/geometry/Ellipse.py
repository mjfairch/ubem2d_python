import numpy as np
from .Body import Body
from .Orientation import Orientation

__all__ = ['Ellipse']

class Ellipse(Body):
    def __init__(self, npanels=50, semi_x=1, semi_y=.5, x0=0, y0=0):
        npanels = int(npanels)
        if (npanels < 4):
            raise ValueError('Number of panels must be at least 4')
        self._x0 = float(x0)
        self._y0 = float(y0)
        self._semi_x = float(semi_x)
        self._semi_y = float(semi_y)
        if (semi_x <= 0 or semi_y <= 0):
            raise ValueError('Ellipse must have positive semi-axes')
        t = np.linspace(0,2*np.pi,npanels+1)
        x = self._x0 + semi_x*np.cos(t)
        y = self._y0 + semi_y*np.sin(t)
        super().__init__(x,y,Orientation.CW)
    
    @property
    def semi_x(self):
        return self._semi_x
    
    @property
    def semi_y(self):
        return self._semi_y
    
    @property
    def x0(self):
        return self._x0
    
    @property
    def y0(self):
        return self._y0
    
    @property
    def center(self):
        return (self._x0,self._y0)

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    n = int(input('Enter number of panels: '))
    a = float(input('Enter semi-x axis: '))
    b = float(input('Enter semi-y axis: '))
    ell = Ellipse(n,a,b)
    for i in range(n+1):
        print('{}\t{}'.format(ell.x[i],ell.y[i]))
    print('\nSemi-x axis:', ell.semi_x)
    print('Semi-y axis:', ell.semi_y)
    print('Center:', ell.center)
    print('Closed:', ell.closed)
    print('Centroid:', ell.centroid)
    print('Perimeter:', ell.perimeter)
    print('Diameter:', ell.diameter)

    plt.plot(ell.x,ell.y,'.-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Ellipse with {} panels'.format(n))
    plt.axis('equal')
    plt.grid(True)
    plt.show()
