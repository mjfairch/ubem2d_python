import numpy as np
from .Body import Body
from .Orientation import Orientation
from .Cylinder import Cylinder

__all__ = ['CircularCylinder']

class CircularCylinder(Cylinder):
    def __init__(self, npanels = 50, radius = 1, x0 = 0, y0 = 0):
        npanels = int(npanels)
        self._radius = float(radius)
        if (self._radius <= 0):
            raise ValueError('Circle must have positive radius')
        super().__init__(lambda th: self._radius, npanels, x0, y0)
    
    @property
    def radius(self):
        return self._radius

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    n = int(input('Enter number of panels: '))
    cyl = CircularCylinder(n)
    for i in range(n+1):
        print('{}\t{}'.format(cyl.x[i],cyl.y[i]))
    print('\nRadius:', cyl.radius)
    print('Center:', cyl.center)
    print('Closed:', cyl.closed)
    print('Centroid:', cyl.centroid)
    print('Perimeter:', cyl.perimeter)
    print('Diameter:', cyl.diameter)

    plt.plot(cyl.x,cyl.y,'.-')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Circular cylinder with {} panels'.format(n))
    plt.axis('equal')
    plt.grid(True)
    plt.show()
