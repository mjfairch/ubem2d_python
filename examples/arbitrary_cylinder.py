import os
import numpy as np
import matplotlib.pyplot as plt
import ubem2d as ubem

def plot_and_save(x,y,title,fname):
    plt.figure()
    plt.plot(x,y)
    plt.axis('equal')
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.savefig(os.path.join(ubem.__plot_dir, fname))

if __name__ == '__main__':
    n = 100
    cyl = ubem.Cylinder(lambda th: np.sin(th)-1., n)
    plot_and_save(cyl.x, cyl.y, 'Cardioid', 'cyl_cardioid.pdf')

    n = 100
    cyl = ubem.Cylinder(lambda th: 3*np.sin(th)+2., n)
    plot_and_save(cyl.x, cyl.y, 'Limacon', 'cyl_limacon.pdf')

    plt.show()
