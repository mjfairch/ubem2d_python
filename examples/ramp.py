import os
import numpy as np
import matplotlib.pyplot as plt
import ubem2d as ubem

if __name__ == '__main__':
    x = np.linspace(-1.,3,100)
    y0 = ubem.exponential_ramp(x)
    y1 = .1 + ubem.finite_ramp(x)
    y2 = ubem.smooth_ramp(x, -.25, .5)
    plt.figure()
    plt.plot(x, y0, x, y1, x, y2)
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(['Exponential ramp', 'Finite ramp', 'Arbitrary ramp'])
    plt.savefig(os.path.join(ubem.__plot_dir, 'ramps.pdf'))
    plt.show()
