import os
import numpy as np
import matplotlib.pyplot as plt
import ubem2d as ubem

x = np.array([0., 2, 2, 1, 1, 0, 0])
y = np.array([0., 0, 2, 2, 1, 1, 0])
plt.figure()
plt.plot(x, y, '.-')
plt.axis('equal')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Turning angle = {:.2f} rad'.format(ubem.turning_angle(x, y, True)))
plt.savefig(os.path.join(ubem.__plot_dir, 'turning_angle.pdf'))
plt.show()
