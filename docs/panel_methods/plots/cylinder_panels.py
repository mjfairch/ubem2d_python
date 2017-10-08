import numpy as np
import matplotlib.pyplot as plt
import os

thisdir = os.path.dirname(__file__)

npan = 8
t = np.linspace(0,2*np.pi,npan+1)
x,y = np.cos(t), np.sin(t)
tx,ty = x[1:]-x[:-1], y[1:]-y[:-1]
edge = np.sqrt(tx*tx+ty*ty)
tx, ty = tx/edge, ty/edge
nx, ny = (-ty, tx)
xmid,ymid = .5*(x[:-1]+x[1:]), .5*(y[:-1]+y[1:])

tt = np.linspace(0,2*np.pi,100)
cyl_x, cyl_y = np.cos(tt), np.sin(tt)

plt.plot(x,y,'k-s')
plt.plot(xmid,ymid,'ro')
plt.plot(cyl_x,cyl_y,'b-')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.axis('equal')
plt.grid(True)
plt.savefig(os.path.join(thisdir,'cylinder_panels.pdf'))
plt.show()
