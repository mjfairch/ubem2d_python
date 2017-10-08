import os
import numpy as np
import matplotlib.pyplot as plt
import ubem2d as ubem

# User input: bodies and mesh
cylinder = ubem.CircularCylinder(50,.5,0,0)
ellipse = ubem.Ellipse(50,.5,.25,-1,1)
rectangle = ubem.Rectangle(.5,.25,10,10,1,-.5)
L_shape = ubem.BrokenLine([0.,1,1,.5,.5,0,0],[1.,1,1.5,1.5,1.25,1.25,1])
bodies = [cylinder,ellipse,rectangle,L_shape]
(X,Y) = ubem.mesh(bodies,50,50,0,True)

# Determine boundary of each body
I = ubem.boundary(bodies,X,Y)

# Plot bodies, mesh points, and interior points
plt.figure()
for body in bodies:
    plt.plot(body.x, body.y)
xmesh = [x for x in X[0,:] for y in Y[:,0]]
ymesh = [y for x in X[0,:] for y in Y[:,0]]
plt.plot(xmesh,ymesh,'k.',markersize=.5)
for ij in I:
    (i,j) = np.where(ij==0)  # (i,j) indices to grid points interior to body
    n = len(i)
    x = [X[i[k],j[k]] for k in range(n)]
    y = [Y[i[k],j[k]] for k in range(n)]
    plt.plot(x, y, 'r.', markersize=2)
plt.axis('equal')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Object boundaries in a mesh via Delaunay triangulation')
plt.savefig(os.path.join(ubem.__plot_dir, 'boundary.pdf'))
plt.show()
