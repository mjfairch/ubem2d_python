import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import ubem2d as ubem

alp = 0.  # Inclination (radians) of onset flow to the x axis
uinf = [np.cos(alp),np.sin(alp)]
xx = np.linspace(-1,1,50)
yy = np.linspace(-1,1,50)
(X,Y) = np.meshgrid(xx,yy)

def plot_and_save(title, fname, *args):
    plt.figure()
    matplotlib.rcParams['contour.negative_linestyle']='solid'
    for arg in args:
        plt.contour(X, Y, arg, 50)#,colors='k')
    plt.axis('equal')
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(title)
    plt.savefig(os.path.join(ubem.__plot_dir, fname))

# Stream function due to a uniform flow
Z = ubem.sf_uniform_flow(uinf,X,Y)
title = 'Stream lines for the uniform flow ({:.3f},{:.3f})'.format(
    uinf[0],uinf[1])
plot_and_save(title,'sf_uniform_flow.pdf',Z)

# Stream function due to a source
Z = ubem.sf_source(1.,0.,0.,X,Y)
title = 'Stream lines for a source at the origin'
plot_and_save(title,'sf_source.pdf',Z)

# Stream function due to a source
Z = ubem.sf_vortex(1.,0.,0.,X,Y)
title = 'Stream lines for a vortex at the origin'
plot_and_save(title,'sf_vortex.pdf',Z)

# Stream function due to a doublet
Z = ubem.sf_doublet(1.,0.,0.,0.,X,Y)
title = 'Stream lines for a doublet at the origin'
plot_and_save(title,'sf_doublet.pdf',Z)

# Stream function for a doublet and a uniform flow
Z = ubem.sf_doublet(-1,0,0,0,X,Y) + ubem.sf_uniform_flow(uinf,X,Y)
title = 'Stream lines for a doublet at the origin and a uniform flow'
plot_and_save(title,'sf_doublet_uniform_flow.pdf',Z)

plt.show()
