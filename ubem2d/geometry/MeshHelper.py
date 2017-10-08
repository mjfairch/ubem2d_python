import numpy as np
from ubem2d.Errors import SizeMismatchError
from scipy.spatial import Delaunay

__all__ = ['boundary', 'boundary_to_string', 'mesh']

def boundary(bodies,X,Y):
    '''
    Return a list I of indicator matrices, one for each body in bodies.  Each
    I[k] has the same size as X and Y, and has values:

    0: if the corresponding mesh point is on or inside body k,
    1: if the corresponding mesh point is outside body k.

    If the first arugment is a single body (i.e. not a list of bodies), its
    indicator matrix, rather than a 1-element list of such, is returned.

    This function works for arbitrary bodies, even non-convex ones, by 
    computing the Delaunay triangulation of the body, which partitions the 
    body into triangular simplicies in a nice way.  Once the body is 
    triangulated, it is a simple matter to check if a given point lies in 
    any of the body's simplices.
    '''
    if (X.shape != Y.shape):
        raise SizeMismatchError()
    list_input = True
    if (type(bodies) not in [list,tuple]):
        list_input = False
        bodies = [bodies]

    mesh_points = np.array([[[x,y] for x in X[0,:]] for y in Y[:,0]])
    nb = len(bodies)
    I = nb*[None]
    for i in range(nb):
        xmid,ymid,nedge = bodies[i].xmid, bodies[i].ymid, bodies[i].nedge
        xy = np.array([[xmid[j],ymid[j]] for j in range(nedge)])
        J = Delaunay(xy).find_simplex(mesh_points)
        J[np.where(J>=0)] = 0
        I[i] = (-1)*J

    if (not list_input):
        I = I[0]
    return I

def boundary_to_string(I, inside='O ', outside='. '):
    '''
    Return a string representation of an indicator matrix, useful for
    visualizing which mesh nodes lie interior or exterior to a body or bodies.

    The default arguments have two characters because in most terminals the
    vertical pixel spacing between adjacent lines is about twice that between 
    adjacent characters on a given line.  So, by using two characters we are
    preserving this aspect ratio.
    '''
    c = len(inside)
    if (len(outside) != c):
        raise SizeMismatchError()
    rows = I.shape[0]
    cols = I.shape[1]
    slist = (rows*(cols+1)-1)*[None]
    k = 0
    for i in range(rows):
        for j in range(cols):
            if (I[i,j] == 0):
                slist[k] = inside
            else:
                slist[k] = outside
            k += 1
        if (i < rows-1):
            slist[k] = '\n'
        k += 1
    return ''.join(slist)

def mesh(bodies, nx=50, ny=50, margin=.5, equal=True):
    '''
    Return a uniform mesh bodies with nx,ny corners in the x,y directions and
    which encloses the given bodies.
    '''
    if (type(bodies) not in [list,tuple]):
        bodies = [bodies]
    Nb = len(bodies)
    # Determine tight bounding box
    xmin = min([min(b.x) for b in bodies])
    ymin = min([min(b.y) for b in bodies])
    xmax = max([max(b.x) for b in bodies])
    ymax = max([max(b.y) for b in bodies])
    dx = xmax-xmin
    dy = ymax-ymin
    if (equal and dx != dy):
        pad = np.abs(dx-dy)/2
        if (dx < dy):
            xmin, xmax, dx = xmin-pad, xmax+pad, dy
        else:
            ymin, ymax, dy = ymin-pad, ymax+pad, dx
    xpad = dx*margin
    ypad = dy*margin
    xmin, xmax = xmin-xpad, xmax+xpad
    ymin, ymax = ymin-ypad, ymax+ypad
    return np.meshgrid(np.linspace(xmin,xmax,nx),np.linspace(ymax,ymin,ny))
