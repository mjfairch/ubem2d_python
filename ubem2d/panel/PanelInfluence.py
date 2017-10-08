'''
This module provides functions to compute the stream function and velocity
velocity field at user-specified mesh points due to panels along which are
distributed singular solutions (e.g. sources/sinks, vortices).
'''
import numpy as np
from ubem2d.Errors import SizeMismatchError
from ubem2d.util.arrayify import arrayify
from ubem2d.fluids.BasicFlows import velocity_source
from ubem2d.fluids.BasicFlows import sf_source
from ubem2d.fluids.BasicFlows import sf_vortex

__all__ = ['sf_source_panel', 'sf_vortex_panel', 'velocity_source_panel',
    'velocity_vortex_panel', 'source_influence_matrices',
    'vortex_influence_matrices']

# ------------------------------------------------------------
# Special panel integrals
# ------------------------------------------------------------
def panel_integral_1(nu,mu,b,c,L):
    '''
    Evalute the panel integral of type 1, i.e. the definite integral

    \int_0^L\frac{nu*s+mu}{s^2+b*s+c}ds.

    nu,L: scalars
    mu,b,c: numpy arrays of arbitrary positive dimension, but of same shape
    Each entry of c must be nonzero.
    '''
    if (mu.shape != b.shape or mu.shape != c.shape):
        raise SizeMismatchError()
    d = 4*c-b*b
    t = 2*mu-b*nu
    z = .5*nu*np.log(np.abs((L*L+b*L+c)/c))

    i = np.where(d>0)
    sdi = np.sqrt(d[i])
    z[i] += t[i]/sdi*(np.arctan((2*L+b[i])/sdi)-np.arctan(b[i]/sdi))
    return z

# ------------------------------------------------------------
# Complex potential
# ------------------------------------------------------------

# ------------------------------------------------------------
# Stream function
# ------------------------------------------------------------
def sf_source_panel(x1,y1,tx,ty,edge,s,X,Y,m=5):
    '''
    Return the stream function at X,Y due to source sheet panels with
    initial corners (x1,y1), tangent vectors (tx,ty), and the given
    edge lengths.  Use m subintervals in the Riemann sum.
    '''
    x1,y1,tx,ty,edge,s = arrayify(x1,y1,tx,ty,edge,s)
    n = len(edge)
    if (len(x1) != n or len(y1) != n or len(tx) != n or len(ty) != n or
        len(s) != n or X.shape != Y.shape):
        raise ValueError('Size mismatch')
    Z = np.zeros(X.shape)
    dl = edge/m
    xx = x1 + (.5*dl)*tx
    yy = y1 + (.5*dl)*ty
    for i in range(m):
        Z += sf_source(s*dl,xx,yy,X,Y)
        xx += dl*tx
        yy += dl*ty
    return Z

def sf_vortex_panel(x1,y1,tx,ty,edge,s,X,Y,m=5):
    x1,y1,tx,ty,edge,s = arrayify(x1,y1,tx,ty,edge,s)
    n = len(edge)
    if (len(x1) != n or len(y1) != n or len(tx) != n or len(ty) != n or
        len(s) != n or X.shape != Y.shape):
        raise ValueError('Size mismatch')
    Z = np.zeros(X.shape)
    dl = edge/m
    xx = x1 + (.5*dl)*tx
    yy = y1 + (.5*dl)*ty
    for i in range(m):
        Z += sf_vortex(s*dl,xx,yy,X,Y)
        xx += dl*tx
        yy += dl*ty
    return Z

# ------------------------------------------------------------
# Velocity potential
# ------------------------------------------------------------

# ------------------------------------------------------------
# Velocity fields
# ------------------------------------------------------------
def velocity_source_panel(x1,y1,tx,ty,edge,s,X,Y):
    '''
    Return the net velocity U,V induced at mesh points X,Y by a source panel
    or panels encoded by x1,y1,tx,ty,edge and with strengths s.
    '''
    x1,y1,tx,ty,edge,s = arrayify(x1,y1,tx,ty,edge,s)
    n = len(edge)
    if (len(x1) != n or len(y1) != n or len(tx) != n or len(ty) != n or
        len(s) != n or X.shape != Y.shape):
        raise ValueError('Size mismatch')
    U = np.zeros(X.shape)
    V = np.zeros(Y.shape)
    for i in range(n):
        dX = X-x1[i]
        dY = Y-y1[i]
        B = -2*(dX*tx[i] + dY*ty[i])
        C = dX**2 + dY**2
        U += (.5*s[i]/np.pi)*panel_integral_1(-tx[i],dX,B,C,edge[i])
        V += (.5*s[i]/np.pi)*panel_integral_1(-ty[i],dY,B,C,edge[i])
    return (U,V)

def velocity_vortex_panel(x1,y1,tx,ty,edge,s,X,Y):
    (U,V) = velocity_source_panel(x1,y1,tx,ty,edge,s,X,Y)
    return (-V,U)

# ------------------------------------------------------------
# Influence matrices
# ------------------------------------------------------------
def source_influence_matrices(x1,y1,tx,ty,nx,ny,edge):
    n = len(edge)
    xmid = x1 + .5*(tx*edge)
    ymid = y1 + .5*(ty*edge)
    At = np.zeros((n,n))
    An = np.zeros((n,n))
    # Build influence matrices column by column
    for j in range(n):
        (u,v) = velocity_source_panel(x1[j],y1[j],tx[j],ty[j],edge[j],1.,
            xmid,ymid)
        At[:,j] = u*tx + v*ty
        An[:,j] = u*nx + v*ny
    # Update diagonal entries (based on hand computation)
    for i in range(n):
        An[i,i] = .5
        At[i,i] = 0
    return (At,An)

def vortex_influence_matrices(x1,y1,tx,ty,nx,ny,edge):
    n = len(edge)
    xmid = x1 + .5*(tx*edge)
    ymid = y1 + .5*(ty*edge)
    Bt = np.zeros((n,n))
    Bn = np.zeros((n,n))
    # Build influence matrices column by column
    for j in range(n):
        (u,v) = velocity_vortex_panel(x1[j],y1[j],tx[j],ty[j],edge[j],1.,
            xmid,ymid)
        Bt[:,j] = u*tx + v*ty
        Bn[:,j] = u*nx + v*ny
    # Update diagonal entries (based on hand computation)
    for i in range(n):
        Bt[i,i] = .5
        Bn[i,i] = 0
    return (Bt,Bn)
