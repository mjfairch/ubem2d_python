''' 
Four basic solutions to Laplace's equation for potential flow are:

1. Uniform flow
2. Flow due to a point source/sink
3. Flow due to a doublet (a limit of a source/sink pair)
4. Flow due to a point vortex

This module provides functions to compute the complex potential, velocity
potential, stream function, and velocity field at user-specified mesh points
due a uniform flow or a collection of sources/sinks, doublets, or point 
vortices.

Throughout this module:
uinf:   tuple defining the onset flow
s:      array of strengths
x,y:    coordinates of sources, vortices, doublets, etc.
X,Y:    coordinates of field points at which to evaluate influence
alpha:  for doublets, array defining axis (as an angle) of each doublet
'''

import numpy as np
import numpy.linalg as nla
from ubem2d.util.arrayify import arrayify
from ubem2d.Errors import SizeMismatchError

__all__ = ['cp_uniform_flow', 'cp_source', 'cp_doublet', 'cp_vortex',
    'sf_uniform_flow', 'sf_source', 'sf_doublet', 'sf_vortex',
    'vp_uniform_flow', 'vp_source', 'vp_doublet', 'vp_vortex',
    'velocity_uniform_flow', 'velocity_source', 'velocity_doublet',
    'velocity_vortex']

def setup(s,x,y,X,Y,alpha=None):
    '''
    This function has logic factored out of the basic flow functions below.
    
    s: strengths of sources, vortices, doublets, etc.
    (x,y): coordinates of sources, vortices, doublets, etc.
    (X,Y): a mesh over which flow quantities are to be computed
    
    Convert s,x,y (and alpha, if provided) from scalar, tuple, or list type 
    into numpy arrays.  Ensure that s,x,y  (and alpha, if provided) have the
    equal shape.  Ensure that X,Y have equal shape.
    
    Return updated values of s,x,y (and alpha, if provided).
    '''
    s,x,y = arrayify(s,x,y)
    if (s.shape != x.shape or x.shape != y.shape or X.shape != Y.shape):
        raise SizeMismatchError()
    if (alpha is not None):
        alpha = arrayify(alpha)
        if (alpha.shape != s.shape):
            raise SizeMismatchError()
        return (s,x,y,len(s),alpha)
    return (s,x,y,len(s))
        
# ------------------------------------------------------------
# Complex potential
# ------------------------------------------------------------
def cp_uniform_flow(uinf,X,Y):
    if (X.shape != Y.shape):
        raise SizeMismatchError()
    Z = X + 1j*Y
    U = nla.norm(uinf)
    alpha = np.arctan2(uinf[1],uinf[0])
    return (U*np.exp(-1j*alpha))*Z

def cp_source(s,x,y,X,Y):
    s,x,y,n = setup(s,x,y,X,Y)
    Z = X + 1j*Y
    W = (0.+0j)*np.zeros(X.shape)
    for i in range(n):
        W += s[i]*np.log(Z-(x[i]+1j*y[i]))
    return (.5/np.pi)*W

def cp_doublet(s,x,y,alpha,X,Y):
    s,x,y,n,alpha = setup(s,x,y,X,Y,alpha)
    Z = X + 1j*Y
    W = (0.+0j)*np.zeros(X.shape)
    for i in range(n):
        W += (s[i]*np.exp(1j*alpha[i]))/(Z-(x[i]+1j*y[i]))
    return (-.5/np.pi)*W

def cp_vortex(s,x,y,X,Y):
    return (-1.j)*cp_source(s,x,y,X,Y)

# ------------------------------------------------------------
# Stream function
# ------------------------------------------------------------
def sf_uniform_flow(uinf,X,Y):
    if (X.shape != Y.shape):
        raise SizeMismatchError()
    return uinf[0]*Y - uinf[1]*X

def sf_source(s,x,y,X,Y):
    s,x,y,n = setup(s,x,y,X,Y)
    Z = np.zeros(X.shape)
    for i in range(n):
        Z += s[i]*np.arctan2(Y-y[i],X-x[i])
    return (.5/np.pi)*Z

def sf_doublet(s,x,y,alpha,X,Y):
    s,x,y,n,alpha = setup(s,x,y,X,Y,alpha)
    Z = np.zeros(X.shape)
    for i in range(n):
        dX = X-x[i]
        dY = Y-y[i]
        Z += s[i]*(np.cos(alpha[i])*dY - np.sin(alpha[i])*dX)/(dX**2 + dY**2)
    return (.5/np.pi)*Z

def sf_vortex(s,x,y,X,Y):
    s,x,y,n = setup(s,x,y,X,Y)
    Z = np.zeros(X.shape)
    for i in range(n):
        dX = X-x[i]
        dY = Y-y[i]
        Z += s[i]*np.log(dX**2 + dY**2)
    return (-.25/np.pi)*Z

# ------------------------------------------------------------
# Velocity potential
# ------------------------------------------------------------
def vp_uniform_flow(uinf,X,Y):
    if (X.shape != Y.shape):
        raise SizeMismatchError()
    return uinf[0]*X + uinf[1]*Y

def vp_source(s,x,y,X,Y):
    return -sf_vortex(s,x,y,X,Y)

def vp_doublet(s,x,y,alpha,X,Y):
    s,x,y,n,alpha = setup(s,x,y,X,Y,alpha)
    Z = np.zeros(X.shape)
    for i in range(n):
        dX = X-x[i]
        dY = Y-y[i]
        Z += s[i]*(np.cos(alp)*dX + np.sin(alp)*dY)/(dX**2 + dY**2)
    return (-.5/np.pi)*Z

def vp_vortex(s,x,y,X,Y):
    return sf_source(s,x,y,X,Y)

# ------------------------------------------------------------
# Velocity fields
# ------------------------------------------------------------
def velocity_uniform_flow(uinf,X,Y):
    if (X.shape != Y.shape):
        raise SizeMismatchError()
    return (uinf[0]*np.ones(X.shape), uinf[1]*np.ones(Y.shape))

def velocity_source(s,x,y,X,Y):
    s,x,y,n = setup(s,x,y,X,Y)
    Vx,Vy = np.zeros(X.shape), np.zeros(Y.shape)
    for i in range(n):
        dX = X-x[i]
        dY = Y-y[i]
        distance = dX**2 + dY**2
        Vx += s[i]*dX/distance
        Vy += s[i]*dY/distance
    return ((.5/np.pi)*Vx, (.5/np.pi)*Vy)

def velocity_doublet(s,x,y,alpha,X,Y):
    s,x,y,n,alpha = setup(s,x,y,X,Y,alpha)
    Vx,Vy = np.zeros(X.shape), np.zeros(Y.shape)
    for i in range(n):
        dX = X-x[i]
        dY = Y-y[i]
        cosalp = np.cos(alpha)
        sinalp = np.sin(alpha)
        diff_term = dX**2 - dY**2
        product_term = 2*dX*dY
        distance = dX**2 + dY**2
        Vx += s[i]*(diff_term*cosalp + product_term*sinalp)/distance
        Vy += s[i]*(product_term*cosalp - diff_term*sinalp)/distance
    return ((.5/np.pi)*Vx, (.5/np.pi)*Vy)

def velocity_vortex(s,x,y,X,Y):
    (Vx,Vy) = velocity_source(s,x,y,X,Y)
    return (-Vy,Vx)
