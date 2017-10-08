import numpy as np
import numpy.linalg as nla
from ubem2d.geometry.Orientation import Orientation
from ubem2d.Errors import SizeMismatchError

__all__ = ['forces', 'force_coefficients', 'drag_lift_vectors',
    'drag_lift_coefficients', 'moments', 'moment_coefficient', 'body_cdclcm',
    'airfoil_cdclcm']

def forces(edge, nx, ny, cp):
    '''
    Compute the force on each panel given the pressure coefficient at the
    panel midpoint.  Closed bodies are assumed to have outward-pointing
    normal vectors, which explains the minus sign in the formulas for Fx
    and Fy below.
    '''
    n = len(edge)
    if (len(cp) != n or len(nx) != n or len(ny) != n):
        raise SizeMismatchError()
    Fx = -cp*edge*nx
    Fy = -cp*edge*ny
    return (Fx,Fy)

def force_coefficients(Fx, Fy, char_len):
    '''
    Compute the force coefficient given a force distribution and a
    characteristic length.
    '''
    n = len(Fx)
    if (len(Fy) != n):
        raise SizeMismatchError()
    CFx = Fx.sum()/char_len
    CFy = Fy.sum()/char_len
    return (CFx,CFy)  # scalar, scalar

def drag_lift_vectors(uinf):
    '''
    Return unit vectors D and L in the drag and lift directions.  The lift
    direction L is obtained by rotating the drag direction D by 90 degrees 
    CCW/CW according as to whether D points rightward/leftward. (The 
    underlying assumption is that the local gravitational acceleration is
    directed along the -y axis.)
    '''
    D = uinf/nla.norm(uinf)
    if (D[0] >= 0):
        L = (-D[1],D[0])
    else:
        L = (D[1],-D[0])
    return (D,L)

def drag_lift_coefficients(CFx, CFy, uinf):
    '''
    Compute the aerodynamic drag and lift coefficients given the force
    coefficients and a vector uinf indicating the drag direction.
    '''
    # Construct unit vectors D and L in the drag and lift directions.
    (D,L) = drag_lift_vectors(uinf)
    CD = CFx*D[0] + CFy*D[1] # == dot((CFx,CFy),D)
    CL = CFx*L[0] + CFy*L[1] # == dot((CFx,CFy),L)
    return (CD,CL)

def moments(Fx, Fy, x, y, x0 = 0., y0 = 0., sense = Orientation.CCW):
    '''
    Compute the moments about the axis through (x0,y0) given forces (Fx,Fy)
    applied at the points (x,y).  Moments are taken positive in the orientation
    specified by the last argument.
    '''
    if (sense == Orientation.CCW):
        mvec = (x-x0)*Fy - (y-y0)*Fx
    elif (sense == Orientation.CW):
        mvec = -(x-x0)*Fy + (y-y0)*Fx
    else:
        raise TypeError('Invalidation orientation type')
    return mvec

def moment_coefficient(mvec, char_len):
    '''
    Compute the moment coefficient given a vector of moments and a
    characteristic length.
    '''
    return mvec.sum()/char_len**2

def body_cdclcm(uinf, body, cp, x0 = 0., y0 = 0., moment_orientation = 
    Orientation.CCW):
    '''
    Compute the drag, lift, and and moment coefficients for a body given the 
    onset flow and the pressure distribution around the body.  Moments are
    computed about (x0,y0) and are taken positive in the sense specified by
    the last argument.  The body's diameter is taken as the characteristic
    length.
    '''
    (Fx,Fy) = forces(body.edge, body.nx, body.ny, cp)
    (CFx,CFy) = force_coefficients(Fx, Fy, body.diameter)
    (CD,CL) = drag_lift_coefficients(CFx, CFy, uinf)
    mvec = moments(Fx, Fy, body.xmid, body.ymid, x0, y0, moment_orientation)
    CM = moment_coefficient(mvec, body.diameter)
    return (CD,CL,CM)

def airfoil_cdclcm(uinf, foil, cp, pp = 0.):
    '''
    Compute the drag, lift, and and moment coefficients for an airfoil given
    the onset flow and the pressure distribution around the body.  Moments
    are computed about the point (x0,y0) on the chord line specified by pp,
    with pp=0 giving leading edge and pp=1 giving the trailing edge.  Moments
    that induce a pitch-up rotation are taken as positive.  The lift
    direction is determined as in the function body_cdclcm.  
    '''
    (x0,y0) = foil.chord_point(pp)
    return body_cdclcm(uinf, foil, cp, x0, y0, foil.pitch_up)
