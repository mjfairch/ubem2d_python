import numpy as np

__all__ = ['turning_angle']

def turning_angle(x, y, closed = False):
    '''
    Compute the turning angle of the broken line whose corners are at the
    points (x,y).  Intuitively, the turning angle is the sum of the changes
    in heading angle as one moves along the broken line.

    closed: If True, add a final turn at the end to return to the initial
            heading.
    '''
    # Compute heading angle of each edge
    th = np.arctan2(np.diff(y),np.diff(x))
    if (closed):
        th = np.concatenate([th,[th[0]]])
    # Difference in heading angle between adjacent edges
    dth = np.diff(th)
    # Deal with heading changes that cross the branch cut of arctan2
    bc_ccw = np.where(dth < -np.pi)
    bc_cw = np.where(dth > np.pi)
    dth[bc_ccw] = 2*np.pi + dth[bc_ccw]
    dth[bc_cw] = dth[bc_cw] - 2*np.pi
    # The sum of the (adjusted) changes in heading angle is the turning angle
    return np.sum(dth)
