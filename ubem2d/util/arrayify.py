import numpy as np

__all__ = ['arrayify']

def arrayify(*args):
    '''
    Convert each of a variable number of scalars, tuples, lists, or numpy
    arrays to numpy arrays, returning a list of numpy arrays of the same
    length as the number of arguments.
    '''
    n = len(args)
    ret = n*[None]
    for i in range(n):
        argi = args[i]
        if (np.isscalar(argi)):
            ret[i] = np.array([argi])
        elif (type(argi) in [list,tuple]):
            ret[i] = np.array(argi)
        elif (type(argi) is np.ndarray):
            ret[i] = argi
        else:
            raise ValueError('Argument must be scalar, list, tuple, or ' \
                'numpy array')
    if (n == 1):
        return ret[0]
    else:
        return ret
