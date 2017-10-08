import numpy as np
from ubem2d.util.arrayify import arrayify
from ubem2d.Errors import SizeMismatchError

__all__ = ['FourierSeries']

class FourierSeries():
    '''
    This class represents the following Fourier series

    y = \sum_{n=1}^N amplitudes[i]*sin(2*pi*n*f1*t - phases[n]),

    where f1 is the fundamental frequency, and the list of amplitudes and
    phases must have the same shape.
    '''
    def __init__(self, base_frequency, amplitudes = 0., phases = None):
        '''
        base_frequency: nonzero scalar frequency of fundamental mode (Hz)
        amplitudes: scalar, list, tuple, or 1-dimensional numpy array of mode
            amplitudes.
        phases: scalar, list, tuple, or 1-dimensional numpy array of mode 
            phases (radians) and of same shape as amplitudes.  If None, it is
            initialized to a zero array of the same shape as amplitudes.
        '''
        if (not np.isscalar(base_frequency)):
            raise ValueError('Base frequency must be a scalar')
        if (base_frequency == 0):
            raise ValueError('Base frequency cannot be zero')
        if (amplitudes is None):
            amplitudes = np.array([0.])
        else:
            amplitudes = arrayify(amplitudes)
        if (phases is not None):
            phases = arrayify(phases)
        else:
            phases = np.zeros(amplitudes.shape)
        if (amplitudes.ndim != 1):
            raise ValueError('Bad shape for mode amplitudes')
        if (phases.shape != amplitudes.shape):
            raise SizeMismatchError('Phase and amplitude shapes must match')
        self._base_frequency = base_frequency
        self._amplitudes = amplitudes
        self._phases = phases

    def __len__(self):
        '''
        Return the number of modes (even counting those with zero amplitude) 
        in this Fourier series
        '''
        return len(self._amplitudes)
    
    @property
    def base_frequency(self):
        return self._base_frequency
    
    @property
    def amplitudes(self):
        return self._amplitudes

    @property
    def phases(self):
        return self._phases
    
    @property
    def period(self):
        return 1./self._base_frequency
    
    def __call__(self, t, modes = None):
        '''
        Sum modes of this Fourier series at the scalar, list, tuple, or numpy 
        ndarray t.  If modes is None, all modes are evaluated.  Otherwise,
        modes is expected to be a list of mode indices to include in the sum.
        '''
        if (np.isscalar(t)):
            y = 0.
        else:
            t = arrayify(t)
            y = np.zeros(t.shape)
        if (not modes):
            modes = range(len(self))
        for n in modes:
            if (not n==int(n) or n < 0 or n >= len(self)):
                raise ValueError('Bad mode index: {}'.format(n))
            y += self._amplitudes[n]*np.sin(
                2*np.pi*(n+1)*self._base_frequency*t - self._phases[n])
        return y
