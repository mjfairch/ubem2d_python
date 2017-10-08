import numpy as np
from ubem2d.util.coroutines import function_stepper
from ubem2d.math.FourierSeries import FourierSeries
from ubem2d.util.arrayify import arrayify
from ubem2d.unsteady.time_step import time_stepper

__all__ = ['fourier_pitch', 'fourier_heave', 'fourier_pitch_heave']

def fourier_pitch(time, freq, pamps, pphases):
    return fourier_pitch_heave(time, freq, pamps, pphases)

def fourier_heave(time, freq, hamps, hphases):
    return fourier_pitch_heave(time, freq, hamps=hamps, hphases=hphases)

def fourier_pitch_heave(time, freq, pamps=None, pphases=None, hamps=None,
    hphases=None):
    '''
    This returns a generator which yields time and pitch,heave data.

    freq:       Frequency in Hz
    pamps:      Pitch amplitudes (degrees)
    pphases:    Pitch phases (degrees)
    hamps:      Heave amplitudes
    hphases:    Heave phases (degrees)
    '''
    if (pamps is not None):
        pamps = arrayify(pamps)
    if (pphases is not None):
        pphases = arrayify(pphases)*np.pi/180
    if (hamps is not None):
        hamps = arrayify(hamps)
    if (hphases is not None):
        hphases = arrayify(hphases)*np.pi/180
    pitch = FourierSeries(freq, pamps, pphases)
    heave = FourierSeries(freq, hamps, hphases)
    return function_stepper(time, lambda t: (pitch(t), heave(t)))

if __name__ == '__main__':
    import numpy as np
    import matplotlib.pyplot as plt
    freq = 1.
    pamps = [5, 2.5]
    pphases = [0,-45]
    hamps = [.3, .15]
    hphases = [-90, 0]
    dt = .02/freq
    ncycles = 2

    nsteps = int(1+ncycles/(freq*dt))
    data = np.zeros((nsteps,3))
    gait = fourier_pitch_heave(time_stepper(dt, nsteps), freq, pamps, pphases,
        hamps, hphases)
    for i, (t, (alp, y)) in enumerate(gait):
        data[i,:] = t, alp, y
        print(('{:>5d} '+3*'{:>7.3f} ').format(i, t, alp, y))

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(data[:,0], data[:,1], '.-')
    plt.grid(True)
    plt.ylabel('Pitch')
    plt.title('Fourier gait')
    plt.subplot(2, 1, 2)
    plt.plot(data[:,0], data[:,2], '.-')
    plt.grid(True)
    plt.ylabel('Heave')
    plt.xlabel('Time')
    plt.show()
