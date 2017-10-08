import math
import ubem2d as ubem

if __name__ == '__main__':
    # ----------------------------------------------------------------------
    # User input
    # ----------------------------------------------------------------------
    save_path = None  # If not none, save (rather than display) the animation
    bitrate = 57600   # Animation bitrate when saving to disk

    foil = ubem.naca4('0010', 50)
    pp = 0
    freq = 1
    dt = .01/freq
    pamps = 5
    pphases = -90
    hamps = .2
    hphases = 0
    ncycles = 3
    nsteps = 1 + math.ceil(ncycles/(freq*dt))
    time = ubem.time_stepper(dt, nsteps)
    motion = ubem.fourier_pitch_heave(time, freq, pamps, pphases, hamps,
        hphases)

    # ----------------------------------------------------------------------
    # Animate!
    # ----------------------------------------------------------------------
    fps = 1./(freq*dt)  # so that animation rate == kinematic rate
    if (save_path is not None):
        print('Writing animation to file:', save_path)
    ubem.animate_motion(foil, pp, motion, fps, save_path, nsteps, bitrate)
