import os
import math
import numpy as np
import ubem2d as ubem
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# User input: body, flow, kinematics
# ------------------------------------------------------------
code = '0012'    # NACA 4-digit code
uinf = (1,0)     # Onset flow
npan = 50        # Number of panels
spc = 50         # Steps per cycle
pp = 0           # Pitch position (0=LE, .5=midchord, 1=TE, etc.)
k = 1            # Reduced frequency
pamps = 2        # Pitch amplitude (degrees)
pphases = 0      # Pitch phase (degrees)
hamps = 0        # Heave amplitude (multiple of chord)
hphases = 0      # Heave phase (degrees)
Tmax = 4         # Fraction of oscillation period

movie = False    # If true, show movie.  If false, show plots.
animation_path = os.path.join(ubem.__movie_dir, 'unsteady.mp4')

# ------------------------------------------------------------
# Motion
# ------------------------------------------------------------
foil = ubem.naca4(code, npan)      # Build airfoil
spdinf = np.linalg.norm(uinf)      # Speed of onset flow
tau = foil.chord/spdinf            # Convective time
freq_hz = k/tau                    # Frequency in Hz of oscillation
Tosc = 1./freq_hz                  # Period of oscillation
dt = 1./(spc*freq_hz)              # Time step
Tfinal = Tmax*Tosc                 # Max simulation time
ncycles = math.floor(Tfinal/Tosc)  # Total number of cycles in simulation
nsteps = 1 + math.ceil(Tfinal/dt)  # Total number of steps in simulation
motion = ubem.fourier_pitch_heave(ubem.time_stepper(dt, nsteps), freq_hz,
    pamps, pphases, hamps, hphases)

# ------------------------------------------------------------
# Solve/animate loop
# ------------------------------------------------------------
wake = ubem.PointVortexWake(eps=1.e-6)
t, alp, y, y_te = np.zeros((4,nsteps))
CT, CL, CM, Ein, Eout, bcirc = np.zeros((6,nsteps))
frames = ubem.airfoil_stepper(foil, motion, pp, wake, uinf)
def process_step(i, kin, out):
    print(('{:>5d} '+6*'{:>10.3g} ').format(i, *kin, *out[0:3]))
    t[i], alp[i], y[i], y_te[i] = *kin, foil.y[0]
    CT[i], CL[i], CM[i], Ein[i], Eout[i], bcirc[i] = out

if (not movie):
    for i, (kin, out) in enumerate(frames):
        process_step(i, kin, out)
else:
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(-.5, 3.5)
    ax.set_ylim(-2, 2)
    ax.set_axis_off()
    foil_line, = ax.plot([], [], 'k-')
    scatter = ax.scatter([], [], c=[], cmap='bwr', vmin=-dt, vmax=dt,
        edgecolors='none')
    def draw(frame):
        i, (kin, out) = frame[0], frame[1]
        process_step(i, kin, out)
        foil_line.set_data(foil.x, foil.y)
        scatter.set_offsets(np.array([wake.x, wake.y]).T)
        scatter.set_array(wake.gam)
        return foil_line, scatter
    ani = animation.FuncAnimation(fig, draw, enumerate(frames),
        interval=round(1000*dt), save_count=nsteps, repeat=False, blit=True)
    if (animation_path is None):
        plt.show()
    else:
        writer = animation.writers['ffmpeg'](fps=freq_hz*spc, bitrate=57600)
        ani.save(animation_path, writer=writer)

# ------------------------------------------------------------
# Post processing
# ------------------------------------------------------------
i0, i1 = (ncycles-1)*spc, ncycles*spc + 1  # Range indices of final cycle
CTF, CLF, CMF = np.mean(CT[i0:i1]), np.mean(CL[i0:i1]), np.mean(CM[i0:i1])
EinF, EoutF = np.sum(Ein[i0:i1]), np.sum(Eout[i0:i1])
effF = EoutF/EinF
St = freq_hz*(np.max(y_te[i0:i1])-np.min(y_te[i0:i1]))/np.linalg.norm(uinf)
print('\n'+'k: {}'.rjust(20).format(k))
print('steps/cycle: {}'.rjust(20).format(spc))
print('Tfinal/Tosc: {}'.rjust(20).format(Tmax))
print('dt: {}'.rjust(20).format(dt))
print('CT: {}'.rjust(20).format(CTF))
print('CL: {}'.rjust(20).format(CLF))
print('CM: {}'.rjust(20).format(CMF))
print('Strouhal: {}'.rjust(20).format(St))
print('Efficiency: {}'.rjust(20).format(effF))
if (not movie):
    tf = t*freq_hz  # For plotting, normalize time by the period of oscillation
    plotdir = ubem.__plot_dir
    ubem.plot_wake(foil, wake, freq_hz, St, CTF, CLF, CMF, effF, save_path=
        os.path.join(plotdir, 'unsteady_wake.pdf'))
    ubem.plot_circulation(tf, 'Time (t/T)', wake, bcirc, os.path.join(plotdir,
        'unsteady_circulation.pdf'))
    ubem.plot_kinematics(tf, 'Time (t/T)', alp, y, os.path.join(plotdir,
        'unsteady_kinematics.pdf'))
    ubem.plot_aerodynamics(tf*Tosc, 'Time (t)', CT, CL, CM, os.path.join(plotdir,
        'unsteady_aerodynamics.pdf'), i0=2)  # i0=2: omit impulsive loading
    plt.show()
