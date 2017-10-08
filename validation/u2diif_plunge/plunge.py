import os
import numpy as np
import matplotlib.pyplot as plt
import ubem2d as ubem

this_dir = os.path.dirname(os.path.abspath(__file__))
foil = ubem.naca4('0015', 100)

uinf = (1,0)
tau = foil.chord/np.linalg.norm(uinf)
k = 2.15/np.pi        # Reduced pitch frequency (k=1: 1osc. per tau)
freq = k/tau          # Frequency (Hz)
pp = 0                # Pitch/moment axis (0=LE, .5=midchord, 1=TE, etc.)
hamp = -.018          # Heave amplitude (multiple of chord)
res = 50              # Resolution: # of steps at fastest time scale
T = 1./freq           # Fundamental oscillation period (sec)
(ts, spc) = ubem.time_step(res, tau, T, T)
ncycles = 4
nsteps = 1 + ncycles*spc
time = ubem.time_stepper(ts, nsteps)
gait = ubem.fourier_pitch_heave(time, freq, hamps=hamp)

# ------------------------------------------------------------
# Solve/animate loop
# ------------------------------------------------------------
wake = ubem.PointVortexWake(eps=1.e-6)
t, alp, y, y_te = np.zeros((4, nsteps))
CT, CL, CM, Ein, Eout, bcirc = np.zeros((6, nsteps))
frames = ubem.airfoil_stepper(foil, gait, pp, wake, uinf)
for i, (kin, out) in enumerate(frames):
    print(('{:>5d} '+6*'{:>10.3g} ').format(i, *kin, *out[0:3]))
    t[i], alp[i], y[i], y_te[i] = *kin, foil.y[0]
    CT[i], CL[i], CM[i], Ein[i], Eout[i], bcirc[i]  = out

# ------------------------------------------------------------
# Post processing
# ------------------------------------------------------------
u2diif_50CD = ubem.read_data(os.path.join(this_dir, 'u2diif_50CD.txt'))
u2diif_CL = ubem.read_data(os.path.join(this_dir, 'u2diif_CL.txt'))
u2diif_CM = ubem.read_data(os.path.join(this_dir, 'u2diif_CM.txt'))
tT = t/T

ubem.plot_wake(foil, wake, k, save_path=os.path.join(this_dir,
    'plunge_wake.pdf'))
ubem.plot_circulation(tT, 'Time (t/T)', wake, bcirc, os.path.join(this_dir,
    'plunge_circulation.pdf'))
ubem.plot_kinematics(tT, 'Time (t/T)', alp, y, os.path.join(this_dir,
    'plunge_kinematics.pdf'))

i0 = np.min(np.where(tT >= 2.))
plt.figure()
plt.plot(tT[i0:], -50*CT[i0:], 'g-')  # CD=-CT
plt.plot(u2diif_50CD[:,0], u2diif_50CD[:,1],'g.')
plt.plot(tT[i0:], CL[i0:], 'b-')
plt.plot(u2diif_CL[:,0], u2diif_CL[:,1], 'b.')
plt.plot(tT[i0:], CM[i0:], 'r-')
plt.plot(u2diif_CM[:,0], u2diif_CM[:,1], 'r.')
plt.xlabel('Time (t/T)')
plt.legend(['50 CD', '50 CD (U2DIIF)', 'CL', 'CL (U2DIIF)', 'CM',
    'CM (U2DIIF)'])
plt.title('Aerodynamic coefficients')
plt.grid(True)
plt.savefig(os.path.join(this_dir, 'plunge_aerodynamics.pdf'))

# Show all plots and block
plt.show()
