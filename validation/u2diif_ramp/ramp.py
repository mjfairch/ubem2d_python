import os
import numpy as np
import matplotlib.pyplot as plt
import ubem2d as ubem

# ------------------------------------------------------------
# User input: body, flow, gait
# ------------------------------------------------------------
this_dir = os.path.dirname(os.path.abspath(__file__))
foil = ubem.Airfoil.load_from_file(os.path.join(this_dir,'von_mises.txt'))
uinf = (1,0)

tau = foil.chord/np.linalg.norm(uinf)  # Convective time
pp = .5          # Pitch position (0=LE, .5=midchord, 1=TE, etc.)
dalp = .1*180/np.pi        # Change in pitch angle (radians)
trise = 1.5*tau  # Rise time for the ramp motion
Tmax = 6*tau
nsteps = 501
ts = Tmax/nsteps

# ------------------------------------------------------------
# Solve for steady-state lift in final configuration
# ------------------------------------------------------------
dalp_rad = dalp*np.pi/180
ss_soln = ubem.solve_hess_smith_body((np.cos(dalp_rad), np.sin(dalp_rad)),
    foil)
ss_CD, ss_CL, ss_CM = ubem.airfoil_cdclcm(uinf, foil, ss_soln.cp, pp)
print('Steady-state drag coefficient (CD): {:<.3f}'.format(ss_CD))
print('Steady-state lift coefficient (CL): {:<.3f}'.format(ss_CL))
print('Steady-state moment coefficient (CM): {:<.3f}'.format(ss_CM))

# ------------------------------------------------------------
# Kinematics
# ------------------------------------------------------------
def pitch_ramp(time, dalp, trise):
    while (True):
        t = next(time)
        if (t <= 0):
            alp = 0.
        elif (t <= trise):
            alp = dalp*(3.-2*t/trise)*(t/trise)**2
        else:
            alp = dalp
        yield t, (alp, 0.)

tau = foil.chord/np.linalg.norm(uinf)     # Convective time
time = ubem.time_stepper(ts, nsteps)
gait = pitch_ramp(time, dalp, trise)

# ------------------------------------------------------------
# Solve/animate loop
# ------------------------------------------------------------
wake = ubem.PointVortexWake(eps=1.e-6)
t, alp, y = np.zeros((3, nsteps))
CT, CL, CM, Ein, Eout, bcirc = np.zeros((6, nsteps))
frames = ubem.airfoil_stepper(foil, gait, pp, wake, uinf)
for i, (kin, out) in enumerate(frames):
    print(('{:>5d} '+6*'{:>10.3g} ').format(i, *kin, *out[0:3]))
    t[i], alp[i], y[i] = kin
    CT[i], CL[i], CM[i], Ein[i], Eout[i], bcirc[i] = out

# ------------------------------------------------------------
# Post processing
# ------------------------------------------------------------
plotdir = this_dir
ttau = t/tau
ubem.plot_wake(foil, wake, save_path=os.path.join(plotdir, 'ramp_wake.pdf'))
ubem.plot_circulation(ttau, 'Time (t/tau)', wake, bcirc, os.path.join(
    plotdir, 'ramp_circulation.pdf'))
ubem.plot_kinematics(ttau, 'Time (t/tau)', alp, y, os.path.join(plotdir,
    'ramp_kinematics.pdf'))

u2diif_CL = ubem.read_data(os.path.join(this_dir, 'u2diif_CL.txt'))
plt.figure()
plt.plot(ttau, CL/ss_CL, '-')
plt.plot(u2diif_CL[:,0], u2diif_CL[:,1], 'k.')
plt.xlabel('Time (t/tau)')
plt.title('Lift coefficient (CL), normalized by steady-state value')
plt.legend(['CL', 'CL (U2DIIF)'])
plt.grid(True)
plt.savefig(os.path.join(this_dir, 'ramp_cl.pdf'))

plt.show()
