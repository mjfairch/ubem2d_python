import os
import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import ubem2d as ubem

if __name__ == '__main__':
    # Set up and solve flow problem
    code = '0015'  # NACA 4-digit code number
    npan = 50      # number of panels
    pp = .25       # pitch position (0 = LE, .5 = midchord, 1 = LE, etc.)
    aoa_deg = 5    # angle of attack (degrees)
    uinf = (1,0)   # Onset flow
    foil = ubem.naca4(code,npan).pitch(aoa_deg*np.pi/180,pp)
    sys = ubem.HessSmithSystem(foil)

    # Compute aerodynamic forces and moments
    soln = sys.solve(uinf)
    cp = sys.pressure_self(uinf, soln)[0]
    CD,CL,CM = ubem.airfoil_cdclcm(uinf, foil, cp, pp)
    print('Drag coefficient, CD:', CD)
    print('Lift coefficient, CL:', CL)
    print('Moment coefficient, CM:', CM)
    print('Moment axis:', foil.chord_point(pp))

    # Plot airfoil geometry
    plt.figure()
    plt.plot(foil.x/foil.chord, foil.y/foil.chord,'.-')
    plt.axis('equal')
    plt.grid(True)
    plt.xlabel('x/c')
    plt.ylabel('y/c')
    plt.title('NACA {} at {} degrees'.format(code, aoa_deg))
    plt.savefig(os.path.join(ubem.__plot_dir, 'naca_geometry.pdf'.format(code)))

    # Plot pressure distribution along boundary.
    plt.figure()
    plt.plot(foil.xmid[0:foil.le+1]/foil.chord, cp[0:foil.le+1], 'bo-')
    plt.plot(foil.xmid[foil.le:]/foil.chord, cp[foil.le:], 'ro-')
    plt.xlabel('x/c')
    plt.ylabel('Pressure coefficient (Cp)')
    plt.gca().invert_yaxis()
    plt.title('NACA {} at {} degrees\nCD={:.3f}, CL={:.3f}, CM={:.3f}, ' \
        'axis/chord={:.3f}'.format(code, aoa_deg, CD, CL, CM, pp))
    plt.legend(['Upper surface','Lower surface'])
    plt.grid(True)
    plt.savefig(os.path.join(ubem.__plot_dir, 'naca_cp_boundary.pdf'))

    # Plot pressure field on mesh
    (X,Y) = ubem.mesh(foil, 100, 100, .5)
    U,V = sys.flow_external(uinf, soln, X, Y)
    Cp = sys.pressure_from_flow(uinf, soln, U, V)
    plt.figure()
    plt.contourf(X, Y, Cp, 100)
    plt.fill(foil.x, foil.y)
    plt.colorbar()
    plt.axis('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Pressure field (Cp), NACA {} at {} degrees'.format(code,
        aoa_deg))
    plt.savefig(os.path.join(ubem.__plot_dir,'naca_cp_field.pdf'))

    # Stream plot of velocity field
    plt.figure()
    plt.streamplot(X, Y, U, V, density=1, linewidth=.5, color='k')
    plt.fill(foil.x, foil.y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('NACA {} at {} degrees'.format(code, aoa_deg))
    plt.axis('equal')
    plt.savefig(os.path.join(ubem.__plot_dir, 'naca_streamlines.pdf'))

    # Show all plots and block
    plt.show()
