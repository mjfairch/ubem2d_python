import os
import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import ubem2d as ubem

if __name__ == '__main__':
    '''
    Validate the pressure coefficient around a NACA 2412 airfoil at 8 degrees
    angle of attack vs. that published in Kuethe & Chow's book, the latter
    using a linear-strength vortex panel method.
    '''
    # Set up and solve flow problem
    code = '2412'  # NACA 4-digit code number
    npan = 120     # number of panels
    pp = 0         # pitch position (0 = LE, .5 = midchord, 1 = LE, etc.)
    aoa_deg = 8    # angle of attack (degrees)
    uinf = (1,0)   # Onset flow
    foil = ubem.naca4(code,npan).pitch(aoa_deg*np.pi/180,pp)
    soln = ubem.solve_hess_smith_body(uinf,foil)

    # Compute aerodynamic forces and moments
    (CD,CL,CM) = ubem.airfoil_cdclcm(uinf, foil, soln.cp, pp)
    print('Drag coefficient, CD:',CD)
    print('Lift coefficient, CL:',CL)
    print('Moment coefficient, CM:',CM)
    print('Moment axis:',foil.chord_point(pp))

    # Plot pressure distribution along boundary.
    this_dir = os.path.dirname(os.path.abspath(__file__))
    kc_CL = ubem.read_data(os.path.join(this_dir, 'kc_2412_8deg_CL.txt'))
    plt.figure()
    plt.plot(foil.xmid, soln.cp,'k-')
    plt.plot(kc_CL[:,0], kc_CL[:,1], 'bs', markerfacecolor='none')
    plt.xlabel('x/c')
    plt.ylabel('Pressure coefficient (Cp)')
    plt.gca().invert_yaxis()
    plt.title('NACA {} with {} panels at {} degrees\nCD={:.3f}, CL={:.3f}, ' \
        'CM={:.3f}, axis/chord={:.3f}'.format(code,npan,aoa_deg,CD,CL,CM,pp))
    plt.legend(['UBEM2D', 'Kuethe-Chow (digitized)'])
    plt.grid(True)
    plt.savefig(os.path.join(this_dir,'cp_boundary_2412_8deg.pdf'))

    # Plot pressure field on mesh
    (X,Y) = ubem.mesh(foil,100,100,.5)
    (Us,Vs) = ubem.velocity_source_body(foil, soln.sigma, X, Y)
    (Uv,Vv) = ubem.velocity_vortex_body(foil, soln.gamma*np.ones(foil.nedge),
        X, Y)
    U = Us+Uv+uinf[0]
    V = Vs+Vv+uinf[1]
    Z = 1. - (U**2 + V**2)/nla.norm(uinf)**2
    plt.figure()
    plt.contourf(X,Y,Z,100)
    plt.fill(foil.x,foil.y)
    plt.colorbar()
    plt.axis('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Pressure field (Cp), NACA 2412 at 8 degrees')
    plt.savefig(os.path.join(this_dir,'cp_field_2412_8deg.pdf'))

    # Show all plots and block
    plt.show()
