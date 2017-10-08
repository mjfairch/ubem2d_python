import os
import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import ubem2d as ubem

if __name__ == '__main__':
    this_dir = os.path.dirname(os.path.abspath(__file__))

    # Set up and solve flow problem
    code = '0001'  # NACA 4-digit code number
    npan = 50      # number of panels
    pp = 0         # pitch position (0 = LE, .5 = midchord, 1 = LE, etc.)
    aoa_deg = 10   # angle of attack (degrees)
    uinf = (1,0)   # Onset flow
    foil = ubem.naca4(code,npan).pitch(aoa_deg*np.pi/180,pp)
    soln = ubem.solve_hess_smith_body(uinf,foil)

    # Compute aerodynamic forces and moments
    (CD,CL,CM) = ubem.airfoil_cdclcm(uinf, foil, soln.cp, pp) 
    print('Drag coefficient, CD:',CD)
    print('Lift coefficient, CL:',CL)
    print('Moment coefficient, CM:',CM)
    print('Moment axis:',foil.chord_point(pp))

    # Plot airfoil geometry
    plt.figure()
    plt.plot(foil.x/foil.chord, foil.y/foil.chord,'.k-')
    plt.axis('equal')
    plt.grid(True)
    plt.xlabel('x/c')
    plt.ylabel('y/c')
    plt.title('NACA {} with {} panels at {} degrees'.format(code,npan,aoa_deg))
    plt.savefig(os.path.join(this_dir,'naca{}_geom.pdf'.format(code)))

    # Plot pressure distribution along boundary.
    plt.figure()
    plt.plot(foil.xmid[0:foil.le+1]/foil.chord,soln.cp[0:foil.le+1],'bo-')
    plt.plot(foil.xmid[foil.le:]/foil.chord,soln.cp[foil.le:],'ro-')
    plt.axis([0,1,-10,2])
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.xlabel('x/c')
    plt.ylabel('Pressure coefficient, Cp')
    plt.title('NACA {} with {} panels at {} degrees\nCD={:.3f}, CL={:.3f}, ' \
        'CM={:.3f}, axis/chord={:.3f}'.format(code,npan,aoa_deg,CD,CL,CM,pp))
    plt.legend(['Upper surface','Lower surface'])
    plt.savefig(os.path.join(this_dir,'naca0001_cp.pdf'))

    # Plot pressure field on mesh
    (X,Y) = ubem.mesh(foil,100,100,.5)
    (Us,Vs) = ubem.velocity_source_body(foil,soln.sigma,X,Y)
    (Uv,Vv) = ubem.velocity_vortex_body(foil,soln.gamma*np.ones(foil.nedge),X,Y)
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
    plt.title('Pressure field (Cp)')
    plt.savefig(os.path.join(this_dir,'naca0001_cp_field.pdf'))

    # Stream plot of velocity field
    plt.figure()
    plt.streamplot(X[0,:],Y[:,0],Us+Uv+uinf[0],Vs+Vv+uinf[1],density=1,
        linewidth=.5)
    plt.fill(foil.x,foil.y,'k')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('NACA {} with {} panels at {} degrees'.format(code,npan,aoa_deg))
    plt.axis('equal')
    plt.savefig(os.path.join(this_dir,'naca0001_stream.pdf'))

    # Show all plots and block
    plt.show()
