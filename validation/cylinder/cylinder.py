import os
import numpy as np
import numpy.linalg as nla
import matplotlib
import matplotlib.pyplot as plt
import ubem2d as ubem

def solve_cylinder_with_circulation(uinf, circ, cyl, axis, At, An):
    # Solve flow problem
    tx,ty,nx,ny = cyl.tx, cyl.ty, cyl.nx, cyl.ny
    (x0,y0) = cyl.center
    (u,v) = ubem.velocity_vortex(circ,x0,y0,cyl.xmid,cyl.ymid)
    rhs = -(uinf[0]*nx + uinf[1]*ny + u*nx + v*ny)
    sigma = nla.solve(An,rhs)

    # Compute aerodynamic forces and moments
    qt = At.dot(sigma) + uinf[0]*tx + uinf[1]*ty + u*tx + v*ty
    cp = 1. - (qt/nla.norm(uinf))**2
    (CD,CL,CM) = ubem.body_cdclcm(uinf, cyl, cp, *axis, ubem.Orientation.CCW)

    return (sigma,qt,cp,CD,CL,CM)

def solve_and_plot(uinf, circ, cyl, axis, At, An, X, Y, I):
    (sigma,qt,cp,CD,CL,CM) = solve_cylinder_with_circulation(uinf, circ, cyl,
        axis, At, An)

    # Define/display mesh; show aerodynamic coefficients
    print('Circulation:',circ)
    print('Drag coefficient, CD:',CD)
    print('Lift coefficient, CL:',CL)
    print('Theoretical (Kutta-Joukowsi) lift coefficient:',-circ)
    print('Moment coefficient, CM:',CM)
    print('Moment axis:',axis)

    # Directory in which to save plots
    plotdir = os.path.dirname(os.path.abspath(__file__))

    # Plot pressure on boundary.  For the theoretical solution, refer to
    # section 4.9 of Currie's book "Fundamental Mechanics of Fluids", 2nd Ed.
    theta = np.arctan2(cyl.ymid,cyl.xmid)
    theta_exact = np.linspace(-np.pi,np.pi,100)
    ii = theta.argsort()
    plt.figure()
    plt.plot(theta[ii],cp[ii],'o')
    qt_exact = -2*nla.norm(uinf)*np.sin(theta_exact)+.5*circ/(np.pi)
    cp_exact = 1.-(qt_exact/nla.norm(uinf))**2
    plt.plot(theta_exact, cp_exact,'k-')
    plt.xlabel('Angle (radians)')
    plt.ylabel('Pressure coefficient (Cp)')
    plt.legend(['Panel method','Theoretical solution'])
    plt.title('Flow past cylinder; circulation = {:.3f}'.format(circ))
    plt.grid(True)
    plt.savefig(os.path.join(plotdir,'cylinder_cp_boundary.pdf'))

    # Plot pressure field on mesh
    cylcolor='.75'
    (x0,y0) = cyl.center
    (Us,Vs) = ubem.velocity_source_body(cyl, sigma, X, Y)
    (Uv,Vv) = ubem.velocity_vortex(circ,x0,y0,X,Y)
    U = Us + Uv + uinf[0]
    V = Vs + Vv + uinf[1]
    Z = 1. - (U**2+V**2)/nla.norm(uinf)**2
    plt.figure()
    lvls = np.linspace(np.min(cp_exact),np.max(cp_exact),100)
    plt.contourf(X,Y,Z,lvls)
    plt.fill(cyl.x,cyl.y,color=cylcolor)
    plt.colorbar()
    plt.axis('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Pressure field (Cp); circulation = {:.3f}'.format(circ))
    plt.savefig(os.path.join(plotdir,'cylinder_cp_field.pdf'))

    # Contour plot of stream function
    plt.figure()
    Z = ubem.sf_source_body(cyl,sigma,X,Y)
    Z += ubem.sf_vortex(circ,0,0,X,Y)
    Z += ubem.sf_uniform_flow(uinf,X,Y)
    interior_nodes = np.where(I == 0)
    exterior_nodes = np.where(I == 1)
    Z[interior_nodes] = 0
    m,M = np.min(Z[exterior_nodes]), np.max(Z[exterior_nodes])
    lvls = np.linspace(m,M,50)
    lvls = lvls[np.where(np.abs(lvls) > 0.01*(M-m))]
    CS = plt.contourf(X,Y,Z,lvls,cmap=matplotlib.colormaps['winter'])
    plt.contour(CS, levels=CS.levels[::2],colors='k')
    xy = np.array([[x,y] for x in X[0,:] for y in Y[:,0]])
    plt.plot(xy[:,0],xy[:,1],'o',markersize=.25)
    plt.fill(cyl.x,cyl.y,color=cylcolor)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Stream lines; circulation = {:.3f}'.format(circ))
    plt.axis('equal')
    plt.savefig(os.path.join(plotdir,'cylinder_stream_lines.pdf'))

    # Show all plots and block
    plt.show()

if __name__ == '__main__':
    # Set up and solve flow problem
    npan = 50
    uinf = (1,0)  # Onset flow
    cyl = ubem.CircularCylinder(npan)
    axis = cyl.center
    (At,An) = ubem.source_influence_matrices_body(cyl)
    (X,Y) = ubem.mesh(cyl,50,50,1.5)
    I = ubem.boundary(cyl,X,Y)

    circulation = 0
    solve_and_plot(uinf, circulation, cyl, axis, At, An, X, Y, I)
