import os
import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import ubem2d as ubem

if __name__ == '__main__':
    bodies = [
        ubem.naca4('0012',50),
        ubem.naca4('0012',50).heave(1),
        ubem.naca4('0012',50).heave(1).surge(1.5),
        ubem.naca4('0012',50).surge(1.5)
    ]
    sys = ubem.HessSmithSystem(bodies)

    # Solve for steady flow past the system of bodies
    uinf = (1,0)
    soln = sys.solve(uinf)
    
    # Plot pressure distribution along body boundaries
    cp = sys.pressure_self(uinf, soln)
    plt.figure()
    for k,body in enumerate(bodies):
        plt.plot(body.xmid-body.x[body.le], cp[k])
    plt.gca().invert_yaxis()
    plt.legend(['Body {}'.format(k) for k in range(len(bodies))])
    plt.xlabel('$x/c$')
    plt.ylabel('Pressure coefficient, $C_p$')
    plt.title('Pressure distribution')
    plt.grid(True)

    # Plot pressure field on a mesh
    X,Y = ubem.mesh(bodies,100,100,.25)
    U,V = sys.flow_external(uinf, soln, X, Y)
    Cp = sys.pressure_from_flow(uinf, soln, U, V)
    plt.figure()
    plt.contourf(X,Y,Cp,100)
    for body in bodies:
        plt.fill(body.x, body.y)
    plt.colorbar()
    plt.legend(['body {}'.format(k) for k in range(len(bodies))])
    plt.axis('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Pressure field')

    # Stream plot of velocity field
    plt.figure()
    sp = plt.streamplot(X[0,:], Y[:,0], U, V, color='k', density=1,
        linewidth=.5)
    for body in bodies:
        plt.fill(body.x, body.y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Stream lines of flow')
    plt.axis('equal')

    # Print force and moment coefficients
    for k,body in enumerate(bodies):
        cd,cl,cm = ubem.airfoil_cdclcm(uinf, body, cp[k])
        print('body {:>2d}, CD={:>9.3g}, CL={:>9.3g}, CM(LE)={:>9.3g}'.format(k,
            cd,cl,cm))

    plt.show()
