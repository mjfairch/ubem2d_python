import numpy as np
import numpy.linalg as nla
import matplotlib.pyplot as plt
import ubem2d as ubem

if __name__ == '__main__':
    # Set up and solve flow problem
    uinf = (1,0)  # Onset flow
    body = ubem.Rectangle(1,.5,20,10).rotate(np.pi/6)
    soln = ubem.solve_source_body(uinf, body)

    # Compute aerodynamic forces and moments
    cp = 1. - (soln.qt/nla.norm(uinf))**2
    CD,CL,CM = ubem.body_cdclcm(uinf, body, cp, *body.center)

    # Define/display mesh; show aerodynamic coefficients
    (X,Y) = ubem.mesh(body, 40, 40, .5, equal=False)
    I = ubem.boundary(body, X, Y)
    print(ubem.boundary_to_string(I))
    print('Drag coefficient, CD:', CD)
    print('Lift coefficient, CL:', CL)
    print('Moment coefficient, CM:', CM)
    print('Moment axis:', body.center)

    # Stream plot of velocity field
    plt.figure()
    Us,Vs = ubem.velocity_source_body(body, soln.sigma, X, Y)
    U,V = Us+uinf[0], Vs+uinf[1]
    plt.streamplot(X[0,:], Y[:,0], U, V, density=1, color='0', linewidth=.5)
    plt.plot(body.x, body.y, 'k-')
    plt.fill(body.x, body.y, color='.75')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Stream plot of velocity field')
    plt.axis('equal')

    # Show all plots and block
    plt.show()
