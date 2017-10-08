import os
import numpy as np
import matplotlib.pyplot as plt
import ubem2d as ubem

# Create a triangular body whose motion we'll specify and graph
body = ubem.BrokenLine([-1,1,0,-1],[-1,-1,1,-1]).scale(.1).center()

# Define the rigid motion (i.e. path in SE(2))
ftheta = lambda t: t
fx = lambda t: np.cos(t)
fy = lambda t: np.sin(t)
motion = ubem.RigidMotion(ftheta, fx, fy)
t = np.linspace(0, 2*np.pi, 17)

# Move body to its initial position and plot it
body.glide(motion(t[0]))

# To update from g to h in SE(2), act on the left by h*g.inv():
update = lambda t0,t: motion(t)*(motion(t0).inv())

# March body along the specified path in SE(2), plotting at each step
plt.figure()
plt.plot(body.x, body.y)
for i in range(1, len(t)):
    body.glide(update(t[i-1], t[i]), 0, 0)
    plt.plot(body.x, body.y)
    plt.text(fx(t[i]), fy(t[i]), str(i))
tt = np.linspace(t[0], t[-1], 100)
plt.plot(fx(tt), fy(tt))
plt.axis('equal')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('SE(2) kinematics')
plt.savefig(os.path.join(ubem.__plot_dir, 'rigid_motion.pdf'))

# Show all plots and block
plt.show()
