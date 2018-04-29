import numpy as np
from ode import RK4f, sqrnorm
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

plt.style.use('ggplot')

def trace_rays(theta,phi,h=0.01):
    phi = np.linspace(-np.deg2rad(phi),np.deg2rad(phi),4)
    theta = np.linspace(-np.deg2rad(theta),np.deg2rad(theta),4)
    x_p = np.cos(theta)*np.cos(phi)
    y_p = np.sin(theta)*np.cos(phi)
    z_p = np.sin(phi)
    velocities = np.vstack(np.meshgrid(x_p,y_p,z_p)).reshape(3,-1).T

    fig = plt.figure()
    ax = plt.axes(projection='3d')

    for velocity in velocities:
        point = np.array([-3.0,0.0,0.0])
        h2 = sqrnorm(np.cross(point,velocity))
        pos = []
        for i in range(500):
            #simple step size control
            rkstep = 0.02

            # standard Runge-Kutta
            y = np.zeros(6)
            y[0:3] = point
            y[3:6] = velocity
            k1 = RK4f( y, h2)
            k2 = RK4f( y + 0.5*h*k1, h2)
            k3 = RK4f( y + 0.5*h*k2, h2)
            k4 = RK4f( y + h*k3, h2)

            increment = rkstep/6. * (k1 + 2*k2 + 2*k3 + k4)
            if (np.linalg.norm(increment[3:6])) > 5:
                break
            velocity += increment[3:6]

            point += increment[0:3]
            pos.append(y[0:3])


        pos = np.array(pos)
        ax.plot3D(pos[:,0],pos[:,1],pos[:,2])
    plt.xlim(-3,1)
    plt.ylim(-2,2)
    ax.set_zlim(-2,2)
    plt.title('Ray Tracing Black Hole')
    plt.savefig('../figures/ray_tracing_bh_3d.png',dpi=300)
