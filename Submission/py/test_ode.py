import numpy as np
from ode import rk4
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

plt.style.use('ggplot')

def sqrnorm(vec):
    return np.einsum('...i,...i',vec,vec)

def RK4f(y,h2):
    f = np.zeros(y.shape)
    f[0:3] = y[3:6]
    f[3:6] = - 1.5 * h2 * y[0:3] / np.power(sqrnorm(y[0:3]),2.5)
    return f

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
            y = np.append(point,velocity)
            y += rk4(y,RK4f,h2,h)
            if np.linalg.norm(y[3:6]) > 5:
                break
            pos.append(y[0:3])
        pos = np.array(pos)
        ax.plot3D(pos[:,0],pos[:,1],pos[:,2])
    plt.xlim(-3,1)
    plt.ylim(-2,2)
    ax.set_zlim(-2,2)
    plt.show()
    #plt.title('Ray Tracing Black Hole')
    #plt.savefig('../figures/ray_tracing_bh_3d.png',dpi=300)

trace_rays(30,30,h=0.001)
