'''

Ray tracing algorithm

'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def ray_tracing(w=400,h=300):
    R = np.linspace(0,10,3)
    FOV = 90
    FOV = np.deg2rad(FOV)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in np.linspace(-FOV/2,FOV/2,w):
        for j in np.linspace(-FOV/2,FOV/2,h):
            x = R*np.cos(i)*np.cos(j)
            y = R*np.sin(i)*np.cos(j)
            z = R*np.sin(j)
            ax.plot(x,y,z)
            print (i,j)
    ax.view_init(0, 90)
    plt.savefig('fov_side.png',dpi=300)

ray_tracing(10,8)
