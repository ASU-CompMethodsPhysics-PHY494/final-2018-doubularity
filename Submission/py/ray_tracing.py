'''

Ray tracing algorithm

'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import time
from ode import rk4

class Plane:
    '''

    Parameters
    ----------
    path : string
      Relative path to an image being added to the scene

    r : array
      Position of the image in spherical coordinates

    '''
    def __init__(self,path='../figures/sagrada-familia.jpg',r=[10,0,0],normal=[1,0,0],scale=10):
        self.path = path
        self.r = np.array(r)
        self.normal = np.array(normal)
        self.img = misc.imread(self.path)/255
        ratio = self.img.shape[0]/self.img.shape[1]
        w = scale
        h = scale*ratio
        self.p1 = np.array((self.r[1]-w/2,self.r[2]-h/2))
        self.p2 = np.array((self.r[1]+w/2,self.r[2]+h/2))
        self.p_x = np.linspace(self.p1[0],self.p2[0],self.img.shape[1])
        self.p_y = np.linspace(self.p1[1],self.p2[1],self.img.shape[0])

    def find_nearest(self,array,value):
        idx = (np.abs(array-value)).argmin()
        return (idx)

    def get_color(self,x,y):
        x_ind = self.find_nearest(self.p_x,x)
        y_ind = self.find_nearest(self.p_y,y)
        return (self.img[y_ind][x_ind])

def sqrnorm(vec):
    return np.einsum('...i,...i',vec,vec)

def RK4f(y,h2):
    f = np.zeros(y.shape)
    f[0:3] = y[3:6]
    f[3:6] = - 1.5 * h2 * y[0:3] / np.power(sqrnorm(y[0:3]),2.5)
    return f

def trace_ray(pos,theta,phi,h=0.1):
    v_x = np.cos(theta)*np.cos(phi)
    v_y = np.sin(theta)*np.cos(phi)
    v_z = np.sin(phi)
    velocity = np.array([v_x,v_y,v_z])
    point = np.array(pos)
    h2 = sqrnorm(np.cross(point,velocity))

    while np.linalg.norm(point) <= 15:
        y = np.zeros(6)
        y[0:3] = point
        y[3:6] = velocity
        increment = rk4(y,RK4f,h2,h)
        color = np.zeros(3)
        if np.linalg.norm(point+increment[0:3]) < 0.01:
            break
        if ((obj.p1[0] <= point[1]+increment[1] <= obj.p2[0]) and
             (obj.p1[1] <= point[2]+increment[2] <= obj.p2[1]) and
             (obj.r[0]-0.1 <= point[0]+increment[0] <= obj.r[0]+0.1)):
            color = obj.get_color(point[1]+increment[1],point[2]+increment[2])
            break
        point += increment[0:3]
        velocity += increment[3:6]
    return (color)


def ray_cast(w=640,h=480,cam=[-10.0,0.0,0.0],quad=None):
    FOV_w=np.deg2rad(40)
    FOV_h=np.deg2rad(30)
    img = np.zeros((h,w,3))
    pix = w*h
    count = 0
    if quad == 1:
        t_ang = np.linspace(0,FOV_h/4,h/2)
        p_ang = np.linspace(0,FOV_w/4,w/2)
    elif quad == 2:
        t_ang = np.linspace(0,FOV_h/4,h/2)
        p_ang = np.linspace(-FOV_w/4,0,w/2)
    elif quad == 3:
        t_ang = np.linspace(-FOV_h/4,0,h/2)
        p_ang = np.linspace(0,FOV_w/4,w/2)
    elif quad == 4:
        t_ang = np.linspace(-FOV_h/4,0,h/2)
        p_ang = np.linspace(-FOV_w/4,0,w/2)
    else:
        t_ang = np.linspace(-FOV_h/2,FOV_h/2,h)
        p_ang = np.linspace(-FOV_w/2,FOV_w/2,w)
    for i,phi in enumerate(t_ang):
        for j,theta in enumerate(p_ang):
            color = trace_ray(cam,theta,phi)
            img[i][j] = color
            if count/pix*100%5 == 0:
                print ('{}%...'.format(count/pix*100))
            count += 1
    print ('100%\nDone.')
    return (img)
    '''
    plt.imshow(img,interpolation='nearest')
    plt.title('Single Image Ray Tracing')
    print ('Saving Image')
    plt.show()
    '''

from multiprocessing.dummy import Pool as ThreadPool

obj = Plane()
CAMPOS = [-10.0,0.0,0.0]
FOV_w = np.deg2rad(40)
FOV_h = np.deg2rad(30)
RES = [16,9]
ray_cast(w=RES[0],h=RES[1])

in_vals = [[RES,CAMPOS,1],[RES,CAMPOS,2],
            [RES,CAMPOS,3],[RES,CAMPOS,4]]

pool = ThreadPool(4)
results = zip(*pool.map(ray_cast, in_vals))
print (results)
