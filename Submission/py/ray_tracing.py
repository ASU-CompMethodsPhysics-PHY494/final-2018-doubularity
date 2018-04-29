'''

Ray tracing algorithm

'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import misc
import time
import curses

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

def intersect_plane(O, D, P, N):
    # Return the distance from O to the intersection of the ray (O, D) with the
    # plane (P, N), or +inf if there is no intersection.
    # O and P are 3D points, D and N (normal) are normalized vectors.
    denom = np.dot(D, N)
    if np.abs(denom) < 1e-6:
        return np.inf
    d = np.dot(P - O, N) / denom
    if d < 0:
        return np.inf
    return d

def trace(theta,phi):
    R = 1
    x = R*np.cos(theta)*np.cos(phi)
    y = R*np.sin(theta)*np.cos(phi)
    z = R*np.sin(phi)
    t = np.inf
    rayO = CAMPOS
    rayD = np.array((x,y,z))
    t_obj = intersect_plane(rayO, rayD, obj.r,obj.normal)
    if t_obj < t:
        t = t_obj
    if t == np.inf:
        return (np.zeros(3),np.zeros(3))
    M = rayO + rayD * t
    if (obj.p1[0] <= M[1] <= obj.p2[0]) and (obj.p1[1] <= M[2] <= obj.p2[1]):
        color = obj.get_color(M[1],M[2])
    else:
        color = np.zeros(3)
    return (color, M)



def ray_cast(w=640,h=480,FOV_w=np.deg2rad(40),FOV_h=np.deg2rad(30)):
    img = np.zeros((h,w,3))
    pix = w*h
    count = 0
    for i,phi in enumerate(np.linspace(-FOV_h/2,FOV_h/2,h)):
        for j,theta in enumerate(np.linspace(-FOV_w/2,FOV_w/2,w)):
            color, (x, y, z) = trace(theta,phi)
            img[i][j] = color
            if count/pix*100%5 == 0:
                print ('{}%...'.format(count/pix*100))
            count += 1
    print ('100%\nDone.')
    plt.imshow(img,interpolation='nearest')
    plt.title('Single Image Ray Tracing')
    print ('Saving Image')
    plt.savefig('../figures/ray_traced_img.png',dpi=1000)
    
obj = Plane()
R = np.linspace(0,5,5)
CAMPOS = [-10,0,0]
FOV_w = np.deg2rad(40)
FOV_h = np.deg2rad(30)
RES = [960,540]
ray_cast(w=RES[0],h=RES[1])
