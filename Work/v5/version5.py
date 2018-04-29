import numpy as np
import matplotlib.pyplot as plt

def cart2Pol(x):
    """Translates cartesian coordinats to spherical coordinates
    x - <array> the [x,y,z] coordinates
    """
    r = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    if r > 0:
        theta = np.arccos(x[2]/r)
    else:
        theta = 0
    if theta != 0:
        phi = np.arccos(x[0]/(r*np.sin(theta)))
    else:
        phi = 0

    return r,theta,phi

def pol2Cart(r):
    """Translates spherical coordinates to cartesian coordinates
    r - <array> the [r,theta,phi] coordinates
    """
    x = r[0] * np.sin(r[1]) * np.cos(r[2])
    y = r[0] * np.sin(r[1]) * np.sin(r[2])
    z = r[0] * np.cos(r[1])

    return x,y,z

def verlet(y,f,t,h):
    """integrates the function, returns numpy array of y
    y         - <numpy.array> [r,theta,phi,r_dot,M,L]
    f         - <numpy.array> derive function of y
    t         - <float> time
    h         - <float> time step
    """
    r_0 = np.array(y[:3])
    v_0 = np.array(f(t,y)[:3])

    v_1_2 = v_0 + h/2 * f(t,y)[3:]
    r_1 = r_0 + h * v_1_2

    y1 = np.append(r_1,v_0)

    v_1 = v_1_2 + h/2 * f(t + h,y1)[3:]

    return np.append(r_1,v_1)

def render(res = [16,9],angle = 30.,D = 50.,R = 5.,Z = -10.,thiccness = 2.,shape = "square",l = 65.,kR = 100,m = 0,h = 0.01):
    """Renders the setup environment
    res       - <array> [x,y] resolution of the image
    angle     - <float> angle covered by the x axis of the image
    D         - <float> distance of camera from center
    Z         - <float> z position of the collision plane
    thiccness - <float> thickness of the collision plane
    shape     - <string> which shape is being used
    l         - <float> standard length of shape
    kR        - <float> maximum particle distance
    m         - <float> mass of black hole
    h         - <float> time step length
    """

    def inShape(r,Z,thiccness,shape,l):
        """Returns True of False if the particle is inside the shape
        r         - <numpy.array> [r,theta,phi] position of the particle
        Z         - <float> z position of the collision plane
        thiccness - <float> thickness of the collision plane
        shape     - <string> which shape is being used
        l         - <float> standard length of shape
        """
        x,y,z = pol2Cart(r)

        inside = False
        if abs(Z - z) <= thiccness:
            if shape == "circle":
                if x**2 + y**2 <= l**2:
                    inside = True
            elif shape == "square":
                if -l<x<l and -l<y<l:
                    inside = True

        return inside

    def f(t,y):
        """Derives the state vector y = [r,theta,phi,r_dot,theta_dot,phi_dot]
        t - <float> time
        y - <numpy.array> state vector
        """
        a_r = y[0]*(y[4]**2 + np.sin(y[1])**2 * y[5]**2)
        a_theta = np.sin(y[1])*np.sin(y[1])*y[2]**2 - 2*y[3]*y[4]/y[0]
        a_phi = -2*(np.cos(y[1])/np.sin(y[1])*y[4]*y[5] + y[3]*y[5]/y[1])

        return np.array([y[3],y[4],y[5],a_r,a_theta,a_phi])

    #Setting up the camera:
    alpha = angle * np.pi/180
    beta = alpha * res[1]/res[0]
    CAM = np.zeros((res[1],res[0]))

    per = -1/(res[0]*res[1])
    m = -1
    for a in np.arange(-alpha/2,alpha/2,alpha/res[0]):
        m += 1
        n = -1
        for b in np.arange(-beta/2,beta/2,beta/res[1]):
            per += 1/(res[0]*res[1])
            n += 1

            #Setting up particle position:
            x0 = [R*np.sin(a),R*np.cos(a)*np.sin(b),D - R*np.cos(a)*np.cos(b)]
            r0 = np.asarray(cart2Pol(x0))
            #Setting up particle velocity:
            r0_dot = np.sqrt(1 - (np.sin(a)**2 + np.cos(a)**2 * np.sin(b)**2))
            theta0_dot = np.sqrt(x0[1]**2 + x0[2]**2)/r0[0]
            v0 = np.array([r0_dot,theta0_dot,0])

            y = np.append(r0,v0)

            t = 0
            positions = [[t,y[0],y[1],y[2]]]
            while 1:
                y = verlet(y,f,t,h)
                y += h
                positions.append([t,y[0],y[1],y[2]])

                if y[0] > kR:
                    break
                elif y[0] <= 3*m:
                    break
                elif inShape(y[:3],Z,thiccness,shape,l):
                    CAM[n,m] += 1
                    break
            print(pol2Cart(positions[0][1:])[2],pol2Cart(positions[-1][1:])[2])
            CAM[n,m] += 1
        print(100*per,"%")

    plt.imshow(CAM)
    plt.show()
