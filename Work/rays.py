import ode
import numpy as np

def runTrace(r0,v0,L,M,R = 150.,Z = -100.):
    """This runs light traces and finds their end position
    r0 - <array> of [r,theta,phi] starting positions
    v0 - <float> of the initial radial velocity
    L - <float> of the azimuthal momentum
    M - <float> the initial angular momentum
    R - <float> the kill distance for the light rays
    Z - the Z position of the test plane
    """

def cam(res,angle,R = 1,Z = 100):
    """This sets up the camera for the scene
    res - <array> number of [x,y] pixels
    R - <float> proportional to the focal length of the camera
    angle - <array> angle of the [x,y] dimensions
    Z - <float> Height of camera from the central coordinate
    """

    def shape(P,Z,t = 1,l = 5):
        """It's an F"""

        inPlane = False

        r = P[1]
        theta = P[2]
        phi = P[3]

        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)

        if abs(z-Z) <= t:
            if -1.5<x/l<1.5 and -2.5<y/l<2.5:
                return True
            elif -1.5<x/l<1.5 and 1.5<y/l<2.5:
                return True
            elif -0.5<x/l<0.5 and 0.5<y/l<1.5:
                return True
            else:
                return False
        else:
            return False
    
    p = np.arange(-angle[0],angle[0],res[0])
    q = np.arange(-angle[1],angle[1],res[1])

    CAM = np.zeros((len(p),len(q)))
    
    n = m = -1
    for alpha in p:
        n += 1
        for beta in q:
            m += 1
            r0 = np.sqrt(2*R**2*(np.sin(alpha)**2 + np.sin(beta)**2) + Z**2 - 2*R*Z*np.sqrt(np.sin(alpha)**2 + np.sin(beta)**2))
            v0 = 
            Ps = runTrace(r0,v0,0,M)
            if inShape(Ps[-1]):
                CAM[alpha,beta] = 1
