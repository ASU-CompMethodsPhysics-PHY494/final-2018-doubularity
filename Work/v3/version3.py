import matplotlib.pyplot as plt
import numpy as np

def render(res,angle,R,D,Z):
    """Renders the image in the defined space
    res   - <list> resolution, number of x/y pixels
    angle - <float> angle covered by x distance
    R     - <float> radius of the "lens" of the camera
    D     - <float> distance of camera from center
    Z     - <float> Z position of the object
    """
    t = 0

    def inShape(r,l,Z,shape = "circle",thickness = 1):
        """The function checking if the coordinates are in the defined shape
        r         - <list> coordinates of the particle
        l         - <float> the standard length of the shape
        Z         - <float> Z position of the shape plane
        shape     - <string> name of the target shape
        thickness - <float> the thickness of the target plane
        """

        #We start by converting the position to cartesian coordinates
        x = r[0] * np.sin(r[1]) * np.cos(r[2])
        y = r[0] * np.sin(r[1]) * np.sin(r[2])
        z = r[0] * np.cos(r[1])

        inside = False

        if abs(z - Z) <= thickness:
            if shape == "circle":
                if np.sqrt(x**2 + y**2) <= l:
                    inside = True
            if shape == "square":
                if -l < 2*x < l and -l < 2*y < l:
                    inside = True

        return inside

    def f(t,q):
        """The derivative of the position-velocity list, q
        t - <float> time
        q - <list> position-velocity list [r,theta,phi,r_dot,M,L]
        """
        a_r = 1/q[0]**3 * (q[4]**2 + q[5]**2/np.sin(q[1])**2)
        M_dot = np.cos(q[1])/np.sin(q[1])**3 * q[5]**2/q[0]**2
        L_dot = 0

        r_dot = q[3]
        theta_dot = q[4]/q[0]**2
        phi_dot = q[5]/(q[0]*np.sin(q[1]))**2

        return [r_dot,theta_dot,phi_dot,a_r,M_dot,L_dot]

    #setting up camera properties
    max_alpha = angle/2
    max_beta = angle/2 * res[1]/res[0]
    alpha = np.arange(-max_alpha,max_alpha,2*max_alpha/res[0])
    beta = np.arange(-max_beta,max_beta,2*max_beta/res[1])

    for a in alpha:
        for b in beta:
            #setting up initial particle position/velocity
            r0 = np.sqrt(R**2 + D**2 - 2*R*D*np.sqrt(1 - (np.sin(a)**2 + np.sin(b)**2)/2))
            
