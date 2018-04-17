import ode
import numpy as np

# y = [r,theta,phi,dr/dt,M]
# r0 -> (r,theta,phi) coordinates for starting position
# v0 -> radial velocity at starting position
# L -> azimuthal angular momentum
# M -> angular momentum
# m -> mass of the black hole
# h -> time step

def tracePath(r0,v0,L,M,m,h = 0.01):
    y = r0
    y.append(v0)
    y.append(M)

    def f(t,y):
        ar_rad = 2*m/(y[0]-2*m)*y[3]**2/y[0] #The acceleration due to the radial velocity
        ar_cent = (r-3m)*(r-2*m)**2/y[0]**6*(M**2+L**2/np.sin(y[1])**2) #The acceleration due to the angular velocity
        ar = ar_rad + ar_cent

        M^dot = np.cos(y[1])/np.sin(y[1])**3 * (y[0]-2*m)/y[0]**3) * L**2 #The change in angular momentum about the theta angle

        theta^dot = (y[0]-2*m)/y[0]**3*y[4]
        phi^dot = (y[0]-2*m)/y[0]**3*y[4]/np.sin(y[1])**2

        return [y[3],theta^dot,phi^dot,ar,M^dot]
        
