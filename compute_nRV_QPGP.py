from imports import *
from sympy.matrices import *

t = np.linspace(0,100,101)
theta = 1, 100, .5, 10, 0

def covariance_matrix(theta, t, sig=0):
    a, l, G, P, s = theta
    K = np.zeros ((t.size, t.size))
    K += -.5/l**2 * (np.tile(t, (t.size,1)) - \
                     np.tile(t, (t.size,1)).T)**2  # :)
    K += -G**2 * np.sin(np.pi/P * abs(np.tile(t, (t.size,1)) - \
                                      np.tile(t, (t.size,1)).T))**2
    K = a*a*np.exp(K)
    return K
