import numpy as np

global raNorth, decNorth, raSouth, decSouth, radius
raNorth, decNorth = 270., 66.56070889
raSouth, decSouth = 90., -66.56070889
radius = 5.

def is_star_in_CVZ(alpha_deg, delta_deg):
    alpha_deg, delta_deg = np.ascontiguousarray(alpha_deg), \
			   np.ascontiguousarray(delta_deg)
    # Check Northern CVZ
    draN  = alpha_deg - raNorth
    ddecN = delta_deg - decNorth
    distfromCVZcentreN = np.sqrt(draN**2 + ddecN**2)
    
    # Check Southern CVZ
    draS  = alpha_deg - raSouth
    ddecS = delta_deg - decSouth
    distfromCVZcentreS = np.sqrt(draS**2 + ddecS**2)

    # is star in a CVZ?
    inCVZ = np.zeros(alpha_deg.size)
    inCVZ[(distfromCVZcentreN <= radius) | (distfromCVZcentreS <= radius)] = 1.
    return inCVZ
