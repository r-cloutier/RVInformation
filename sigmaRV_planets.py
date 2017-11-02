import numpy as np
import rvs


def draw_FGK_planet(Pmin, nplanets):
    '''
    Draw planets from around FGK stars from the occurrence rates derived in 
    Fressin et al 2013. Ensure that each pair in the multi-planet system is 
    Lagrange stable given circular orbits.
    '''
    assert nplanets > 1
    Pgrid = np.array([.8,2.,3.4,5.9,10.,17.,29.,50.,85.,145.,245.,418.])
    rpgrid = np.array([.8,1.25,2.,4.,6.,22.])
    occurrence = np.array([[.015,.067,.17,.18,.27,.23,.35,.71,1.25,.94,1.05],
                           [.004,.006,.11,.091,.29,.32,.49,.66,.43,.53,.24],
                           [.035,.18,.73,1.93,3.67,5.29,6.45,5.25,4.31,3.09,1.29],
                           [.17,.74,1.49,2.90,4.3,4.49,5.29,3.66,6.54,4.65,3.01],
                           [.18,.61,1.72,2.7,2.7,2.93,4.08,3.46,4.55,5.25,4.3]]).T
    occurrence /= norm

    Ps, mps = np.zeros(nplanets-1), np.zeros(nplanets-1)
    for i in range(nplanets-1):
        for j in range(Pgrid.size-1):
            for k in range(rpgrid.size-1):
                if np.random.rand() <= occurrence[j,k]:
                    
    


# rvs.is_Lagrangestable(Ps, Ms, mps, eccs)
def get_sigmaRV_planets(P, Ms, mult):
    return None
