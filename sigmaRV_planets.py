import numpy as np
import rvs
import massradius as mr


def draw_FGK_planet(P, rp, nplanets, Ms):
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
    occurrence *= 1e-2

    # sample planets until we get the correct multiplicity
    Ps, mp = np.array([P]), mr.radF2mass(rp)
    while Ps.size < nplanets:
        
        Ps, rps = np.array([P]), np.array([rp])
        mps = mr.radF2mass(rps)
        for i in range(Pgrid.size-1):
            for j in range(rpgrid.size-1):
                
                if np.random.rand() <= occurrence[i,j] and Pgrid[i] > P:

                    # draw a system
                    Ps  = np.append(Ps, np.random.uniform(Pgrid[i], Pgrid[i+1]))
                    rps = np.append(rps, np.random.uniform(rpgrid[j], rpgrid[j+1]))
                    mps = mr.radF2mass(rps)

        sort = np.argsort(Ps)
        Ps, rps, mps = Ps[sort], rps[sort], mps[sort]

        # get stable pairs
        if Ps.size >= nplanets:
            lagrangestable = rvs.is_Lagrangestable(Ps, Ms, mps,
                                                   np.zeros(Ps.size))
            if lagrangestable[-1] == 1:
                lagrangestable = np.append(lagrangestable.astype(bool), True)
            else:
                lagrangestable = np.append(lagrangestable.astype(bool), False)
            Ps, rps, mps = Ps[lagrangestable], rps[lagrangestable], \
                           mps[lagrangestable]

    # compute semi-amplitudes
    Ks = rvs.RV_K(Ps, Ms, mps)
            
    return Ps, rps, mps, Ks


def _compute_sigmaRV_planets(Ks):
    rms = np.sin(np.linspace(0,2*np.pi,100)).std()
    return np.sqrt(np.sum((Ks*rms)**2))


# rvs.is_Lagrangestable(Ps, Ms, mps, eccs)
def get_sigmaRV_planets(P, rp, Ms, mult):
    _,_,_,Ks = draw_FGK_planet(P, rp, mult, Ms)
    return _compute_sigmaRV_planets(Ks)

