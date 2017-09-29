'''
Compute nRV for all TESS targets using the code cloned from
https://github.com/r-cloutier/Transiting_RVcalculator.git
'''
import numpy as np
from get_tess_data import get_TESS_data
from compute_sigmaRV import compute_sigmaRV, get_reduced_spectrum
import pylab as plt
import rvs
from uncertainties import unumpy as unp


global G, rhoEarth, R, aperture, QE
G, rhoEarth, R, aperture, QE = 6.67e-11, 5.51, 75e3, 3.58, .1


def estimate_nRV(planetindex, band_strs, Z=0, R=75e3):
    '''
    Estimate the number of RVs to measure the mass of a transiting planet at 
    a given detection signifiance.
    '''
    # Read-in TESS data for this planetary system
    ra,dec,rp,P,S,K,Rs,Teff,Vmag,Imag,Jmag,Kmag,dist,_,_,snr,_ = get_TESS_data()
    nplanets = ra.size
    assert 0 <= planetindex < nplanets
    rp, P, S, K, Rs, Teff, Vmag, Imag, Jmag, Kmag, dist = rp[planetindex], \
                                                          P[planetindex], \
                                                          S[planetindex], \
                                                          K[planetindex], \
                                                          Rs[planetindex], \
                                                          Teff[planetindex], \
                                                          Vmag[planetindex], \
                                                          Imag[planetindex], \
                                                          Jmag[planetindex], \
                                                          Kmag[planetindex], \
                                                          dist[planetindex]

    mp = get_planet_mass(rp)
    Ms = get_stellar_mass(P, mp, K)
    logg = np.log10(G * rvs.Msun2kg(Ms) / rvs.Rsun2m(Rs)**2 * 1e2)

    # round Teff and logg
    Teffs = np.append(np.arange(23e2,7e3,1e2), np.arange(7e3,121e2,2e2))
    Teff_round = Teffs[abs(Teffs-Teff) == np.min(abs(Teffs-Teff))]
    loggs = np.arange(0,6.1,.5)
    logg_round = loggs[abs(loggs-logg) == np.min(abs(loggs-logg))]

    # compute sigmaRV in each band
    sigmaRVs = np.zeros(len(band_strs))
    for i in range(sigmaRVs.size):
    
        wl, spec = get_reduced_spectrum(Teff_round, logg_round, Z, vsini,
                                        band_strs[i], R)
        texp = exposure_time_calculator_per_band(mags[i], band_strs[i],
                                                 aperture_m, QE, R)
        sigmaRVs[i] = compute_sigmaRV(wl, spec, mags[i], band_strs[i], texp)
    
    # Compute nRV
    sigmaRV = 1. / np.sqrt(np.sum(1./sigmaRVs**2))
    sigmaK_target = get_sigmaK_target_v1(mp, rp)
    N = 2 * (sigmaRV / sigmaK_target)**2
    
    return N



def get_planet_mass(rp):
    '''Compute the TESS planet mass from its reported radius.'''
    if rp < 1.5:
        return .44*rp**3 + .614*rp**4
    else:
        return 2.69*rp**.93
    

def get_stellar_mass(P, mp, K):
    '''Compute the stellar mass from the orbital period, planet mass, and RV 
    semi-amplitude.'''
    return rvs.kg2Msun(np.sqrt(2*np.pi*G / rvs.days2sec(P) * \
                               (rvs.Mearth2kg(mp) / K)**3))


def get_sigmaK_target_v1(rp, K, P, Ms, sigP=5e-5,
                         fracsigrp=1e-1, fracsigMs=1e-1, fracsigrho_target=.2):
    '''Compute the K detection significance that must be achieved in order to 
    measure a planet's density at a specific precision.'''
    urp, uP, uMs = unp.uarray(rp, rp*fracsigrp), unp.uarray(P, sigP), \
                   unp.uarray(Ms, Ms*fracsigMs)

    fracsigKs = np.logspace(-2,1,100)
    fracsigrho = np.zeros(fracsigKs.size)
    
    for i in range(fracsigKs.size):
        uK = unp.uarray(K, K*fracsigKs[i])
        ump = rvs.RV_mp(uP, uMs, uK)
        urho = rhoEarth * ump / urp**3
        fracsigrho[i] = unp.nominal_values(urho) / unp.std_devs(urho)

    return fracsigrho
