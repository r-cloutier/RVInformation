'''
Compute nRV for all TESS targets using the code cloned from
https://github.com/r-cloutier/Transiting_RVcalculator.git
'''
import numpy as np
from get_tess_data import get_TESS_data
from compute_sigmaRV import *
import pylab as plt
import rvs
from uncertainties import unumpy as unp


global G, rhoEarth, R, aperture, QE
G, rhoEarth, R, aperture, QE = 6.67e-11, 5.51, 75e3, 3.58, .1


def estimate_Nrv_TESS(planetindex, band_strs, R, aperture_m, QE, Z=0,
                      optical=False, nIR=True, protseed=None):
    '''
    Estimate the number of RVs to measure the mass of a transiting planet at 
    a given detection significance with a particular instrument.
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

    # get vsini
    Prot = _get_prot(Teff, seed=protseed)
    I = np.arccos(np.random.uniform(-1,1))
    vsini = 2*np.pi * rvs.Rsun2m(Rs)*1e-3 * np.sin(I) / rvs.days2sec(Prot)
    
    # get optical and nIR magnitudes from those given in the Sullivan sample
    known_mags = [Vmag, Imag, Jmag, Kmag]
    mags = _get_magnitudes(band_strs, known_mags, Teff_round, logg_round, Z,
                           optical=optical, nIR=nIR)

    # Estimate Nrv for this TESS planet
    startheta = mags, float(Teff_round), float(logg_round), Z, vsini
    planettheta = rp, mp, K
    instrumenttheta = band_strs, R, aperture_m, QE
    Nrv = _estimate_Nrv(startheta, planettheta, instrumenttheta)

    return Nrv


def _estimate_Nrv(startheta, planettheta, instrumenttheta):
    '''
    Estimate the number of RVs required to measure the semiamplitude K at a 
    given signficance of a particular planet around a particular star.
    '''
    mags, Teff_round, logg_round, Z, vsini = startheta
    rp, mp, K = planettheta
    band_strs, R, aperture_m, QE = instrumenttheta
    
    # compute sigmaRV in each band
    sigmaRVs = np.zeros(len(mags))
    for i in range(sigmaRVs.size):
        wl, spec = get_reduced_spectrum(Teff_round, logg_round, Z, vsini,
                                        band_strs[i], R)
        texp = exposure_time_calculator_per_band(mags[i], band_strs[i],
                                                 aperture_m, QE, R)
        sigmaRVs[i] = compute_sigmaRV(wl, spec, mags[i], band_strs[i], texp,
                                      aperture_m, QE, R)
        print band_strs[i], texp, sigmaRVs[i]
        
    # Compute the effective sigmaRV
    sigmaRV = 1. / np.sqrt(np.sum(1./sigmaRVs**2))
    sigmaRV_activity = 0.
    sigmaRV_eff = np.sqrt(sigmaRV**2 + sigmaRV_activity**2)

    # Compute Nrv to measure K at a given significance
    ##sigmaK_target = get_sigmaK_target_v1(mp, rp)
    sigmaK_target = .2 * K
    Nrv = np.round(2 * (sigmaRV_eff / sigmaK_target)**2)

    return Nrv


def _get_prot(Teff, seed=None):
    '''
    Draw a stellar rotation period based on the measured distribution from
    McQuillan+2014 (2014ApJS..211...24M).
    '''
    Teffs, Prots = np.loadtxt('input_data/asu.tsv', skiprows=37).T
    # Isolate range of effective temperatures
    dT = 1e2
    if Teff > Teffs.max():
        g = Teffs >= Teffs.max()-dT
    elif Teff < Teffs.min():
        g = Teffs <= Teffs.min()+dT
    else:
        g = (Teffs >= Teff-dT) & (Teffs <= Teff+dT)
    # Set seed
    if seed != None:
        np.random.seed(int(seed))
    return np.random.choice(Prots[g]) + np.random.randn() * .1

 
def _get_magnitudes(band_strs, known_mags, Teff, logg, Z,
                    optical=False, nIR=False):
    '''
    Get the stellar apparent magnitude in all bands of interest and normalized 
    by the known magnitudes from Sullivan et al 2015.
    '''
    # Get the full spectrum
    wl = get_wavelengthgrid()
    _, spectrum = get_full_spectrum(Teff, logg, Z)

    # Get reference magnitude
    if optical:
        nIR = False
    elif not nIR:
        raise ValueError("One of `optical' or `nIR' keyword arguments " + \
                         "must be set to True.")
    Vmag, Imag, Jmag, Kmag = known_mags
    if optical:  # optical bands
        if 'V' in band_strs:
            ref_band, ref_mag = 'V', Vmag
        elif 'i' in band_strs:
            ref_band, ref_mag = 'i', Imag
        else:
            raise ValueError('Do not have an optical reference magnitude ' + \
                             "in `band_strs'.")

    else:  # nIR bands
        if 'J' in band_strs:
            ref_band, ref_mag = 'J', Jmag
        elif 'K' in band_strs:
            ref_band, ref_mag = 'K', Kmag
        else:
            raise ValueError('Do not have a nIR reference magnitude ' + \
                             "in `band_strs'.")
        
    # Integrate the spectrum over each band of interest to get flux
    fluxes = np.zeros(len(band_strs))
    for i in range(fluxes.size):
        # Get band transmission
        #wl_band, transmission,_ = get_band_transmission(band_strs[i])
	wlmin, wlmax,_ = get_band_range(band_strs[i])
	wl_band = np.linspace(wlmin, wlmax, 10)
	transmission = np.ones(wl_band.size)

        # Get total flux over the bandpass
        fint = interp1d(wl_band, transmission)
        g = (wl >= wl_band.min()) & (wl <= wl_band.max())
        wl2, spectrum2 = wl[g], spectrum[g]
        transmission2 = fint(wl2)
        fluxes[i] = np.sum(spectrum2 * transmission2)

        # Get reference flux
        if band_strs[i] == ref_band:
            ref_flux = fluxes[i]

    # Convert to magnitudes
    mag  = -2.5*np.log10(fluxes / ref_flux) + ref_mag    
    return mag
        

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

    return None
