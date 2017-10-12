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


global G, rhoEarth
G, rhoEarth = 6.67e-11, 5.51


def estimate_Nrv_TESS(planetindex, band_strs, R, aperture_m, QE,
                      Z=0, sigmaRV_activity=0., sigmaRV_noisefloor=.5,
                      protseed=None, verbose=True):
    '''
    Estimate the number of RVs required to measure the mass of a transiting 
    TESS planet at a particular detection significance with a particular 
    instrument.

    Parameters
    ----------
    `planetindex': scalar
        Index of the TESS planet from the results of the predicted TESS planets 
        from Sullivan et al 2015
    `band_strs': list of strs
        A list of the spectral bands that span the wavelength coverage of the 
        spectrograph used to measure the radial velocities of the TESS star. 
        All band_strs entries must be in ['u','g','r','i','z','Y','J','H','K'] 
        and at least one of ['i','J','K'] must be included for scaling of the 
        TESS stars from Sullivan
    `R': scalar
        The spectral resolution of the spectrograph (lambda / d_lambda)
    `aperture_m': float
        The telescope's aperture diameter in meters
    `QE': scalar
        The quantum efficiency of the detector (0<QE<=1)
    `Z': scalar
        The assumed metallicity ([Fe/H]) of the TESS star in solar units
    `sigmaRV_activity': scalar
        An additive source of RV uncertainty from RV activity or jitter in m/s. 
        To be added in quadrature to the photon-noise RV precision derived for 
        the TESS star
    `protseed': scalar
        Seed for the random number generator used to draw the stellar rotation 
        period which is not know a-priori for the TESS stars from Sullivan
    `verbose': boolean
        If True, results for this star are printed to screen

    Returns
    -------
    `Nrv': int
        The number of radial velocity measurements required to detect the TESS 
        planet's mass at a given significance
    `texp': float
        The exposure time in minutes for each radial velocity measurement
    `sigmaRV_phot': float
        The photon-noise limit on the RV measurement precision in m/s
    `sigmaRV_eff': float
        The effective RV measurement precision from the combined effects of 
        photon-noise, instrument stability, multiple planets, and stellar 
        activity in m/s

    '''
    # Read-in TESS data for this planetary system
    ra,dec,rp,P,S,K,Rs,Teff,Vmag,Imag,Jmag,Kmag,dist,_,_,snr,_ = get_TESS_data()
    nplanets, planetindex = ra.size, int(planetindex)
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

    # Compute parameters of interest
    mp = get_planet_mass(rp)
    Ms = get_stellar_mass(P, mp, K)
    logg = np.log10(G * rvs.Msun2kg(Ms) / rvs.Rsun2m(Rs)**2 * 1e2)

    # Round Teff and logg
    Teffs = np.append(np.arange(23e2,7e3,1e2), np.arange(7e3,121e2,2e2))
    Teff_round = Teffs[abs(Teffs-Teff) == np.min(abs(Teffs-Teff))]
    loggs = np.arange(0,6.1,.5)
    logg_round = loggs[abs(loggs-logg) == np.min(abs(loggs-logg))]

    # Compute vsini
    ## TEMP
    Prot = _get_prot(Teff, seed=protseed)
    I = abs(np.arccos(np.random.uniform(-1,1)))
    vsini = 2*np.pi * rvs.Rsun2m(Rs)*1e-3 * np.sin(I) / rvs.days2sec(Prot)

    # optical or nIR bands?
    optical = True if np.any([band.islower() for band in band_strs]) else False
    nIR = not optical

    # Get the stellar magnitudes in the desired bands scaled to the results
    # from Sullivan
    known_mags = [Vmag, Imag, Jmag, Kmag]
    mags = _get_magnitudes(band_strs, known_mags, Teff_round, logg_round, Z,
                           optical=optical, nIR=nIR)

    # Estimate Nrv for this TESS planet
    startheta = mags, float(Teff_round), float(logg_round), Z, vsini
    planettheta = rp, mp, K
    instrumenttheta = band_strs, R, aperture_m, QE
    texp, sigmaRV_phot, sigmaRV_eff, Nrv = _estimate_Nrv(startheta, planettheta,
                                                         instrumenttheta,
                                        sigmaRV_activity=sigmaRV_activity,
                                        sigmaRV_noisefloor=sigmaRV_noisefloor)

    if verbose:
        print '\n%35s = %.3f minutes'%('Exposure time', texp)
        print '%35s = %.3f m/s'%('Photon-noise limited RV uncertainty',
                                 sigmaRV_phot)
        print '%35s = %.3f m/s'%('Effective RV uncertainty', sigmaRV_eff)
        print '%35s = %i'%('Number of RVs', Nrv)
    
    return Nrv, texp, sigmaRV_phot, sigmaRV_eff



def _estimate_Nrv(startheta, planettheta, instrumenttheta,
                  sigmaRV_activity=0., sigmaRV_noisefloor=.5):
    '''
    Estimate the number of RVs required to measure the semiamplitude K at a 
    given signficance of a particular planet around a particular star.
    
    Parameters
    ----------
    `startheta': tuple (5 entries)
        Containing the
            - list of stellar magnitudes
            - effective temperature (in K) included in the PHOENIX models
            - logg (in cgs units) included in the PHOENIX models
            - metallicity ([Fe/H] in solar units) included in the PHOENIX models
            - projected rotation velocity (in km/s)
    `planettheta': tuple (3 entries)
        Containing the
            - planet radius (in Earth radii)
            - planet mass (in Earth masses)
            - RV semi-amplitude (in m/s)
    `instrumenttheta': tuple (4 entries)
        Containing the
            - list of spectral bands
            - spectral resolution (lambda / d_lambda)
            - aperture diameter (in m)
            - quantum efficiency (fraction of 1)
    `sigmaRV_activity': scalar
        An additive source of RV uncertainty from RV activity or jitter in m/s. 
        To be added in quadrature to the photon-noise RV precision derived for 
        the TESS star

    Returns
    -------
    `texp': float
        The exposure time in minutes for each radial velocity measurement
    `sigmaRV_phot': float
        The photon-noise limit on the RV measurement precision in m/s
    `sigmaRV_eff': float
        The effective RV measurement precision from the combined effects of 
        photon-noise, instrument stability, multiple planets, and stellar 
        activity in m/s
    `Nrv': int
        The number of radial velocity measurements required to detect the TESS 
        planet's mass at a given significance

    '''
    mags, Teff_round, logg_round, Z, vsini = startheta
    rp, mp, K = planettheta
    band_strs, R, aperture_m, QE = instrumenttheta
    
    # compute texp in each band
    texps = np.zeros(len(mags))
    for i in range(texps.size):
        wl, spec = get_reduced_spectrum(Teff_round, logg_round, Z, vsini,
                                        band_strs[i], R)
        texps[i] = exposure_time_calculator_per_band(mags[i], band_strs[i],
                                                     aperture_m, QE, R)
    texp = np.median(texps)

    # compute sigmaRV in each band for a fixed texp
    sigmaRVs = np.zeros(texps.size)
    for i in range(sigmaRVs.size):
        wl, spec = get_reduced_spectrum(Teff_round, logg_round, Z, vsini,
                                        band_strs[i], R)
        sigmaRVs[i] = compute_sigmaRV(wl, spec, mags[i], band_strs[i], texp,
                                      aperture_m, QE, R)

    # Compute the effective sigmaRV
    sigmaRV_phot = 1. / np.sqrt(np.sum(1./sigmaRVs**2))
    sigmaRV_phot = sigmaRV_phot if sigmaRV_phot > sigmaRV_noisefloor \
                   else float(sigmaRV_noisefloor)
    sigmaRV_eff = np.sqrt(sigmaRV_phot**2 + sigmaRV_activity**2)

    # Compute Nrv to measure K at a given significance
    ##sigmaK_target = get_sigmaK_target_v1(mp, rp)
    sigmaK_target = get_sigmaK_target_v2(K)
    Nrv = int(np.round(2 * (sigmaRV_eff / sigmaK_target)**2))

    # Return exposure time in minutes, photon-noise limit on sigmaRV,
    # the effective sigmaRV (plus activity and planets), and the number of RVs
    return texp, sigmaRV_phot, sigmaRV_eff, Nrv


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
        if 'v' in band_strs:
            ref_band, ref_mag = 'v', Vmag
        elif 'i' in band_strs:
            ref_band, ref_mag = 'i', Imag
        else:
            raise ValueError('Do not have an optical reference magnitude ' + \
                             "in `band_strs'. Must include wither v or i.")

    else:  # nIR bands
        if 'J' in band_strs:
            ref_band, ref_mag = 'J', Jmag
        elif 'K' in band_strs:
            ref_band, ref_mag = 'K', Kmag
        else:
            raise ValueError('Do not have a nIR reference magnitude ' + \
                             "in `band_strs'. Must include either J or K.")
        
    # Integrate the spectrum over each band of interest to get flux
    fluxes = np.zeros(len(band_strs))
    for i in range(fluxes.size):

        # Get total flux over the bandpass
        wlmin, wlmax,_ = get_band_range(band_strs[i])
        g = (wl >= wlmin) & (wl <= wlmax)
        wl2, spectrum2 = wl[g], spectrum[g]
        fluxes[i] = np.sum(spectrum2)

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


def get_sigmaK_target_v2(K):
    return .13 * K


# GJ1132: mags=array([ 16.44372851,  13.53664024,  13.04561156,  12.3])
def TEST_estimate_Nrv_TESS(mags=[16.4, 13.5, 13.0, 12.3],
                           band_strs=['u','g','r','i'],
                           Teff_round=3300, logg_round=5, vsini=.01, rp=1.1,
                           mp=1.6, K=2.8, R=1e5, aperture_m=3.6, QE=.1, Z=0,
                           sigmaRV_activity=0., protseed=None):
    startheta = mags, float(Teff_round), float(logg_round), Z, vsini
    planettheta = rp, mp, K
    instrumenttheta = band_strs, R, aperture_m, QE
    texp, sigmaRV, sigmaRV_eff, Nrv = _estimate_Nrv(startheta, planettheta,
                                                    instrumenttheta,
                                            sigmaRV_activity=sigmaRV_activity)
    print '\nExposure time = %.3f min'%texp
    print 'Photon-noise limited RV uncertainty = %.3f m/s'%sigmaRV
    print 'Effective RV uncertainty = %.3f m/s'%sigmaRV_eff
    print 'Number of RVs = %i'%Nrv


# GJ 1214:
# mags=np.array([ 16.3387182 ,  13.31858571,  12.71599573,  11.5       ])
# Teff_round, rp, mp, K = 3000, 2.7, 6.55, 12.2
