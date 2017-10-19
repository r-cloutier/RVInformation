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


global G, rhoEarth, c, h
G, rhoEarth, c, h = 6.67e-11, 5.51, 299792458., 6.62607004e-34


def estimate_Nrv_TESS(planetindex, band_strs, R, aperture_m, QE,
                      Z=0, sigmaRV_activity=0., sigmaRV_noisefloor=.5,
                      protseed=None, testplanet_sigmaK=0, verbose=True):
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
        All band_strs entries must be in 
        ['U','B','V','R','I','Z','Y','J','H','K']
        and at least one of ['V','I','J','K'] must be included for scaling of 
        the TESS stars from Sullivan
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
    `testplanet_sigmaK': scalar
        If 0, assume we are calculating a TESS planet and use its mass to 
        determine the required constraint on the K measurement uncertainty 
        (sigmaK). Otherwise, set this value to the fractional K measurement 
        uncertainty as is done for test cases (e.g. 0.33 for GJ1132)  
    `verbose': boolean
        If True, results for this star are printed to screen

    Returns
    -------
    `Nrv': int
        The number of radial velocity measurements required to detect the TESS 
        planet's mass at a given significance
    `texp': float
        The exposure time in minutes for each radial velocity measurement
    `tobserving': float
        The total observing time in hours including a constant estimate of 
        overheads
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

    # Get the stellar magnitudes in the desired bands scaled to the results
    # from Sullivan
    known_mags = [Vmag, Imag, Jmag, Kmag]
    mags = _get_magnitudes(band_strs, known_mags, Teff_round, logg_round, Z, Ms)
    
    # Estimate Nrv for this TESS planet
    startheta = mags, float(Teff_round), float(logg_round), Z, vsini
    planettheta = rp, mp, K
    instrumenttheta = band_strs, R, aperture_m, QE
    Nrv, texp, tobserving, sigmaRV_phot, sigmaRV_eff = \
                        estimate_Nrv(startheta, planettheta, instrumenttheta,
                                     sigmaRV_activity=sigmaRV_activity,
                                     sigmaRV_noisefloor=sigmaRV_noisefloor,
                                     testplanet_sigmaK=testplanet_sigmaK)

    if verbose:
        print '\n%35s = %.3f m/s'%('Photon-noise limited RV uncertainty',
                                   sigmaRV_phot)
        print '%35s = %.3f m/s'%('Effective RV uncertainty', sigmaRV_eff)
        print '%35s = %i'%('Number of RVs', Nrv)
        print '%35s = %.3f minutes'%('Exposure time', texp)
        print '%35s = %.3f hours'%('Total observing time', tobserving)
    
    return Nrv, texp, tobserving, sigmaRV_phot, sigmaRV_eff


def estimate_Nrv(startheta, planettheta, instrumenttheta,
                 sigmaRV_activity=0., sigmaRV_noisefloor=.5, texpmin=10,
                 texpmax=60, testplanet_sigmaK=0):
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
    `texpmin': scalar
        The minimum exposure time in minutes. Required to mitigate the effects 
        of stellar pulsations and granulation
    `texpmax': scalar
        The maximum exposure time in minutes. Used to moderate the limit the 
        observational time that can be dedicated to a single star
    `testplanet_sigmaK': scalar
        If 0, assume we are calculating a TESS planet and use its mass to 
        determine the required constraint on the K measurement uncertainty 
        (sigmaK). Otherwise, set this value to the fractional K measurement 
        uncertainty as is done for test cases (e.g. 0.33 for GJ1132)  
    
    Returns
    -------
    `Nrv': int
        The number of radial velocity measurements required to detect the TESS 
        planet's mass at a given significance
    `texp': float
        The exposure time in minutes for each radial velocity measurement
    `tobserving': float
        The total observing time in hours including a constant estimate of 
        overheads
    `sigmaRV_phot': float
        The photon-noise limit on the RV measurement precision in m/s
    `sigmaRV_eff': float
        The effective RV measurement precision from the combined effects of 
        photon-noise, instrument stability, multiple planets, and stellar 
        activity in m/s

    '''
    mags, Teff_round, logg_round, Z, vsini = startheta
    rp, mp, K = planettheta
    band_strs, R, aperture_m, QE = instrumenttheta
    mags, band_strs = np.ascontiguousarray(mags), \
                      np.ascontiguousarray(band_strs)
    assert mags.size == band_strs.size
    
    # compute texp in a reference band (either V or J)
    texp = exposure_time_calculator_per_band(mags, band_strs, aperture_m, QE,
                                             R, texpmin=texpmin,
                                             texpmax=texpmax)

    # compute sigmaRV in each band for a fixed texp
    sigmaRVs = np.zeros(len(mags))
    for i in range(sigmaRVs.size):
        wl, spec = get_reduced_spectrum(Teff_round, logg_round, Z, vsini,
                                        band_strs[i], R, pltt=False)
        sigmaRVs[i] = compute_sigmaRV(wl, spec, mags[i], band_strs[i], texp,
                                      aperture_m, QE, R)

    # Compute the effective sigmaRV
    sigmaRV_phot = 1. / np.sqrt(np.sum(1./sigmaRVs**2))
    sigmaRV_phot = sigmaRV_phot if sigmaRV_phot > sigmaRV_noisefloor \
                   else float(sigmaRV_noisefloor)
    sigmaRV_eff = np.sqrt(sigmaRV_phot**2 + sigmaRV_activity**2)

    # Compute Nrv to measure K at a given significance
    if testplanet_sigmaK != 0:     # use for testing
        sigmaK_target = get_sigmaK_target_v2(K, testplanet_sigmaK)
    else:                          # use for TESS planets
        sigmaK_target = get_sigmaK_target_v1(mp, rp)
    Nrv = int(np.round(2 * (sigmaRV_eff / sigmaK_target)**2))

    toverhead = 5.
    tobserving = (texp+toverhead)*Nrv / 6e1
    
    return Nrv, texp, tobserving, sigmaRV_phot, sigmaRV_eff



def _get_prot(Teff, seed=None):
    '''
    Draw a stellar rotation period based on the measured distribution from
    McQuillan+2014 (2014ApJS..211...24M).

    Parameters
    ----------
    `Teff': scalar
        The effective temperature of the star whose rotation period is being 
        sampling
    `protseed': scalar
        Seed for the random number generator used to draw the stellar rotation 
        period which is not know a-priori for the TESS stars from Sullivan
    
    Returns
    -------
    `protseed': float
        The star's sampled rotation period in days

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


 
def _get_magnitudes(band_strs, known_mags, Teff, logg, Z, Ms):
    '''
    Get the stellar apparent magnitude in all bands of interest and normalized 
    by the known magnitudes from Sullivan et al 2015.

    Parameters
    ----------
    `band_strs': list of strs
        A list of the spectral bands that span the wavelength coverage of the 
        spectrograph used to measure the radial velocities of the TESS star. 
        All band_strs entries must be in 
        ['U','B','V','R','I','Z','Y','J','H','K']
        and at least one of ['V','I','J','K'] must be included for scaling of 
        the TESS stars from Sullivan
    `known_mags': list of scalars
        A list of the V, I, J, and K stellar apparent magnitudes from Sullivan
    `Teff': scalar
        The stellar effective temperature in Kelvin
    `logg': scalar
        The stellar logg in cgs units
    `Z': scalar
        The stellar metallicity [Fe/H] in solar units
    
    Returns
    -------
    `mags': numpy.array
        Array of the stellar magnitudes for the star based on its spectrum and 
        scaled to the known magnitudes from Sullivan

    '''
    # Use isochrone colours to compute mags in each band of interest
    # solar metallicity at a fixed age of 10^9 yrs
    Mu,Mb,Mv,Mr,Mi,Mj,Mh,Mk = _get_absolute_stellar_magnitudes(Ms)

    Vmag, Imag, Jmag, Kmag = known_mags
    if 'V' in band_strs:
        ref_band, ref_mag, ref_absmag = 'V', Vmag, Mv
    elif 'J' in band_strs:
        ref_band, ref_mag, ref_absmag = 'J', Jmag, Mj
    elif 'I' in band_strs:
        ref_band, ref_mag, ref_absmag = 'I', Imag, Mi
    elif 'K' in band_strs:
        ref_band, ref_mag, ref_absmag = 'K', Kmag, Mk
    else:
        raise ValueError('Do not have a reference magnitude ' + \
                         "in `band_strs'. Must include one of V, I, J or K.")

    mags = np.zeros(len(band_strs))
    for i in range(mags.size):
        if band_strs[i] == 'U':
            absmag = Mu
        elif band_strs[i] == 'B':
            absmag = Mb
        elif band_strs[i] == 'V':
            absmag = Mv
        elif band_strs[i] == 'R':
            absmag = Mr
        elif band_strs[i] == 'I':
            absmag = Mi
        elif band_strs[i] == 'J':
            absmag = Mj
        elif band_strs[i] == 'H':
            absmag = Mh
        elif band_strs[i] == 'K':
            absmag = Mk
        else:
            raise ValueError('Unknown passband: %s'%band_strs[i])
        mags[i] = ref_mag - ref_absmag + absmag

    return mags


def _get_absolute_stellar_magnitudes(Ms, logage=9):
    '''
    Get the absolute magnitudes of a star with a given stellar mass at a 
    given age using the isochrones from 2005A&A...436..895G

    Parameters
    ----------
    `Ms': scalar
        The stellar mass in MSun
    `logage': scalar
        The log base 10 of the age of the star in years

    Returns
    -------
    `mags': numpy.array
        The absolute magnitudes of the star 

    '''
    logages,Mss,Mus,Mbs,Mvs,Mrs,Mis,Mjs,Mhs,Mks = \
                                np.loadtxt('input_data/isoc_z019.dat',
                                usecols=(0,1,7,8,9,10,11,12,13,14)).T
    g = logages == 9
    Mss,Mus,Mbs,Mvs,Mrs,Mis,Mjs,Mhs,Mks = Mss[g],Mus[g],Mbs[g],Mvs[g],Mrs[g], \
                                          Mis[g],Mjs[g],Mhs[g],Mks[g]
    g = abs(Mss-Ms) == np.min(abs(Mss-Ms))
    Mu,Mb,Mv,Mr,Mi,Mj,Mh,Mk = Mus[g],Mbs[g],Mvs[g],Mrs[g],Mis[g],Mjs[g], \
                              Mhs[g],Mks[g]
    return Mu, Mb, Mv, Mr, Mi, Mj, Mh, Mk 
        

def get_planet_mass(rp):
    '''
    Compute the TESS planet mass from its reported radius using the 
    deterministic mass-radius relation from Weiss & Marcy 2014.

    Parameters
    ----------
    `rp': scalar
        The planet radius in Earth radii to be converted into a mass

    Returns
    -------
    `mp': scalar
        The planet mass in Earth masses

    '''
    if rp < 1.5:
        return .44*rp**3 + .614*rp**4
    else:
        return 2.69*rp**.93


def get_stellar_mass(P, mp, K):
    '''
    Compute the stellar mass from the orbital period, planet mass, and RV 
    semi-amplitude provided in Sullivan et al 2015.

    Parameters
    ----------
    `P': scalar
        The orbital period (in days) of the planet orbiting the star whose mass 
        is being measured
    `mp': scalar
        The planet mass (in Earth masses) of the planet orbiting the star 
        whose mass is being measured
    `K': scalar
        The RV semi-amplitude (in m/s) of the planet orbiting the star whose 
        mass is being measured

    Returns
    -------
    `Ms': float
        The stellar mass in solar masses
    
    '''
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


def get_sigmaK_target_v2(K, fracsigmaK):
    print 'fracsigmaK = %.3f'%fracsigmaK
    return fracsigmaK * K
