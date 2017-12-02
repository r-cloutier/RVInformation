'''
Compute nRV for all TESS targets using the code cloned from
https://github.com/r-cloutier/Transiting_RVcalculator.git
'''
import numpy as np
from get_tess_data import get_TESS_data
from compute_sigmaRV import *
from sigmaRV_activity import *
from sigmaRV_planets import *
import pylab as plt
import glob, rvs, sys, os
from uncertainties import unumpy as unp


global G, rhoEarth, c, h
G, rhoEarth, c, h = 6.67e-11, 5.51, 299792458., 6.62607004e-34


def estimate_Nrv_TESS(planetindex, band_strs, R, aperture_m,
                      QE=.1, Z=0, sigmaRV_activity=0, sigmaRV_planets=0,
                      sigmaRV_noisefloor=.5, testingseed=False,
                      testplanet_sigmaKfrac=0, systnum=0, verbose=True):
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
        ['U','B','V','R','I','Y','J','H','K']
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
    `sigmaRV_planets': scalar
        An additive source of RV uncertainty from unseen planets in m/s. 
        To be added in quadrature to the photon-noise RV precision derived for 
        the TESS star
    `testingseed': boolean
        If True, use a seed for the random number generator used to draw the stellar 
        rotation period, RV activity, and RV from planets. Useful for testing
    `testplanet_sigmaKfrac': scalar
        If 0, assume we are calculating a TESS planet and use its mass to 
        determine the required constraint on the K measurement uncertainty 
        (sigmaK). Otherwise, set this value to the fractional K measurement 
        uncertainty as is done for test cases (e.g. 0.33 for GJ1132)
    `systnum': scalar
        System index used when saving the simulation results
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
    ra,dec,rp,P,S,K,Rs,Teff,Vmag,Imag,Jmag,Kmag,dist,_,_,_,mult = \
                                                            get_TESS_data()
    nplanets, planetindex = ra.size, int(planetindex)
    assert 0 <= planetindex < nplanets
    ra, dec = ra[planetindex], dec[planetindex]
    rp,P,S,K,Rs,Teff,Vmag,Imag,Jmag,Kmag,dist,mult = rp[planetindex], \
                                                     P[planetindex], \
                                                     S[planetindex], \
                                                     K[planetindex], \
                                                     Rs[planetindex], \
                                                     Teff[planetindex], \
                                                     Vmag[planetindex], \
                                                     Imag[planetindex], \
                                                     Jmag[planetindex], \
                                                     Kmag[planetindex], \
                                                     dist[planetindex], \
                                                     mult[planetindex]
    mult = 1. if mult == 0 else mult
    
    # Compute parameters of interest
    mp = get_planet_mass(rp)
    Ms = get_stellar_mass(P, mp, K)
    logg = np.log10(G * rvs.Msun2kg(Ms) / rvs.Rsun2m(Rs)**2 * 1e2)

    # Round Teff and logg
    Teffs = np.append(np.arange(23e2,7e3,1e2), np.arange(7e3,121e2,2e2))
    Teff_round = Teffs[abs(Teffs-Teff) == np.min(abs(Teffs-Teff))][0]
    loggs = np.arange(0,6.1,.5)
    logg_round = loggs[abs(loggs-logg) == np.min(abs(loggs-logg))][0]

    # Get the stellar magnitudes in the desired bands scaled to the results
    # from Sullivan
    known_mags = [Vmag, Imag, Jmag, Kmag]
    band_strs_tmp = list(np.append(band_strs, 'B'))
    mags = _get_magnitudes(band_strs_tmp, known_mags, Teff_round, logg_round,
                           Z, Ms)
    mags, Bmag = mags[:-1], float(mags[-1])
    B_V = Bmag - Vmag
    
    # compute vsini
    if testingseed:
        np.random.seed(1)
    fname = 'Results/star%.4d/TESSplanet%.4d_%s_%.4d.dat'%(planetindex,
                                                           planetindex,
                                                           ''.join(band_strs),
                                                            systnum)
    Prot, vsini = _draw_prot(Ms, Rs, fname)

    # Estimate Nrv for this TESS planet
    startheta = mags, float(Teff_round), float(logg_round), Z, vsini, Ms, \
                Prot, B_V
    planettheta = rp, mp, K, P, mult
    instrumenttheta = band_strs, R, aperture_m, QE
    Nrv,texp,tobs,sigK_target,sig_phot,sig_act,sig_planets,sig_eff = \
                            estimate_Nrv(startheta, planettheta,
                                         instrumenttheta, fname=fname,
                                         sigmaRV_activity=sigmaRV_activity,
                                         sigmaRV_planets=sigmaRV_planets,
                                         sigmaRV_noisefloor=sigmaRV_noisefloor,
                                    testplanet_sigmaKfrac=testplanet_sigmaKfrac)

    if verbose:
        print '\n%35s = %.3f m/s'%('Photon-noise limited RV uncertainty',
                                   sigmaRV_phot)
        print '%35s = %.3f m/s'%('Effective RV uncertainty', sig_eff)
        print '%35s = %.3f'%('Target fractional K uncertainty', sigK_target/K)
        print '%35s = %i'%('Number of RVs', Nrv)
        print '%35s = %.3f minutes'%('Exposure time', texp)
        print '%35s = %.3f hours'%('Total observing time', tobs)

    # Save values
    save_results(planetindex, band_strs, mags, ra, dec, P, rp, mp, K, S, Ms,
                 Rs, Teff, dist, Prot, vsini, Z, sig_act,
                 sig_planets, R, aperture_m, QE, sigK_target/K,
                 sig_phot, sig_eff, texp, tobs, Nrv, fname)


def estimate_Nrv(startheta, planettheta, instrumenttheta,
                 fname='', sigmaRV_activity=0., sigmaRV_planets=0.,
                 sigmaRV_noisefloor=.5, texpmin=10,
                 texpmax=60, testplanet_sigmaKfrac=0):
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
            - stellar mass (in solar masses)
            - stellar rotation period (in days) 
            - B-V colour 
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
    `sigmaRV_activity': scalar
        An additive source of RV uncertainty from unseen planets in m/s. 
        To be added in quadrature to the photon-noise RV precision derived for 
        the TESS star
    `texpmin': scalar
        The minimum exposure time in minutes. Required to mitigate the effects 
        of stellar pulsations and granulation
    `texpmax': scalar
        The maximum exposure time in minutes. Used to moderate the limit the 
        observational time that can be dedicated to a single star
    `testplanet_sigmaKfrac': scalar
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
    mags, Teff_round, logg_round, Z, vsini, Ms, Prot, B_V = startheta
    rp, mp, K, P, mult = planettheta
    band_strs, R, aperture_m, QE = instrumenttheta
    mags, band_strs = np.ascontiguousarray(mags), \
                      np.ascontiguousarray(band_strs)
    assert mags.size == band_strs.size
    
    # compute texp in a reference band (either V or J)
    texp = exposure_time_calculator_per_band(mags, band_strs, aperture_m, QE, R,
                                             texpmin=texpmin, texpmax=texpmax)

    # compute sigmaRV in each band for a fixed texp
    sigmaRVs = np.zeros(mags.size)
    for i in range(sigmaRVs.size):
        wl, spec = get_reduced_spectrum(Teff_round, logg_round, Z, vsini,
                                        band_strs[i], R, pltt=False)
        sigmaRVs[i] = compute_sigmaRV(wl, spec, mags[i], band_strs[i], texp,
                                      aperture_m, QE, R)

    # Apply corrections from Artigau+2017 to bands with a known correction
    correctionsYJHK = np.array([.47, .63, 1.59, 1.72])
    sigmaRVs[band_strs == 'Y'] = sigmaRVs[band_strs == 'Y'] / correctionsYJHK[0]
    sigmaRVs[band_strs == 'J'] = sigmaRVs[band_strs == 'J'] / correctionsYJHK[1]
    sigmaRVs[band_strs == 'H'] = sigmaRVs[band_strs == 'H'] / correctionsYJHK[2]
    sigmaRVs[band_strs == 'K'] = sigmaRVs[band_strs == 'K'] / correctionsYJHK[3]

    # Compute sigmaRV over all bands
    sigmaRV_phot = 1. / np.sqrt(np.sum(1./sigmaRVs**2))
    sigmaRV_phot = sigmaRV_phot if sigmaRV_phot > sigmaRV_noisefloor \
                   else float(sigmaRV_noisefloor)

    # estimate sigmaRV due to stellar activity
    fs = np.array(glob.glob('%s_*_%s'%(fname.split('_')[0],
                                       fname.split('_')[-1])))
    if fs.size > 0:
        sigmaRV_activity = np.loadtxt(fs[0])[14]
    else:
        if sigmaRV_activity == 0 and testplanet_sigmaKfrac == 0:
            sigmaRV_activity = get_sigmaRV_activity(Teff_round, Ms, Prot, B_V)
        
    # estimate sigmaRV due to unseen planets
    if fs.size > 0:
        sigmaRV_planets = np.loadtxt(fs[0])[15]
    else:
        if sigmaRV_planets == 0 and testplanet_sigmaKfrac == 0 and mult > 1:
            sigmaRV_planets = get_sigmaRV_planets(P, rp, Teff_round, Ms, mult,
                                                  sigmaRV_phot)
            
    # compute effective sigmaRV
    sigmaRV_eff = np.sqrt(sigmaRV_phot**2 + \
                          sigmaRV_activity**2 + \
                          sigmaRV_planets**2)

    # Compute Nrv to measure K at a given significance
    if testplanet_sigmaKfrac != 0:     # use for testing
        sigmaK_target = get_sigmaK_target_v2(K, testplanet_sigmaKfrac)
    else:                          # use for TESS planets
        sigmaK_target = get_sigmaK_target_v3(P, Ms, K)
    Nrv = 2 * (sigmaRV_eff / sigmaK_target)**2

    toverhead = 5.
    tobserving = (texp+toverhead)*Nrv / 6e1
    
    return Nrv, texp, tobserving, sigmaK_target, sigmaRV_phot, sigmaRV_activity, sigmaRV_planets, sigmaRV_eff

 
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
        ['U','B','V','R','I','Y','J','H','K']
        and at least one of ['V','I','J','K'] must be included for scaling of 
        the TESS stars from Sullivan
    `known_mags': list of scalars
        A list of the v, i, J, and K stellar apparent magnitudes from Sullivan
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
    Mu, Mb, Mv, Mr, Mi, MY, MJ, MH, MK = _get_absolute_stellar_magnitudes(Ms)

    Vmag, Imag, Jmag, Kmag = known_mags
    if 'V' in band_strs:
        ref_band, ref_mag, ref_absmag = 'V', Vmag, Mv
    elif 'J' in band_strs:
        ref_band, ref_mag, ref_absmag = 'J', Jmag, MJ
    elif 'I' in band_strs:
        ref_band, ref_mag, ref_absmag = 'I', Imag, Mi
    elif 'K' in band_strs:
        ref_band, ref_mag, ref_absmag = 'K', Kmag, MK
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
        elif band_strs[i] == 'Y':
            absmag = MY
        elif band_strs[i] == 'J':
            absmag = MJ
        elif band_strs[i] == 'H':
            absmag = MH
        elif band_strs[i] == 'K':
            absmag = MK
        else:
            raise ValueError('Unknown passband: %s'%band_strs[i])
        mags[i] = ref_mag - ref_absmag + absmag

    return mags


def _get_absolute_stellar_magnitudes(Ms):
    '''
    Get the absolute magnitudes of a star with a given stellar mass at a 
    given age using the isochrones from 2005A&A...436..895G

    Parameters
    ----------
    `Ms': scalar
        The stellar mass in MSun
    
    Returns
    -------
    `mags': numpy.array
        The absolute magnitudes of the star 

    '''
    # Appoximate MS lifetime to see at what age to obtain colours
    logtMS_yrs = np.log10(1e10 * (1./Ms)**(2.5))
    logage = round(logtMS_yrs*.77 / 5e-2) * 5e-2   # get ~3/4 through MS

    # First set of isochrones (ubvri)
    logages,Mss,Mus,Mbs,Mvs,Mrs,Mis,Mjs,Mhs,Mks = \
                                np.loadtxt('input_data/isoc_z019_ubvrijhk.dat',
                                usecols=(0,1,7,8,9,10,11,12,13,14)).T
    g = abs(logages-logage) == np.min(abs(logages-logage))
    Mss,Mus,Mbs,Mvs,Mrs,Mis,Mjs,Mhs,Mks = Mss[g],Mus[g],Mbs[g],Mvs[g],Mrs[g], \
                                          Mis[g],Mjs[g],Mhs[g],Mks[g]
    g = abs(Mss-Ms) == np.min(abs(Mss-Ms))
    if g.sum() > 1:
	g = np.where(g)[0][0]
    Mu,Mb,Mv,Mr,Mi,Mj,Mh,Mk = Mus[g],Mbs[g],Mvs[g],Mrs[g],Mis[g],Mjs[g], \
                              Mhs[g],Mks[g]
    # Second set of isochrones (ZYJHK)
    logages2,Mss2,MZs,MYs,MJs,MHs,MKs = \
                                np.loadtxt('input_data/isoc_z019_ZYJHK.dat',
                                usecols=(0,1,7,8,9,10,11)).T
    g = abs(logages2-logage) == np.min(abs(logages2-logage))
    Mss2,MZs,MYs,MJs,MHs,MKs = Mss2[g],MZs[g],MYs[g],MJs[g],MHs[g],MKs[g]
    g = abs(Mss2-Ms) == np.min(abs(Mss2-Ms))
    if g.sum() > 1:
	g = np.where(g)[0][0]
    MZ,MY,MJ,MH,MK = MZs[g],MYs[g],MJs[g],MHs[g],MKs[g]    

    return Mu, Mb, Mv, Mr, Mi, MY, MJ, MH, MK 


def _draw_prot(Ms, Rs, fname):
    '''
    Draw a rotation period and compute vsini or use the existing parameters 
    for this system.
    '''
    # If a simulation of this system exists, read in those parameters
    fs = np.array(glob.glob('%s_*_%s'%(fname.split('_')[0], fname.split('_')[-1])))
    if fs.size > 0:
        Prot, vsini = np.loadtxt(fs[0], usecols=(11,12))
        
    # Otherwise sample new rotation information 
    else:
        Prot = draw_prot_empirical(Ms)
        I = abs(np.arccos(np.random.uniform(-1,1)))
        vsini = 2*np.pi * rvs.Rsun2m(Rs)*1e-3 * np.sin(I) / rvs.days2sec(Prot) 
    return Prot, vsini
    

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


#def get_sigmaK_target_v1(rp, K, P, Ms, sigP=5e-5,
#                         fracsigrp=.05, fracsigMs=.1, fracsigrho_target=.3):
#    '''
#    Compute the K detection significance required to measure a planet's density 
#    at a specific fractional precision.
#    '''
#    urp, uP, uMs = unp.uarray(rp, rp*fracsigrp), unp.uarray(P, sigP), \
#                   unp.uarray(Ms, Ms*fracsigMs)
#
#    fracsigKs = np.logspace(-2,0,100)
#    fracsigrho = np.zeros(fracsigKs.size)
#    
#    for i in range(fracsigKs.size):
#        uK = unp.uarray(K, K*fracsigKs[i])
#        ump = rvs.RV_mp(uP, uMs, uK)
#        urho = rhoEarth * ump / urp**3
#        fracsigrho[i] =  unp.std_devs(urho) / unp.nominal_values(urho)
#
#    # find corresponding fractional sigmaK
#    fint = interp1d(fracsigrho, fracsigKs)
#    return float(fint(fracsigrho_target) * K)


def get_sigmaK_target_v2(K, fracsigmaK):
    return float(fracsigmaK * K)


def get_sigmaK_target_v3(P, Ms, K, sigP=5e-5, fracsigMs=.1):
    '''
    Compute the K detection significance required to measure a planet's mass at 
    3 sigma.
    '''
    P = unp.uarray(P, sigP)
    Ms = unp.uarray(Ms, fracsigMs*Ms)
    sigKs = np.logspace(-2, np.log10(K), 1000)
    detsig_mp = np.zeros(sigKs.size)
    for i in range(sigKs.size):
        Ktmp = unp.uarray(K, sigKs[i])
        mp = rvs.RV_mp(P, Ms, Ktmp)
        detsig_mp[i] = unp.nominal_values(mp) / unp.std_devs(mp)
    # Get where sigK returns a 3sigma mass detection
    fint = interp1d(detsig_mp, sigKs)
    return float(fint(3))
    

def save_results(planetindex, band_strs, mags, ra, dec, P, rp, mp, K, S, Ms,
                 Rs, Teff, dist, Prot, vsini, Z, sigmaRV_activity, sigmaRV_planets,
                 R, aperture_m, QE, fracsigmaK_target, sigmaRV_phot, sigmaRV_eff,
                 texp, tobserving, Nrv, fname):
    '''
    Write the stellar, planet, observatory parameters, and results to a text 
    file.
    '''
    try:
        os.mkdir('Results')
    except OSError:
        pass
    try:
        os.mkdir('Results/star%.4d'%planetindex)
    except OSError:
        pass

    f = open(fname, 'w')
    g = ''
    for i in range(len(mags)):
        g += '# %s = %.3f\n'%(band_strs[i], mags[i])

    g += '# ra (deg), dec (deg), P (days), rp (REarth), mp (MEarth), ' + \
         'K (m/s), S (SEarth), Ms (MSun), Rs (RSun), Teff (K), dist (pc), ' + \
         'Prot (days), vsini (km/s), Z ([Fe/H]), sigmaRV_act (m/s), ' + \
         'sigmaRV_planets (m/s), R (l/dl), aperture (m), QE, ' + \
         'fracsigmaK_target, sigmaRV_phot (m/s), sigmaRV_eff (m/s), texp (min), ' + \
         'tobs_tot (hrs), Nrv\n'
    g += '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%i\t%.3f\t%.3f\t%.3f\t%.1f\t%.3f\t%.3f\t%i\t%.2f\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.1f\n'%(ra, dec, P, rp, mp, K, S, Ms, Rs, Teff, dist, Prot, vsini, Z, sigmaRV_activity, sigmaRV_planets, R, aperture_m, QE, fracsigmaK_target, sigmaRV_phot, sigmaRV_eff, texp, tobserving, Nrv)
    f.write(g)
    f.close()
