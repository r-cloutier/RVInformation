from imports import *
from compute_sigmaRV import *
from astropy.io import fits


def exposure_time_calculator_per_band(mags, band_strs, aperture_m, QE, R,
                                      texpmin=10, texpmax=60, SNRtarget=1e2):
    '''
    Compute the exposure time required to reach a target SNR per resolution 
    element in a particular band. Either V for visible spectrographs or J for 
    nIR spectrographs.
    
    Parameters
    ----------
    `mags': list of scalars
        A list of the stellar magnitudes for the star in each spectral band
    `band_strs': list of strs
        A list of the spectral bands that span the wavelength coverage of the 
        spectrograph used to measure the radial velocities of the TESS star. 
        All band_strs entries must be in 
        ['U','B','V','R','I','Y','J','H','K']
        and at least one of ['V','J'] as a reference band
    `aperture_m': float
        The telescope's aperture diameter in meters
    `QE': scalar
        The quantum efficiency of the detector (0<QE<=1)
    `R': scalar
        The spectral resolution of the spectrograph (lambda / d_lambda)
    `texpmin': scalar
        The minimum exposure time in minutes. Required to mitigate the effects 
        of stellar pulsations and granulation
    `texpmax': scalar
        The maximum exposure time in minutes. Used to moderate the limit the 
        observational time that can be dedicated to a single star
    
    Returns
    -------
    `texp': float
        The minimum exposure time in minutes required to achieve a desired SNR 
        in a particular reference band (either 'V' or 'J')

    '''
    if 'V' in band_strs:
        reference_band = 'V'
    elif 'J' in band_strs:
        reference_band = 'J'
    else:
        raise ValueError("No reference band. Neither 'V' nor 'J' are " + \
                         'included in band_strs.')
    mags, band_strs = np.ascontiguousarray(mags), \
                      np.ascontiguousarray(band_strs)
    reference_mag = float(mags[band_strs == reference_band])
    
    texps = np.arange(texpmin, texpmax+.1, .1)  # minutes
    SNRs = np.zeros(texps.size)
    for i in range(texps.size):
        SNRs[i] = get_snr(reference_mag, reference_band, texps[i],
                          aperture_m, QE, R)

    if SNRs.min() > SNRtarget:
        return float(texpmin)
    elif SNRs.max() < SNRtarget:
        return float(texpmax)
    else:
        fint = interp1d(SNRs, texps)
        return float(fint(SNRtarget))


def _get_band_spectra(band_strs, Teff_round, logg_round, Z=0., vsini=.1, R=7e4):
    
    for i in range(len(band_strs)):
        wl, spec = get_reduced_spectrum(Teff_round, logg_round, Z, vsini,
                                        band_strs[i], R, pltt=False)
        hdu = fits.HDUList()
        hdu.append(fits.ImageHDU(wl))
        hdu.append(fits.ImageHDU(spec))
        hdu.writeto('input_data/band_spectrum_%s'%band_strs[i], overwrite=True)
    
    
def compute_texps(mags, band_strs, aperture_m, QE, R, Teff_round):

    wlTAPAS, transTAPAS = np.loadtxt('input_data/tapas_000001.ipac', \
                                     skiprows=23).T
    wlTAPAS *= 1e-3
    transmission_threshold = .02
    
    SNRtargets = np.linspace(50, 5e2, 100)
    texps, sigmaRV_phots = np.zeros(SNRtargets.size), np.zeros(SNRtargets.size) 

    for i in range(texps.size):
        print float(i) / texps.size
    	texps[i] = exposure_time_calculator_per_band(mags, band_strs,
                                                     aperture_m, QE, R,
                                                     SNRtarget=SNRtargets[i])
        sigmaRVs = np.zeros(mags.size)
        for j in range(mags.size):
            hdu = fits.open('input_data/band_spectrum_%s'%band_strs[j])
            wl_band, spec_band = hdu[0].data, hdu[1].data
            ##
            W = compute_W(wl_band, spec_band)
            if np.any(wl_band <=.8):
                boxsteps = np.arange(wlTAPAS.min(), wlTAPAS.max(), 1e-2)
                transTAPAS_continuum = np.zeros(boxsteps.size-1)
                for i in range(boxsteps.size-1):
                    transTAPAS_continuum[i] = np.max(transTAPAS[(wlTAPAS >= \
                                                                 boxsteps[i]) \
                                                                & (wlTAPAS <= \
                                                                   boxsteps[i+1])])
                fint = interp1d(boxsteps[1:]-np.diff(boxsteps)[0]*.5,
                                transTAPAS_continuum,
                                bounds_error=False, fill_value=.878)
                transTAPAS = transTAPAS - fint(wlTAPAS) + 1.

            # resample the wavelength grid
            fint = interp1d(wlTAPAS, transTAPAS, bounds_error=False,
                            fill_value=0)
            transmission = fint(wl_band)

            # only keep where transmission > 1-threshold
            notellurics = np.where(transmission > (1.-transmission_threshold))[0]
            W_clean =  W[notellurics]
            ##
            g = np.arange(W_clean.size)
            sigmaRV = c / np.sqrt(np.sum(W_clean[g]))
            sigmaRVs[j] = rescale_sigmaRV(sigmaRV, mags[j], band_strs[j],
                                          texps[i], aperture_m, QE, R)
            
        correctionsYJHK = np.array([.47, .63, 1.59, 1.72])
        sigmaRVs[band_strs == 'Y'] = sigmaRVs[band_strs == 'Y'] / \
                                     correctionsYJHK[0]
        sigmaRVs[band_strs == 'J'] = sigmaRVs[band_strs == 'J'] / \
                                     correctionsYJHK[1]
        sigmaRVs[band_strs == 'H'] = sigmaRVs[band_strs == 'H'] / \
                                     correctionsYJHK[2]
        sigmaRVs[band_strs == 'K'] = sigmaRVs[band_strs == 'K'] / \
                                     correctionsYJHK[3]
        sigmaRV_phots[i] = 1. / np.sqrt(np.sum(1./sigmaRVs**2))

    # save
    np.savetxt('testsigmaRVphot.dat',
               np.array([SNRtargets,texps,sigmaRV_phots]).T,
               fmt='%.5e', delimiter='\t')
    return SNRtargets, texps, sigmaRV_phots

        

if __name__ == '__main__':
    band_strs, mags = np.array(['Y','J','H','K']), \
                      np.array([9.99, 9.5, 8.93, 8.695])
    aperture_m, QE, R = 3.6, .1, 7e4
    SNRtargets, texps, sigmaRV_phots = compute_texps(mags, band_strs,
                                                     aperture_m, QE, R)
