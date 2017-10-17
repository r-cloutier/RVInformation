'''
Compute the information content in various spectral bands to estimate the 
RV uncertainty in each band using model spectra and the conversion factors 
from Artigau+2017.
'''
import numpy as np
import pylab as plt
import astropy.io.fits as fits
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import interp1d, UnivariateSpline
from PyAstronomy.pyasl import broadGaussFast, rotBroad


global c, h, SNRreference, bands
c, h, SNRreference = 299792458., 6.62607004e-34, 1e2
bands = ['u','b','v','r','i','Y','J','H','K']


def get_reduced_spectrum(Teff, logg, Z, vsini, band_str, R, pltt=False):
    '''
    Download a PHOENIX stellar model spectrum and reduce the spectrum over 
    a particular spectral bin via convolution with the instrumental resolution 
    and a stellar rotation kernel.

    Parameters
    ----------
    `Teff': scalar
        The stellar effective temperature in Kelvin
    `logg': scalar
        The stellar logg in cgs units
    `Z': scalar
        The stellar metallicity [Fe/H] in solar units
    `vsini': scalar
        The projected stellar rotation velocity in km/s
    `band_str': str
        The letter designating the spectral band under consideration. Must be 
        in ['u','b','v','r','i','Y','J','H','K']
    `R': scalar
        The spectral resolution of the spectrograph (lambda / d_lambda)
    `pltt': boolean
        If True, the fully reduced spectrum is plotted using pylab.show()

    Returns
    -------
    `wl_resamp': numpy.array
        Spectral array of wavelengths in microns
    `spec_resamp': numpy.array
        Stellar spectrum array in Nphotons/s/cm^2/cm

    '''
    wl = get_wavelengthgrid()
    _, spectrum = get_full_spectrum(float(Teff), float(logg), float(Z))
    wl_conv, spec_conv = _convolve_band_spectrum(wl, spectrum, band_str, R,
                                                 pltt=False)
    #TEMP
    #spec_conv = _rotational_convolution(wl_conv, spec_conv, vsini, pltt=pltt)
    wl_resamp, spec_resamp = _resample_spectrum(wl_conv, spec_conv, R)
    spec_scaled = cgs2Nphot(wl, spectrum, wl_resamp, spec_resamp)
    return wl_resamp, spec_scaled


def get_wavelengthgrid():
    '''
    Read-in the wavelength grid for the PHOENIX model spectra and return it.
    
    Returns
    -------
    `wl': numpy.array
        Full spectral array of wavelengths in microns (i.e. full spectral 
        coverage)

    '''
    dat = fits.open('input_data/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')[0]
    return np.ascontiguousarray(dat.data) * 1e-4


def get_full_spectrum(Teff, logg, Z):
    '''
    Read-in model spectrum from the PHOENIX library and return the header and 
    the full spectral data (i.e. no cut in wavelength).

    Parameters
    ----------
    `Teff': scalar
        The stellar effective temperature in Kelvin
    `logg': scalar
        The stellar logg in cgs units
    `Z': scalar
        The stellar metallicity [Fe/H] in solar units

    Returns
    -------
    `header': astropy.io.fits.header.Header
        The fits header of the PHOENIX stellar model
    `data': numpy.array
        Stellar spectrum array in erg/s/cm^2/cm

    '''
    # Define stellar values
    Teffs = np.append(np.arange(23e2,7e3,1e2), np.arange(7e3,121e2,2e2)) 
    assert Teff in Teffs
    loggs = np.arange(0,6.1,.5)
    assert logg in loggs
    Zs = np.append(np.arange(-4,-2,1), np.arange(-2,1.5,.5))
    assert Z in Zs

    # Download the spectrum
    prefix = "ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/" + \
             "PHOENIX-ACES-AGSS-COND-2011/"
    fname = "Z-%.1f/lte%.5d-%.2f-%.1f"%(Z, Teff, logg, Z) + \
            ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
    spec_fits = fits.open(prefix+fname)[0]

    return spec_fits.header, spec_fits.data


def _convolve_band_spectrum(wl_microns, spectrum, band_str, R, pltt=False):
    '''
    Convolve the spectrum in a given band with a Gaussian profile whose FWHM 
    is specified by the spectral resolution of the instrument.

    Parameters
    ----------
    `wl_microns': array-like
        Spectral array of wavelengths in microns
    `spectrum': array-like
        Stellar spectrum array in erg/s/cm^2/cm
    `band_str': str
        The letter designating the spectral band under consideration. Must be 
        in ['u','b','v','r','i','Y','J','H','K']
    `R': scalar
        The spectral resolution of the spectrograph (lambda / d_lambda)
    `pltt': boolean
        If True, the convolved and original spectra are plotted using 
        pylab.show()

    Returns
    -------
    `wl_band': numpy.array 
        Spectral array of wavelengths in microns over the spectral band
    `spec_band': numpy.array 
        Stellar spectrum array in erg/s/cm^2/cm over the spectral band

    '''
    # Create equidistant wl array
    wl_microns2 = np.linspace(wl_microns.min(), wl_microns.max(),
                              wl_microns.size)
    fint = interp1d(wl_microns, spectrum)
    spectrum2 = fint(wl_microns2)
    
    # Isolate wavelength range
    if band_str not in bands:
        raise ValueError('Unknown passband: %s'%band_str)
    wlmin, wlmax, wl_central_microns = get_band_range(band_str)
    in_band = (wl_microns2 >= wlmin) & (wl_microns2 <= wlmax)
    wl_band, spectrum_band = wl_microns2[in_band], spectrum2[in_band]

    # Convolve to instrument resolution
    FWHM_microns = wl_central_microns / float(R)
    sigma_microns = FWHM_microns / (2*np.sqrt(2*np.log(2)))
    print '\nConvolving the stellar spectrum to the instrumental resolution...'
    spectrum_conv = broadGaussFast(wl_band, spectrum_band, sigma_microns)
    
    if pltt:
        plt.plot(wl_band, spectrum_band, 'k-', label='Full')
        plt.plot(wl_band, spectrum_conv, 'b-', label='Convolved')
        plt.xlabel('Wavelength [microns]'), plt.ylabel('Flux [erg/s/cm2/cm]')
        plt.legend(), plt.show()

    return wl_band, spectrum_conv


def _rotational_convolution(wl_band, spec_band, vsini, epsilon=0.6,
                            pltt=False):
    '''
    Convolve the spectrum with rotational kernel based on the star's projected 
    rotation velocity and assuming a constant linear limb-darkening across the 
    stellar disk.

    Parameters
    ----------
    `wl_band': array-like
        Spectral array of wavelengths in microns
    `spec_band': array-like
        Stellar spectrum array in erg/s/cm^2/cm
    `band_str': str
        The letter designating the spectral band under consideration. Must be 
        in ['u','b','v','r','i','Y','J','H','K']
    `R': scalar
        The spectral resolution of the spectrograph (lambda / d_lambda)
    `pltt': boolean
        If True, the convolved and original spectra are plotted using 
        pylab.show()

    Returns
    -------
    `spec_conv': numpy.array 
        Stellar spectrum array in erg/s/cm^2/cm over the spectral band


    '''
    print '\nConvolving the stellar spectrum with a rotational profile...'
    spec_conv = rotBroad(wl_band, spec_band, epsilon, vsini)

    if pltt:
        plt.plot(wl_band, spec_band, 'k-', label='Original')
        plt.plot(wl_band, spec_conv, 'b-', label='Rotationally convolved')
        plt.xlabel('Wavelength [microns]'), plt.ylabel('Flux [erg/s/cm2/cm]')
        plt.legend(), plt.show()

    return spec_conv


def _resample_spectrum(wl, spec, R, pixels_per_element=3):
    '''
    Resample the input wavelength array and spectrum to the width of the 
    resolution element.
    
    Parameters
    ----------
    `wl': array-like
        Spectral array of wavelengths in microns
    `spec': array-like
        Stellar spectrum array in erg/s/cm^2/cm
    `R': scalar
        The spectral resolution of the spectrograph (lambda / d_lambda)
    `pixels_per_element': scalar
        The number of detector pixel in a single resolution element

    '''
    dl = wl.min() / (float(pixels_per_element) * R)
    wl_resamp = np.arange(wl.min(), wl.max(), dl)
    fint = interp1d(wl, spec)
    return wl_resamp, fint(wl_resamp)


def get_band_range(band_str):
    '''
    Define the wavelength range for a given band and its central wavelength.

    Parameters
    ----------
    `band_str': str
        The letter designating the spectral band under consideration. Must be 
        in ['u','b','v','r','i','Y','J','H','K']

    Returns
    -------
    `wlmin': float
        The smallest wavelength in the spectral band in microns
    `wlmax': float
        The largest wavelength in the spectral band in microns
    `wlcentral': float
        The smallest wavelength in the spectral band in microns


    '''
    # bin widths are chosen to have the same number of wavelength elements
    # per band for flux normalization purposes
    if band_str == 'u':
        wlmin, wlmax, wlcentral = 0.31, 0.42, .3656
    elif band_str == 'b':
        wlmin, wlmax, wlcentral = 0.37, 0.59, .4353    
    elif band_str == 'v':
        wlmin, wlmax, wlcentral = 0.47, 0.70, .5477
    elif band_str == 'r':
        wlmin, wlmax, wlcentral = 0.54, 0.87, .6349
    elif band_str == 'i':
        wlmin, wlmax, wlcentral = 0.70, 1.13, .8797
    elif band_str == 'Y':
        wlmin, wlmax, wlcentral = 1.0, 1.1, 1.02
    elif band_str == 'J':
        wlmin, wlmax, wlcentral = 1.17, 1.33, 1.22
    elif band_str == 'H':
        wlmin, wlmax, wlcentral = 1.5, 1.75, 1.63
    elif band_str == 'K':
        wlmin, wlmax, wlcentral = 2.07, 2.35, 2.19
    else:
        raise ValueError('Unknown bandpass: %s'%band_str)

    return wlmin, wlmax, wlcentral
    

def cgs2Nphot(wl_full_microns, spec_full_cgs, wl_band_microns, spec_band_cgs):
    '''
    Convert the input spectrum from cgs units (erg/s/cm^2/cm) to a the number 
    of photons with a fixed SNR at the center of the J band (1.25 microns).

    Parameters
    ----------
    `wl_full_microns': array-like
        Spectral array of wavelengths (in microns) across the full PHOENIX 
        spectral coverage
    `spec_full_cgs': array-like
        Stellar spectrum array (in erg/s/cm^2/cm) across the full PHOENIX 
        spectral coverage
    `wl_band_microns': array-like
        Spectral array of wavelengths (in microns) across a single spectral band
    `spec_band_cgs': array-like
        Stellar spectrum array (in erg/s/cm^2/cm) across a single spectral band
    
    Results
    -------
    `spec_Nphot_scaled': numpy.array
        Stellar spectrum array in Nphotons/s/cm^2/cm in a single spectral band

    '''
    assert wl_full_microns.size == spec_full_cgs.size
    assert wl_band_microns.size == spec_band_cgs.size
    wl_full_cm, wl_band_cm = wl_full_microns * 1e-4, wl_band_microns * 1e-4
    energy_full_erg, energy_band_erg = h*c / (wl_full_cm*1e-2) * 1e7, \
                                       h*c / (wl_band_cm*1e-2) * 1e7
    spec_full_Nphot, spec_band_Nphot = spec_full_cgs / energy_full_erg, \
                                       spec_band_cgs / energy_band_erg
    Jwl = 1.25
    centralJindex = np.where(abs(wl_full_microns-Jwl) == \
                             np.min(abs(wl_full_microns-Jwl)))[0][0]
    # SNR = sqrt(Nphot)
    norm = SNRreference**2 / \
           np.max(spec_full_Nphot[centralJindex-3:centralJindex+3])
    spec_Nphot_scaled = norm * spec_band_Nphot
    return spec_Nphot_scaled
                    

def _rescale_sigmaRV(sigmaRV, mag, band_str, texp_min, aperture_m, QE, R):
    '''
    Rescale sigmaRV from SNR=100 per resolution element to whatever SNR is 
    acheived in the input band over a given integration time.

    Parameters
    ----------
    `sigmaRV': scalar
        The photon-noise-limited RV precision in m/s
    `mag': scalar
        The stellar magnitude in spectral band given in `band_str'
    `band_str': str
        The letter designating the spectral band under consideration. Must be 
        in ['u','b','v','r','i','Y','J','H','K']
    `texp_min': scalar
        The integration time in minutes
    `aperture_m': float
        The telescope's aperture diameter in meters
    `QE': scalar
        The quantum efficiency of the detector (0<QE<=1)
    `R': scalar
        The spectral resolution of the spectrograph (lambda / d_lambda)

    Returns
    -------
    `SNR': float
        The rescaled SNR of the spectrum

    '''
    snr = _get_snr(mag, band_str, texp_min, aperture_m, QE, R)
    return sigmaRV * np.sqrt(SNRreference / snr)
    

def _get_snr(mag, band_str, texp_min, aperture_m, QE, R):
    '''
    Compute the SNR of the spectrum in a certain band (e.g. 'J').

    Parameters
    ----------
    `mag': scalar
        The stellar magnitude in spectral band given in `band_str'
    `band_str': str
        The letter designating the spectral band under consideration. Must be 
        in ['u','b','v','r','i','Y','J','H','K']
    `texp_min': scalar
        The integration time in minutes
    `aperture_m': float
        The telescope's aperture diameter in meters
    `QE': scalar
        The quantum efficiency of the detector (0<QE<=1)
    `R': scalar
        The spectral resolution of the spectrograph (lambda / d_lambda)

    Returns
    -------
    `SNR': float
        The SNR of the spectrum in a certain band

    '''
    # Define constants
    area_cm2 = np.pi * (aperture_m*1e2/2)**2
    res_kms = c/R
    texp_s = texp_min * 60.
    c_As = c * 1e10
    c_kms = c_As * 1e-13
    h_ergss = h * 1e7

    # Get flux density zeropoint (for m=0) in ergs s^-1 A^-1 cm^-2,
    # wavelength in angstroms and bandwidth in microns
    if band_str == 'u':
        Fl = 3.639e-9
        l  = 3594.93
    elif band_str == 'g':
        Fl = 5.521e-9
        l  = 4640.42
    elif band_str == 'r':
        Fl = 2.529e-9
        l  = 6122.33
    elif band_str == 'i':
        Fl = 1.409e-9
        l  = 7439.49
    elif band_str == 'z':
        Fl = 8.501e-10
        l  = 8897.06
    elif band_str == 'Y':
        Fl = 6.063e-10
        l  = 10174.
    elif band_str == 'J':
        Fl = 3.143e-10
        l  = 12350.
    elif band_str == 'H':
        Fl = 1.144e-10
        l  = 16620.
    elif band_str == 'K':
        Fl = 4.306e-11
        l  = 21590.
    else:
        raise ValueError('Unknown passband: %s'%band_str)

    fl = Fl * 10**(-.4 * mag)
    Ephot = h_ergss * c_As / l
    Nphot_A = fl * texp_s * area_cm2 * QE / Ephot
    
    # Get # photons per resolution element (dl/l)
    Nphot_res = Nphot_A * (l / R)

    # Add photon noise and dark current to noise budget
    darkcurrent, footprint = 1e-2, 12     # electrons/s, pixels
    SNR = Nphot_res / np.sqrt(Nphot_res + darkcurrent*footprint*texp_s)

    return SNR


def _remove_tellurics_from_W(wl_band, W, transmission_threshold=.02):
    '''
    Remove wavelengths that are sampled at wavelengths affected by tellurics at the 
    level more than a specified threshold.

    Parameters
    ----------
    `wl_band': array-like
        Spectral array of wavelengths in microns
    `W': array-like
        Spectral weighting function from Eq. X in Bouchy et al 2001 in Nphot/s/cm^2/cm
    `transmission_threshold': scalar
        Maximum fraction absorption from tellurics. Only keep where transmission is 
        greater than 1-`transmission_threshold'

    Returns
    -------
    `Wout': numpy.array()
        Spectral weighting function only where the atmospheric transmission is 
        favourable

    '''
    assert wl_band.size == W.size
    wlTAPAS, transTAPAS = np.loadtxt('input_data/tapas_000001.ipac', \
                                     skiprows=23).T
    wlTAPAS *= 1e-3

    # remove rayleigh continuum via boxcar smoothing if passband is in that regime
    if np.any(wl_band <=.8):
        boxsteps = np.arange(wlTAPAS.min(), wlTAPAS.max(), 1e-2)
        transTAPAS_continuum = np.zeros(boxsteps.size-1)
    	for i in range(boxsteps.size-1):
	    transTAPAS_continuum[i] = np.max(transTAPAS[(wlTAPAS >= boxsteps[i]) \
                                                        & (wlTAPAS <= boxsteps[i+1])])
	fint = interp1d(boxsteps[1:]-np.diff(boxsteps)[0]*.5, transTAPAS_continuum, \
                        bounds_error=False, fill_value=.878)
	transTAPAS = transTAPAS - fint(wlTAPAS) + 1.

    # resample the wavelength grid
    fint = interp1d(wlTAPAS, transTAPAS)
    transmission = fint(wl_band)

    # only keep where transmission > 1-threshold
    notellurics = np.where(transmission > (1.-transmission_threshold))[0]
    return W[notellurics]


def exposure_time_calculator_per_band(mag, band_str, aperture_m, QE, R,
                                      texpmin=10, texpmax=60, SNRtarget=150):
    '''
    Compute the exposure time to recieve a SNR target per resolution element in 
    a particular band.
    '''
    texps = np.arange(texpmin, texpmax+.1, .1)  # minutes
    SNRs = np.zeros(texps.size)
    for i in range(texps.size):
        SNRs[i] = _get_snr(mag, band_str, texps[i], aperture_m, QE, R)

    g = abs(SNRs-SNRtarget) == abs(SNRs-SNRtarget).min()
    return texps[g].min()



def _compute_W(wl_band, spec_band):
    '''
    Compute the weighting function F from the input spectrum using the Eqs. 
    from Bouchy+2001.
    '''
    # over-sample the spectrum
    wl2 = np.linspace(wl.min(), wl.max(), wl.size*10)
    fint = interp1d(wl, spec)
    A = fint(wl2)

    # compute W
    dwl = np.diff(wl2)[0]
    W = wl2**2 * np.gradient(A, dwl)**2 / A
    
    # resample W
    fint = interp1d(wl2, W)
    return fint(wl)


def compute_sigmaRV(wl, spec, mag, band_str, texp_min=5, aperture_m=3.58,
                    QE=.15, R=75e3):
    W = _compute_W(wl, spec)
    # remove tellurics
    W_clean = _remove_tellurics_from_W(wl, W)
    g = np.arange(4, W_clean.size-4, dtype=int)
    sigmaRV = c / np.sqrt(np.sum(W_clean[g]))
    sigmaRV_scaled = _rescale_sigmaRV(sigmaRV, mag, band_str, texp_min,
                                      aperture_m, QE, R)
    return sigmaRV_scaled
