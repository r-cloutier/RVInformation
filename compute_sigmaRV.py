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


global c, h, bands
c, h = 299792458., 6.62607004e-34
bands = ['u','g','r','i','z','Y','J','H','K']


def get_reduced_spectrum(Teff, logg, Z, vsini, band_str, R, pltt=False):
    '''
    Download a PHEONIX stellar model spectrum and reduce the spectrum over 
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
        in ['u','g','r','i','z','Y','J','H','K']
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
                                                 pltt=pltt)
    #TEMP
    #spec_conv = _rotational_convolution(wl_conv, spec_conv, vsini, pltt=pltt)
    wl_resamp, spec_resamp = _resample_spectrum(wl_conv, spec_conv, R)
    spec_scaled = _cgs2Nphot(wl, spectrum, wl_resamp, spec_resamp)
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
        The fits header of the PHEONIX stellar model
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
    Convolve the spectrum in a given band with a Gaussian profile with a 
    FWHM specified by the spectral resolution of the instrument.
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


def _rotational_convolution(wl, spec, vsini, epsilon=0.6, pltt=False):
    '''
    Convolve the spectrum with Gaussian rotational profile based on the 
    star's projected rotation velocity and assuming a constant linear 
    limb-darkening across the stellar disk.
    '''
    print '\nConvolving the stellar spectrum with a rotational profile...'
    spec_conv = rotBroad(wl, spec, epsilon, vsini)
    if pltt:
        plt.plot(wl, spec, 'k-', label='Original')
        plt.plot(wl, spec_conv, 'b-', label='Rotationally convolved')
        plt.xlabel('Wavelength [microns]'), plt.ylabel('Flux [erg/s/cm2/cm]')
        plt.legend(), plt.show()
    return spec_conv


def _resample_spectrum(wl, spec, R):
    '''
    Resample the input spectrum to the size of the resolution element.
    '''
    pixels_per_element = 3.
    dl = wl.min() / (pixels_per_element * R)
    wl_resamp = np.arange(wl.min(), wl.max(), dl)
    fint = interp1d(wl, spec)
    return wl_resamp, fint(wl_resamp)


def get_band_range(band_str):
    '''
    Define wavelength range for a given band and return the indices that lie 
    within the band.
    '''
    if band_str == 'u':
        wlmin, wlmax, wlcentral = .31, .42, .3656
    elif band_str == 'b':
        wlmin, wlmax, wlcentral = .37, .54, .4353    
    elif band_str == 'v':
        wlmin, wlmax, wlcentral = .47, .69, .5477
    elif band_str == 'r':
        wlmin, wlmax, wlcentral = .54, .85, .6349
    elif band_str == 'i':
        wlmin, wlmax, wlcentral = .70, 1.10, .8797
    elif band_str == 'Y':
        wlmin, wlmax, wlcentral = 1.0, 1.1, 1.01743
    elif band_str == 'J':
        wlmin, wlmax, wlcentral = 1.17, 1.33, 1.2350
    elif band_str == 'H':
        wlmin, wlmax, wlcentral = 1.50, 1.75, 1.6620
    elif band_str == 'K':
        wlmin, wlmax, wlcentral = 2.07, 2.35, 2.1590
    else:
        raise ValueError('Unknown bandpass: %s'%band_str)

    return wlmin, wlmax, wlcentral
    

def _cgs2Nphot(wl_full_microns, spec_full_cgs, wl_band_microns, spec_band_cgs,
               SNR=1e2, Jwl=1.25):
    '''
    Get the conversion factor to convert the spectrum from cgs units to a 
    number of photons with a SNR=100 at the center of the J band (1.25 microns).
    '''
    assert wl_full_microns.size == spec_full_cgs.size
    assert wl_band_microns.size == spec_band_cgs.size
    wl_full_cm, wl_band_cm = wl_full_microns * 1e-4, wl_band_microns * 1e-4
    energy_full_erg, energy_band_erg = h*c / (wl_full_cm*1e-2) * 1e7, \
                                       h*c / (wl_band_cm*1e-2) * 1e7
    spec_full_Nphot, spec_band_Nphot = spec_full_cgs / energy_full_erg, \
                                       spec_band_cgs / energy_band_erg
    centralJindex = np.where(abs(wl_full_microns-Jwl) == \
                             np.min(abs(wl_full_microns-Jwl)))[0][0]
    # SNR = sqrt(Nphot)
    norm = SNR**2 / np.max(spec_full_Nphot[centralJindex-3:centralJindex+3])
    spec_Nphot_scaled = norm * spec_band_Nphot
    return spec_Nphot_scaled
                    

def _rescale_SNR(sigmaRV, mag, band_str, texp_min, aperture_m, QE, R):
    '''
    Rescale sigmaRV from SNR=100 per resolution element to whatever SNR is 
    acheived in the input band over a given integration time.
    '''
    snr = _get_snr(mag, band_str, texp_min, aperture_m, QE, R)
    return sigmaRV * np.sqrt(1e2 / snr)
    

def _get_snr(mag, band_str, texp_min, aperture_m, QE, R):
    '''
    Compute the SNR of the spectrum from the apparent magnitude of
    the star in a certain band (e.g. 'J'), the exposure time in
    minutes, the aperture of the telescope in meters, detector efficiency,
    the element resolution in km/s.
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


def _remove_tellurics(wl, W):
    '''
    Remove wavelengths samples that are affected by tellurics at the 
    level of > 2%.
    **To be used before computing sigmaRV from W**
    '''
    assert wl.size == W.size
    wlTAPAS, transTAPAS = np.loadtxt('input_data/tapas_000001.ipac', \
                                     skiprows=23).T
    wlTAPAS *= 1e-3
    # remove rayleigh continuum via boxcar smoothing if passband is in
    # that regime
    if np.any(wl <=.8):
        boxsteps = np.arange(wlTAPAS.min(), wlTAPAS.max(), 1e-2)
        transTAPAS_continuum = np.zeros(boxsteps.size-1)
    	for i in range(boxsteps.size-1):
	    transTAPAS_continuum[i] = np.max(transTAPAS[(wlTAPAS >= boxsteps[i]) \
                                                        & (wlTAPAS <= boxsteps[i+1])])
	fint = interp1d(boxsteps[1:]-np.diff(boxsteps)[0]*.5, transTAPAS_continuum, \
                        bounds_error=False, fill_value=.878)
	transTAPAS = transTAPAS - fint(wlTAPAS) + 1.
    # resample wavelength grid
    fint = interp1d(wlTAPAS, transTAPAS)
    transmission = fint(wl)
    # only keep where transmission is > 98%
    notellurics = np.where(transmission > .98)[0]
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



def _compute_W(wl, spec):
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
    W_clean = _remove_tellurics(wl, W)
    g = np.arange(4, W_clean.size-4, dtype=int)
    sigmaRV = c / np.sqrt(np.sum(W_clean[g]))
    sigmaRV_scaled = _rescale_SNR(sigmaRV, mag, band_str, texp_min, aperture_m,
                                  QE, R)
    return sigmaRV_scaled
