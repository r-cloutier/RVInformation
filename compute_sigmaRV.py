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
#from PyAstronomy.pyasl import broadGaussFast


global c, h
c, h = 299792458., 6.62607004e-34


def get_reduced_spectrum(Teff, logg, Z, band_str, mag, texp_min,
                         aperture_m=3.58, QE=.1, R=75e3):
    '''
    Download a stellar spectrum and get the reduced spectrum in the spectral 
    bin of interest.
    '''
    wl = _get_wavelengthgrid()
    _, spectrum = _get_full_spectrum(Teff, logg, Z)
    wl_conv, spec_conv = _convolve_band_spectrum(wl, spectrum, band_str, R=R)
    spec_scaled = _rescale_spectrum(spec_conv, mag, band_str, texp_min,
                                    aperture_m, QE, R)
    return wl_conv, spec_scaled


def _get_wavelengthgrid():
    '''
    Read-in the wavelength grid for the PHOENIX model spectra and return it in
    microns.
    '''
    dat = fits.open('input_data/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')[0]
    return np.ascontiguousarray(dat.data) * 1e-4


def _get_full_spectrum(Teff, logg, Z):
    '''
    Read-in model spectra from the PHOENIX library and return the header and 
    the spectral data.
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


def _convolve_band_spectrum(wl_microns, spectrum, band_str, R=75e3):
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
    wl_band, transmission, wl_central_microns = _get_band_transmission(band_str)
    in_band = (wl_microns2 >= wl_band.min()) & (wl_microns2 <= wl_band.max())
    wl_band, spectrum_band = wl_microns2[in_band], spectrum2[in_band]

    # Convolve to instrument resolution
    FWHM_microns = wl_central_microns / float(R)
    sigma_microns = FWHM_microns / (2*np.sqrt(2*np.log(2)))
    #spectrum_conv = broadGaussFast(wl_band, spectrum_band, sigma_microns)
    spectrum_conv = spectrum_band ## TEMP
    
    
    if False:
        plt.plot(wl_band, spectrum_band, 'k-')
        plt.plot(wl_band, spectrum_conv, 'b-'), plt.show()

    return wl_band, spectrum_conv


def _get_band_transmission(band_str):
    '''
    Read-in the transmission curve of a specific band on the same wavelength 
    grid as the PHOENIX spectra.
    '''
    if band_str == 'J':
        fname = '2MASS-2MASS.J.dat'
    elif band_str == 'H':
        fname = '2MASS-2MASS.H.dat'
    elif band_str == 'Ks':
        fname = '2MASS-2MASS.Ks.dat'
    else:
        raise ValueError('Do not know what input bandpass %s is.'%band_str)

    wl, transmission = np.loadtxt('input_data/%s'%fname).T
    wl *= 1e-4   # angstrom -> microns

    return wl, transmission, np.average(wl, weights=transmission)


def _rescale_spectrum(spec, mag, band_str, texp_min, aperture_m, QE, R):
    '''
    Rescale the spectrum such that its sum over the band equals the SNR.
    spec` = k*spec where k = snr^2 / sum(spec) since sum(spec) = Nphot = snr^2
    '''
    snr = _get_snr(mag, band_str, texp_min, aperture_m, QE, R)
    factor = snr**2 / np.sum(spec)
    return factor * spec


def _get_snr(mag, band, texp_min, aperture_m, QE, R):
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
    if band == 'Y':
        Fl = 6.063e-10
        l  = 10174.
    elif band == 'J':
        Fl = 3.143e-10
        l  = 12350.
    elif band == 'H':
        Fl = 1.144e-10
        l  = 16620.
    elif band == 'K':
        Fl = 4.306e-11
        l  = 21590.
    else:
        raise ValueError('Not sure which band to use')

    fl = Fl * 10**(-.4 * mag)
    Ephot = h_ergss * c_As / l
    Nphot_A = fl * texp_s * area_cm2 * QE / Ephot
    
    # Get # photons per resolution element (dl/l)
    Nphot_res = Nphot_A * (l / R)

    # Add photon noise and dark current to noise budget
    darkcurrent, footprint = 1e-2, 12     # electrons/s, pixels
    SNR = Nphot_res / np.sqrt(Nphot_res + darkcurrent*footprint*texp_s)

    return SNR



###############################################################################
###############################################################################
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


def compute_sigRV(wl, spec):
    A = spec
    W = wl**2 * np.gradient(A, np.diff(wl)[0])**2 / A
    g = np.arange(4, W.size-4, dtype=int)
    return c / np.sqrt(np.sum(W[g]))



