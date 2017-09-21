'''
Compute the information content in various spectral bands to estimate the 
RV uncertainty in each band using model spectra and the conversion factors 
from Artigau+2017.
'''
import numpy as np
import pylab as plt
import astropy.io.fits as fits
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import interp1d
from PyAstronomy.pyasl import broadGaussFast


global c, h
c, h = 299792458., 6.62607004e-34


def get_wavelengthgrid():
    '''
    Read-in the wavelength grid for the PHOENIX model spectra and return it in
    microns.
    '''
    dat = fits.open('input_data/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')[0]
    return np.ascontiguousarray(dat.data) * 1e-4


def get_spectrum(Teff, logg, Z):
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
    spectrum_conv = broadGaussFast(wl_band, spectrum_band, sigma_microns)

    if False:
        plt.plot(wl_band, spectrum_band, 'k-')
        plt.plot(wl_band, spectrum_conv, 'b-'), plt.show()

    return wl_band, spectrum_conv


def get_reduced_spectrum(Teff, logg, Z, band_str, R=75e3):
    '''
    Download a stellar spectrum and get the reduced spectrum in the spectral 
    bin of interest.
    '''
    wl = get_wavelengthgrid()
    _, spectrum = get_spectrum(Teff, logg, Z)
    wl_conv, spec_conv = _convolve_band_spectrum(wl, spectrum, band_str, R=R)
    return wl_conv, spec_conv


def _compute_Q(wl, spec):
    '''
    Compute the quality factor Q from the input spectrum using the Eqs. from 
    Bouchy+2001.
    '''
    A = spec
    W = wl**2 * np.gradient(A, np.diff(wl)[0])**2 / A
    return np.sqrt(np.sum(W)) / np.sqrt(np.sum(A))


def _compute_Nelectrons(wl_microns, spec,
                        texp_min=5., diameter_m=3.5, QE=.1):
    '''
    Compute the number of photoelectrons collected.
    '''
    energy_ergs = h*c / (wl_microns*1e-6) * 1e7
    texp_sec, A_cm2 = texp_min*6e1, np.pi*(diameter_m*1e2/2.)**2
    Nphot = spec / energy_ergs * texp_sec * A_cm2 * wl_microns*1e-4
    return Nphot * QE
    

def compute_sigmaRV():
    return None
