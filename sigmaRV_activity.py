import numpy as np
import rvs


def F82sigmaRV(F8, Teff):
    '''
    Convert the F8 flicker to an expected sigmaRV using Eqs. 3 & 4 in Cegla et 
    al. 2014.

    Parameters
    ----------
    `F8': array-like
        The 8-hour flicker or photometric variability on a timescale of 8 hours
        measured in parts-per-thousand.
    `Teff': array-like
        The stellar effective temperature in kelvin

    Returns
    -------
    `sigmaRVs': numpy.array
        The expected RV noise from activity as probed by the 8-hour photometric 
        flicker in m/s

    '''
    F8, Teff = np.ascontiguousarray(F8), np.ascontiguousarray(Teff)
    assert F8.size == Teff.size
    sigmaRVs = np.zeros(F8.size)

    g1 = Teff >= 6e3
    sigmaRVs[g1] = 84.23 * F8[g1] - 3.35

    g2 = Teff < 6e3
    sigmaRVs[g2] = 18.04 * F8[g2] - .98

    return sigmaRVs


def logRhk2sigmaRV(logRhk, Teff, mean_only=False):
    '''
    Convert the logRhk chromospheric index to an expected sigmaRV using Eqs. 2, 
    3, and 4 from Santos et al. 2000.

    Parameters
    ----------
    `logRhk': array-like
        The chromospheric activity index logRhk
    `Teff': array-like
        The stellar effective temperature in kelvin
    `mean_only': boolean
        If True, the residual rms in the reported relation is not used to 
        resample the mean sigmaRV

    Returns
    -------
    `sigmaRVs': numpy.array
        The expected RV noise from activity as probed by logRhk in m/s
    
    '''
    R5, Teff = np.ascontiguousarray(1e5*10**logRhk), np.ascontiguousarray(Teff)
    assert R5.size == Teff.size
    sigmaRVs = np.zeros(R5.size)

    g1 = (Teff < 6650) & (Teff >= 6e3)
    sigmaRVs[g1] = 9.2 * R5[g1]**(.75)
    
    g2 = (Teff >= 52e2) & (Teff < 6e3)
    sigmaRVs[g2] = 7.9 * R5[g2]**(.55)

    g3 = Teff < 52e2
    sigmaRVs[g3] = 7.8 * R5[g3]**(.12)

    if not mean_only:
        sigmaRVs[g1] += np.random.randn(sigmaRVs[g1].size) * 10**(.17)
        sigmaRVs[g2] += np.random.randn(sigmaRVs[g2].size) * 10**(.18)
        sigmaRVs[g3] += np.random.randn(sigmaRVs[g3].size) * 10**(.19)
        
    return sigmaRVs


def Prot2logRhk(Prot, Ms, mean_only=False):
    '''
    Convert the stellar rotation period to the logRhk chromospheric index for 
    stars with masses less than 0.8 solar masses from Astudillo-Defru et al. 
    2017.

    Parameters
    ----------
    `Prot': array-like
        The stellar rotation period in days
    `Ms': array-like
        The stellar mass in solar masses
    `mean_only': boolean
        If True, the residual rms in the reported relation is not used to 
        resample the mean logRhk

    Returns
    -------
    `logRhk': numpy.array
        The chromospheric activity index logRhk
    
    '''
    Prot, Ms = np.ascontiguousarray(Prot), np.ascontiguousarray(Ms)
    assert Prot.size == Ms.size
    if np.any(Ms > .8):
        raise ValueError("Relationship is only valid for stars with masses " + \
                         "less than or equal to 0.8 solar masses.")
    logRhks = np.zeros(Ms.size)

    sigm, sigb, siglogRhk = .007, .02, .093 if not mean_only else 0, 0, 0
    m = -1.509 + np.random.randn(Ms.size) * sigm
    b = -2.550 + np.random.randn(Ms.size) * sigb

    g1 = Prot >= 10    
    logRhks[g1] = m[g1] * np.log10(Prot[g1]) + b[g1]

    g2 = Prot < 10
    logRhks[g2] = -4.045 + np.random.randn(Ms[g2].size) * siglogRhk

    return logRhks



# logRhk for low-mass stars: http://adsabs.harvard.edu/abs/2017A%26A...600A..13A
# lgoRhk for FGK stars: http://adsabs.harvard.edu/abs/2011arXiv1107.5325L 
