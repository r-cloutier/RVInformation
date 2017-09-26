'''
Compute nRV for all TESS targets using the code cloned from
https://github.com/r-cloutier/Transiting_RVcalculator.git
'''
from get_tess_data import get_TESS_data
from compute_sigmaRV import compute_sigRV, get_reduced
import pylab as plt
from MRrelations import WM14
from nRVcalculator import nRVcalculator


global R, aperture, QE
R, aperture, QE = 75e3, 3.58, .1


ra,dec,rp,P,S,K,Rs,Teff,Vmag,Imag,Jmag,Kmag,dist,_,_,snr,_ = get_TESS_data()
nplanets = ra.size

for i in range(nplanets):

    # compute the expected sigmaRV
    wl, spec = get_reduced_spectrum(Teff, logg, Z, band_str, mag, texp_min,
                                    aperture_m=aperture, QE=QE, R=R)
    sigmaRV = compute_sigmaRV(wl, spec)

    startheta = get_startheta()

    planettheta = (P[i], sigP, rp[i])

    jitterRMS = get_jitterRMS()

    self = nRVcalculator(startheta, planettheta, MRfunc,
                         detsigs=[3,5], sigmaRV=sigmaRV,
                         bands=['Y','J','H','K'], texp=texp_min,
                         additiveRVjitter=0,
                         R=R, aperture=aperture, efficiency=QE, SNRtarget=150)


def get_startheta(Rs, Teff, dist, K, P, rp):
    '''
    Return the stellar parameters based on the results from the 
    simulated Sullivan stellar sample.
    '''
    rvs.RV_mp(P, Ms, K)

    # mp from rp & WM14
    # Ms from K, P, mp, e=0, i=90 deg
    # SpT from Ms
    # mags from dist and SED of SpT
    # vsini from Rs and sampled Prot
    
    return mags, SpT, vsini, Ms, sigMs


# set desired sigma_mp based on sigma_rho / rho = 0.2
# compute corresponding sigma_K from (K, P, sigma_P, M, sigma_M, mp, sigma_mp)
# K, P, rp from Sullivan
# mp from rp and WM14 mean relation
# set sigma_P == constant
# set sigma_M = .1 * M
# M from K, P, mp
# compute sigma_eff from sigma_RV plus other planets and activity
# compute sigma_RV from quality factor based on stellar and instrument specs

