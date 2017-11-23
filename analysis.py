from RVInformationResults import *
import runTESS, rvs
from uncertainties import unumpy as unp
from get_tess_data import *


def _get_inds_of_complete_systems(self):
    starnums = np.unique(self.starnums)
    systnums = np.unique(self.systnums)
    files2keep = []
    for i in starnums:
        for j in systnums:
            if self.files[(self.starnums==i) & (self.systnums==j)].size == 3:
                for f in self.files[(self.starnums==i) & (self.systnums==j)]:
                    files2keep.append(f)
    files2keep = np.ascontiguousarray(files2keep)
    return np.in1d(self.files, files2keep)    
    

def TESS_mp_3sigma(self):
    '''
    Report the number of observations required to detect all TESS planet masses 
    at 3sigma. For this, sigK / K should be 0.327. 
    '''
    inds = _get_inds_of_complete_systems(self)
    starnums = np.unique(self.starnums[inds]).astype(int)
    Nharps, Nnirps, Nspirou = np.zeros(starnums.size), \
                              np.zeros(starnums.size), \
                              np.zeros(starnums.size)
    texpharps, texpnirps, texpspirou = np.zeros(starnums.size), \
                                       np.zeros(starnums.size), \
                                       np.zeros(starnums.size)
    tobsharps, tobsnirps, tobsspirou = np.zeros(starnums.size), \
                                       np.zeros(starnums.size), \
                                       np.zeros(starnums.size)

    for i in range(starnums.size):
        
        h = (self.starnums[inds] == starnums[i]) & \
            (self.spectrographs[inds] == 'H')
        Nharps[i]  = np.median(self.Nrvs[inds][h])
        texpharps[i] = self.texps[inds][h][0]
        tobsharps[i] = Nharps[i] * texpharps[i] / 60

        n = (self.starnums[inds] == starnums[i]) & \
            (self.spectrographs[inds] == 'N')
        Nnirps[i]  = np.median(self.Nrvs[inds][n])
        texpnirps[i] = self.texps[inds][n][0]
        tobsnirps[i] = Nnirps[i] * texpnirps[i] / 60
        
        s = (self.starnums[inds] == starnums[i]) & \
            (self.spectrographs[inds] == 'S')
        Nspirou[i] = np.median(self.Nrvs[inds][s])
        texpspirou[i] = self.texps[inds][s][0]
        tobsspirou[i] = Nspirou[i] * texpspirou[i] / 60
        
    # which spectrograph is the most efficient for each planet?
    Nrvs = np.array([Nharps, Nnirps, Nspirou]).T
    min_Nrv = np.min(Nrvs, axis=1)    
    bestspectrograph_Nrv = np.argmin(Nrvs, axis=1)    
    tobss = np.array([tobsharps, tobsnirps, tobsspirou]).T
    min_tobs = np.min(tobss, axis=1)
    bestspectrograph_tobs = np.argmin(tobss, axis=1)

    # save results to file
    hdr = 'TESS planet index\nNrv_HARPS\nNrv_NIRPS\nNrv_SPIROU\ntexp_HARPS ' + \
          '[min]\ntexp_NIRPS [min]\ntexp_SPIROU [min]\ntobs_HARPS [hrs]\n' + \
          'tobs_NIRPS [hrs]\ntobs_SPIROU [hrs]\nmin Nrv\nmin spectrograph ' + \
          '(0=HARPS, 1=NIRPS, 2=SPIROU)\nmin tobs [hrs]\nmin spectrograph ' + \
          '(0=HARPS, 1=NIRPS, 2=SPIROU)'
    output = np.array([starnums, Nharps, Nnirps, Nspirou, texpharps, texpnirps,
                       texpspirou, tobsharps, tobsnirps, tobsspirou, min_Nrv,
                       bestspectrograph_Nrv, min_tobs, bestspectrograph_tobs]).T
    np.savetxt('Results/median_results_3sigma_mp.dat', output, header=hdr, 
	       delimiter=',', fmt='%.4f')


def TESS_mp_5sigma():
    '''
    Compute the number of observations and total observation time required to 
    measure the TESS planets at 5 sigma (e.g. for characterizing the MR 
    relation or for habitable zone planets).
    '''
    starnums, Nharps, Nnirps, Nspirou, texpharps, texpnirps, texpspirou, tobsharps, tobsnirps, tobsspirou, min_Nrv, bestspectrograph_Nrv, min_tobs, bestspectrograph_tobs = np.loadtxt('Results/median_results_3sigma_mp.dat', delimiter=',').T

    # get fractional change in Nrv by decreasing sigmaK from the 3sigma case to
    # the 5sigma case
    # for a 5 sigma mass detection with sigP=5e-5 and sigMs/Ms = 0.1,
    # frac_sigK == .188 (.327 for 3sigma)
    frac_increase_Nrv = (.327 / .188)**2

    Nharps *= frac_increase_Nrv
    Nnirps *= frac_increase_Nrv
    Nspirou *= frac_increase_Nrv

    tobsharps = Nharps * texpharps / 60
    tobsnirps = Nnirps * texpnirps / 60
    tobsspirou = Nspirou * texpspirou / 60

    # save results to file
    hdr = 'TESS planet index\nNrv_HARPS\nNrv_NIRPS\nNrv_SPIROU\ntexp_HARPS ' + \
          '[min]\ntexp_NIRPS [min]\ntexp_SPIROU [min]\ntobs_HARPS [hrs]\n' + \
          'tobs_NIRPS [hrs]\ntobs_SPIROU [hrs]\nmin Nrv\nmin spectrograph ' + \
          '(0=HARPS, 1=NIRPS, 2=SPIROU)\nmin tobs [hrs]\nmin spectrograph ' + \
          '(0=HARPS, 1=NIRPS, 2=SPIROU)'
    output = np.array([starnums, Nharps, Nnirps, Nspirou, texpharps, texpnirps,
                       texpspirou, tobsharps, tobsnirps, tobsspirou, min_Nrv,
                       bestspectrograph_Nrv, min_tobs, bestspectrograph_tobs]).T
    np.savetxt('Results/median_results_5sigma_mp.dat', output, header=hdr, 
	       delimiter=',', fmt='%.4f')


def compute_sigrp(frac_sigRs=.1):
    '''
    Use the TESS parameters to estimate the measurement uncertainty on the 
    planet's radius. See http://adsabs.harvard.edu/abs/2008ApJ...689..499C 
    for equations.
    '''
    # Get the 3sigma results
    starnums, Nharps, Nnirps, Nspirou, texpharps, texpnirps, texpspirou, tobsharps, tobsnirps, tobsspirou, min_Nrv, bestspectrograph_Nrv, min_tobs, bestspectrograph_tobs = np.loadtxt('Results/median_results_3sigma_mp.dat', delimiter=',').T

    # Get TESS parameters including photometric uncertainty
    inds = np.array([2,3,5,6,14])
    rp, P, K, Rs, logsigV = np.ascontiguousarray(get_TESS_data())[inds]

    # Compute transit depth uncertainty
    depth    = _compute_depth(rp, Rs)
    Gamma    = _compute_Gamma()
    T        = _compute_transit_duration(rp, P, K, Rs)
    Q        = _compute_Q(Gamma, T, depth, logsigV)
    sigdepth = _compute_sigdepth(depth, Q)

    # compute corresponding planet radius uncertainty
    depth = unp.uarray(depth, sigdepth)
    Rs2 = rvs.m2Rearth(rvs.Rsun2m(unp.uarray(Rs, frac_sigRs*Rs)))
    rp2 = unp.sqrt(depth)*Rs2
    return rp, unp.std_devs(rp2)
    

def _compute_depth(rp, Rs):
    return (rvs.Rearth2m(rp) / rvs.Rsun2m(Rs))**2

def _compute_Gamma():
    cadence_minutes = 2.
    return 1./ (cadence_minutes / 60 / 24)   # days^-1

def _compute_transit_duration(rp, P, K, Rs, b=0):
    mp = np.array([runTESS.get_planet_mass(i) for i in rp])
    Ms = runTESS.get_stellar_mass(P, mp, K)
    sma = rvs.AU2m(rvs.semimajoraxis(P, Ms, mp))
    Rs2 = rvs.Rsun2m(Rs)
    tau0 = P * Rs2 / (2*np.pi*sma)
    return 2. * tau0 * np.sqrt(1. - b**2)

def _compute_Q(Gamma, T, depth, logsig):
    sig = 10**logsig
    return np.sqrt(Gamma * T) * depth / sig

def _compute_sigdepth(depth, Q):
    return depth / Q
