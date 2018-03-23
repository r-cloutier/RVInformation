from imports import *

global h, c
h, c = 6.62607004e-34, 3e8

def compute_transmission_SNR(self):
    # compute noise
    throughput = .5
    binwidth = 15e-3  # microns
    D, wlcen, F0 = 25e4, 1.25, 2.985e-6  # cm2, microns, ergs/s/cm2/micron
    F = F0 * 10**(-.4*self.Jmags_med)  # ergs/s/cm2/micron
    T14 = rvs.days2sec(rvs.transit_width(self.Ps_med, self.Mss_med, self.Rss_med, self.rps_med, 0, self.mps_med)) # sec
    Ephot = h*c / (wlcen*1e-6) * 1e7  # ergs
    Nphot = F * T14 * D * binwidth / Ephot * throughput
    noise_ppm = 1e6 / np.sqrt(Nphot)   # X ppm = 1e6/SNR, so SNR=1e6 is 1 ppm

    # compute SNR of spectral features
    SNR = self.transmission_ppm / noise_ppm

    return noise, SNR


def compute_transmission_SNR_single(Jmag, P, Ms, Rs, rp, mp, Tp, mu, b=0, binwidth=15):
    # compute noise
    throughput = .5
    binwidth = 15e-3  # microns
    D, wlcen, F0 = 25e4, 1.25, 2.985e-6  # cm2, microns, ergs/s/cm2/micron
    F = F0 * 10**(-.4*Jmag)  # ergs/s/cm2/micron
    T14 = rvs.days2sec(rvs.transit_width(P, Ms, Rs, rp, b, mp)) # sec
    Ephot = h*c / (wlcen*1e-6) * 1e7  # ergs
    Nphot = F * T14 * D * binwidth / Ephot * throughput
    noise_ppm = 1e6 / np.sqrt(Nphot)   # X ppm = 1e6/SNR, so SNR=1e6 is 1 ppm

    # compute SNR of spectral features
    transmission_ppm = rvs.transmission_spectroscopy_depth(Rs, mp, rp, Tp, mu)
    SNR = transmission_ppm / noise_ppm

    return SNR
