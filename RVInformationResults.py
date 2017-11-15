import numpy as np
import glob, os
import cPickle as pickle


class RVInformation:

    def __init__(self):

        self.files = np.array(glob.glob('Results/star*/TESS*'))
        self.nfiles = self.files.size
        self._get_data()
	self._get_median_star_results()
        self._get_median_star_results_per_spectrograph()
	self._pickleobject()


    def _initialize_arrays(self):
        self.starnums, self.systnums = np.zeros(self.nfiles), \
                                       np.zeros(self.nfiles)
        self.spectrographs = np.zeros(self.nfiles, dtype=str)
        self.ras, self.decs = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.Ps, self.rps = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.mps, self.Ks = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.Fs, self.Mss = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.Rss, self.Teffs = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.Bmags, self.Vmags = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.Ymags, self.Jmags = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.Hmags, self.Kmags = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.dists, self.Prots = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.vsinis, self.Zs = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.sigmaRV_phot = np.zeros(self.nfiles)
        self.sigmaRV_acts = np.zeros(self.nfiles)
        self.sigmaRV_planets = np.zeros(self.nfiles)
        self.sigmaRV_eff = np.zeros(self.nfiles)
        self.Rs, self.apertures = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.QEs, self.texps = np.zeros(self.nfiles), np.zeros(self.nfiles)
        self.fracsigmaK_targets = np.zeros(self.nfiles)
        self.tobss, self.Nrvs = np.zeros(self.nfiles), np.zeros(self.nfiles)

        
    def _get_data(self):
        '''
        Compile all the saved data for each simulated system.
        '''
        self._initialize_arrays()

        self.bands, self.mags = [], []
        for i in range(self.nfiles):
            # Get saved properties aside the magnitudes
            self.ras[i], self.decs[i], self.Ps[i], self.rps[i], self.mps[i], self.Ks[i], self.Fs[i], self.Mss[i], self.Rss[i], self.Teffs[i], self.dists[i], self.Prots[i], self.vsinis[i], self.Zs[i], self.sigmaRV_acts[i], self.sigmaRV_planets[i], self.Rs[i], self.apertures[i], self.QEs[i], self.fracsigmaK_targets[i], self.sigmaRV_phot[i], self.sigmaRV_eff[i], self.texps[i], self.tobss[i], self.Nrvs[i] = np.loadtxt(self.files[i])
            self.starnums[i] = int(self.files[i].split('star')[-1].split('/')[0])
            self.systnums[i] = int(self.files[i].split('_')[-1].split('.')[0])

            # Get bands and magnitudes
            f = open(self.files[i], 'r')
            g = f.readlines()
            f.close()
            bands, mags = [], []
            for j in range(len(g)-2):
                bands.append(g[j].split(' ')[1])
                mags.append(float(g[j].split(' ')[3]))
            self.bands.append(bands)
            self.mags.append(mags)

            # get mags
            if 'V' in self.bands[i]:
                self.spectrographs[i] = 'H'
                self.Bmags[i], self.Vmags[i] = self.mags[i][0], self.mags[i][1]
            elif 'K' in self.bands[i]:
                self.spectrographs[i] = 'S'
                self.Kmags[i] = self.mags[i][-1]
            elif 'J' in self.bands[i]:
                self.spectrographs[i] = 'N'
                self.Ymags[i], self.Jmags[i], self.Hmags[i] = self.mags[i][0], \
                                                              self.mags[i][1], \
                                                              self.mags[i][2]
            else:
                raise ValueError("Spectra band do not contain a reference band.")
                    
                

    def _initialize_median_arrays(self):
        self.starnums_med = np.unique(self.starnums)
	self.nstars = self.starnums_med.size
        self.ras_med, self.decs_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Ps_med, self.rps_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.mps_med, self.Ks_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Fs_med, self.Mss_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Rss_med, self.Teffs_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Bmags_med, self.Vmags_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Ymags_med, self.Jmags_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Hmags_med, self.Kmags_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.dists_med, self.Prots_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.vsinis_med, self.Zs_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.sigmaRV_phot_med = np.zeros(self.nstars)
        self.sigmaRV_acts_med = np.zeros(self.nstars)
        self.sigmaRV_planets_med = np.zeros(self.nstars)
        self.sigmaRV_eff_med = np.zeros(self.nstars)
        self.Rs_med, self.apertures_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.QEs_med, self.texps_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.fracsigmaK_targets_med = np.zeros(self.nstars)
        self.tobss_med, self.Nrvs_med = np.zeros(self.nstars), np.zeros(self.nstars)
        # Errors
        self.Prots_emed = np.zeros(self.nstars)
        self.vsinis_emed = np.zeros(self.nstars)
        self.sigmaRV_phot_emed = np.zeros(self.nstars)
        self.sigmaRV_acts_emed = np.zeros(self.nstars)
        self.sigmaRV_planets_emed = np.zeros(self.nstars)
        self.sigmaRV_eff_emed = np.zeros(self.nstars)
        self.texps_emed = np.zeros(self.nstars)
        self.fracsigmaK_targets_emed = np.zeros(self.nstars)
        self.tobss_emed, self.Nrvs_emed = np.zeros(self.nstars), np.zeros(self.nstars)


    def _get_median_star_results(self):
	'''
	Get the median value of each attribute for each star.
	'''
	self._initialize_median_arrays()
	
        for i in range(self.nstars):
	    g = self.starnums == self.starnums_med[i]
	    self.ras_med[i] = np.median(self.ras[g])
 	    self.decs_med[i] = np.median(self.decs[g])
            self.Ps_med[i] = np.median(self.Ps[g])
            self.rps_med[i] = np.median(self.rps[g])
            self.mps_med[i] = np.median(self.mps[g])
            self.Ks_med[i] = np.median(self.Ks[g])
            self.Fs_med[i] = np.median(self.Fs[g])
            self.Mss_med[i] = np.median(self.Mss[g])
            self.Rss_med[i] = np.median(self.Rss[g])
            self.Teffs_med[i] = np.median(self.Teffs[g])
            self.Bmags_med[i] = np.median(self.Bmags[(g) & (self.Bmags != 0)])
            self.Vmags_med[i] = np.median(self.Vmags[(g) & (self.Vmags != 0)])
            self.Ymags_med[i] = np.median(self.Ymags[(g) & (self.Ymags != 0)])
            self.Jmags_med[i] = np.median(self.Jmags[(g) & (self.Jmags != 0)])
            self.Hmags_med[i] = np.median(self.Hmags[(g) & (self.Hmags != 0)])
            self.Kmags_med[i] = np.median(self.Kmags[(g) & (self.Kmags != 0)])
            self.dists_med[i] = np.median(self.dists[g])
            self.Prots_med[i] = np.median(self.Prots[g])
            self.vsinis_med[i] = np.median(self.vsinis[g])
            self.Zs_med[i] = np.median(self.Zs[g])
            self.sigmaRV_acts_med[i] = np.median(self.sigmaRV_acts[g])
            self.sigmaRV_planets_med[i] = np.median(self.sigmaRV_planets[g])
            self.Rs_med[i] = np.median(self.Rs[g])
            self.apertures_med[i] = np.median(self.apertures[g])
            self.QEs_med[i] = np.median(self.QEs[g])
            self.fracsigmaK_targets_med[i] = np.median(self.fracsigmaK_targets[g])
            self.sigmaRV_phot_med[i] = np.median(self.sigmaRV_phot[g])
            self.sigmaRV_eff_med[i] = np.median(self.sigmaRV_eff[g])
            self.texps_med[i] = np.median(self.texps[g])
            self.tobss_med[i] = np.median(self.tobss[g])
            self.Nrvs_med[i] = np.median(self.Nrvs[g])
            # Errors
            self.Prots_emed[i] = MAD(self.Prots[g])
            self.vsinis_emed[i] = MAD(self.vsinis[g])
            self.sigmaRV_acts_emed[i] = MAD(self.sigmaRV_acts[g])
            self.sigmaRV_planets_emed[i] = MAD(self.sigmaRV_planets[g])
            self.fracsigmaK_targets_emed[i] = MAD(self.fracsigmaK_targets[g])
            self.sigmaRV_phot_emed[i] = MAD(self.sigmaRV_phot[g])
            self.sigmaRV_eff_emed[i] = MAD(self.sigmaRV_eff[g])
            self.texps_emed[i] = MAD(self.texps[g])
            self.tobss_emed[i] = MAD(self.tobss[g])
            self.Nrvs_emed[i] = MAD(self.Nrvs[g])


    def _initialize_median_arrays_per_spectrograph(self):
        # Optical spectrographs
        self.Prots_med_O = np.zeros(self.nstars)
        self.vsinis_med_O = np.zeros(self.nstars)
        self.sigmaRV_phot_med_O = np.zeros(self.nstars)
        self.sigmaRV_acts_med_O = np.zeros(self.nstars)
        self.sigmaRV_planets_med_O = np.zeros(self.nstars)
        self.sigmaRV_eff_med_O = np.zeros(self.nstars)
        self.Rs_med_O, self.apertures_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.texps_med_O = np.zeros(self.nstars)
        self.fracsigmaK_targets_med_O = np.zeros(self.nstars)
        self.tobss_med_O, self.Nrvs_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        # Errors
        self.Prots_emed_O = np.zeros(self.nstars)
        self.vsinis_emed_O = np.zeros(self.nstars)
        self.sigmaRV_phot_emed_O = np.zeros(self.nstars)
        self.sigmaRV_acts_emed_O = np.zeros(self.nstars)
        self.sigmaRV_planets_emed_O = np.zeros(self.nstars)
        self.sigmaRV_eff_emed_O = np.zeros(self.nstars)
        self.texps_emed_O = np.zeros(self.nstars)
        self.fracsigmaK_targets_emed_O = np.zeros(self.nstars)
        self.tobss_emed_O, self.Nrvs_emed_O = np.zeros(self.nstars), np.zeros(self.nstars)
        
        # Infrared spectrographs
        self.Prots_med_I = np.zeros(self.nstars)
        self.vsinis_med_I = np.zeros(self.nstars)
        self.sigmaRV_phot_med_I = np.zeros(self.nstars)
        self.sigmaRV_acts_med_I = np.zeros(self.nstars)
        self.sigmaRV_planets_med_I = np.zeros(self.nstars)
        self.sigmaRV_eff_med_I = np.zeros(self.nstars)
        self.Rs_med_I, self.apertures_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.QEs_med_I, self.texps_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.fracsigmaK_targets_med_I = np.zeros(self.nstars)
        self.tobss_med_I, self.Nrvs_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        # Errors
        self.Prots_emed_I = np.zeros(self.nstars)
        self.vsinis_emed_I = np.zeros(self.nstars)
        self.sigmaRV_phot_emed_I = np.zeros(self.nstars)
        self.sigmaRV_acts_emed_I = np.zeros(self.nstars)
        self.sigmaRV_planets_emed_I = np.zeros(self.nstars)
        self.sigmaRV_eff_emed_I = np.zeros(self.nstars)
        self.texps_emed_I = np.zeros(self.nstars)
        self.fracsigmaK_targets_emed_I = np.zeros(self.nstars)
        self.tobss_emed_I, self.Nrvs_emed_I = np.zeros(self.nstars), np.zeros(self.nstars)


    def _get_median_star_results_per_spectrograph(self):
	'''
	Get the median value of each attribute for each star.
	'''
	self._initialize_median_arrays_per_spectrograph()
	
        for i in range(self.nstars):

            opt = (self.starnums == self.starnums_med[i]) & \
                  (self.spectrographs == 'H')
            self.Prots_med_O[i] = np.median(self.Prots[opt])
            self.vsinis_med_O[i] = np.median(self.vsinis[opt])
            self.sigmaRV_acts_med_O[i] = np.median(self.sigmaRV_acts[opt])
            self.sigmaRV_planets_med_O[i] = np.median(self.sigmaRV_planets[opt])
            self.fracsigmaK_targets_med_O[i] = np.median(self.fracsigmaK_targets[opt])
            self.sigmaRV_phot_med_O[i] = np.median(self.sigmaRV_phot[opt])
            self.sigmaRV_eff_med_O[i] = np.median(self.sigmaRV_eff[opt])
            self.texps_med_O[i] = np.median(self.texps[opt])
            self.tobss_med_O[i] = np.median(self.tobss[opt])
            self.Nrvs_med_O[i] = np.median(self.Nrvs[opt])
            # Errors
            self.Prots_emed_O[i] = MAD(self.Prots[opt])
            self.vsinis_emed_O[i] = MAD(self.vsinis[opt])
            self.sigmaRV_acts_emed_O[i] = MAD(self.sigmaRV_acts[opt])
            self.sigmaRV_planets_emed_O[i] = MAD(self.sigmaRV_planets[opt])
            self.fracsigmaK_targets_emed_O[i] = MAD(self.fracsigmaK_targets[opt])
            self.sigmaRV_phot_emed_O[i] = MAD(self.sigmaRV_phot[opt])
            self.sigmaRV_eff_emed_O[i] = MAD(self.sigmaRV_eff[opt])
            self.texps_emed_O[i] = MAD(self.texps[opt])
            self.tobss_emed_O[i] = MAD(self.tobss[opt])
            self.Nrvs_emed_O[i] = MAD(self.Nrvs[opt])

            nir = (self.starnums == self.starnums_med[i]) & \
                  (self.spectrographs == 'N')
            self.Prots_med_I[i] = np.median(self.Prots[nir])
            self.vsinis_med_I[i] = np.median(self.vsinis[nir])
            self.sigmaRV_acts_med_I[i] = np.median(self.sigmaRV_acts[nir])
            self.sigmaRV_planets_med_I[i] = np.median(self.sigmaRV_planets[nir])
            self.fracsigmaK_targets_med_I[i] = np.median(self.fracsigmaK_targets[nir])
            self.sigmaRV_phot_med_I[i] = np.median(self.sigmaRV_phot[nir])
            self.sigmaRV_eff_med_I[i] = np.median(self.sigmaRV_eff[nir])
            self.texps_med_I[i] = np.median(self.texps[nir])
            self.tobss_med_I[i] = np.median(self.tobss[nir])
            self.Nrvs_med_I[i] = np.median(self.Nrvs[nir])
            # Errors
            self.Prots_emed_I[i] = MAD(self.Prots[nir])
            self.vsinis_emed_I[i] = MAD(self.vsinis[nir])
            self.sigmaRV_acts_emed_I[i] = MAD(self.sigmaRV_acts[nir])
            self.sigmaRV_planets_emed_I[i] = MAD(self.sigmaRV_planets[nir])
            self.fracsigmaK_targets_emed_I[i] = MAD(self.fracsigmaK_targets[nir])
            self.sigmaRV_phot_emed_I[i] = MAD(self.sigmaRV_phot[nir])
            self.sigmaRV_eff_emed_I[i] = MAD(self.sigmaRV_eff[nir])
            self.texps_emed_I[i] = MAD(self.texps[nir])
            self.tobss_emed_I[i] = MAD(self.tobss[nir])
            self.Nrvs_emed_I[i] = MAD(self.Nrvs[nir])


    def _pickleobject(self):
	try:
	    os.mkdir('pickles')
	except OSError:
	    pass
	f = open('pickles/RVInformation', 'wb')
	pickle.dump(self, f)
	f.close()


def loadpickle(fname):
    f = open(fname, 'rb')
    self = pickle.load(f)
    f.close()
    return self

def MAD(arr):
    return np.median(abs(arr-np.median(arr)))


if __name__ == '__main__':
    self = RVInformation()
