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
            # 'O': optical & 'I': infrared
            self.spectrographs[i] = 'O' if 'V' in self.bands[i] else 'I'
            

    def _initialize_median_arrays(self):
        self.starnums_med = np.unique(self.starnums)
	self.nstars = self.starnums_med.size
        self.ras_med, self.decs_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Ps_med, self.rps_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.mps_med, self.Ks_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Fs_med, self.Mss_med = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Rss_med, self.Teffs_med = np.zeros(self.nstars), np.zeros(self.nstars)
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


    def _initialize_median_arrays_per_spectrograph(self):
        # Optical spectrographs
        self.ras_med, self.decs_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Ps_med_O, self.rps_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.mps_med_O, self.Ks_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Fs_med_O, self.Mss_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Rss_med_O, self.Teffs_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.dists_med_O, self.Prots_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.vsinis_med_O, self.Zs_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.sigmaRV_phot_med_O = np.zeros(self.nstars)
        self.sigmaRV_acts_med_O = np.zeros(self.nstars)
        self.sigmaRV_planets_med_O = np.zeros(self.nstars)
        self.sigmaRV_eff_med_O = np.zeros(self.nstars)
        self.Rs_med_O, self.apertures_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.QEs_med_O, self.texps_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        self.fracsigmaK_targets_med_O = np.zeros(self.nstars)
        self.tobss_med_O, self.Nrvs_med_O = np.zeros(self.nstars), np.zeros(self.nstars)
        # Infrared spectrographs
        self.ras_med_I, self.decs_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Ps_med_I, self.rps_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.mps_med_I, self.Ks_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Fs_med_I, self.Mss_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.Rss_med_I, self.Teffs_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.dists_med_I, self.Prots_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.vsinis_med_I, self.Zs_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.sigmaRV_phot_med_I = np.zeros(self.nstars)
        self.sigmaRV_acts_med_I = np.zeros(self.nstars)
        self.sigmaRV_planets_med_I = np.zeros(self.nstars)
        self.sigmaRV_eff_med_I = np.zeros(self.nstars)
        self.Rs_med_I, self.apertures_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.QEs_med_I, self.texps_med_I = np.zeros(self.nstars), np.zeros(self.nstars)
        self.fracsigmaK_targets_med_I = np.zeros(self.nstars)
        self.tobss_med_I, self.Nrvs_med_I = np.zeros(self.nstars), np.zeros(self.nstars)


    def _get_median_star_results_per_spectrograph(self):
	'''
	Get the median value of each attribute for each star.
	'''
	self._initialize_median_arrays_per_spectrograph()
	
        for i in range(self.nstars):

            opt = (self.starnums == self.starnums_med[i]) & \
                  (self.spectrographs == 'O')
            self.ras_med_O[i] = np.median(self.ras[opt])
 	    self.decs_med_O[i] = np.median(self.decs[opt])
            self.Ps_med_O[i] = np.median(self.Ps[opt])
            self.rps_med_O[i] = np.median(self.rps[opt])
            self.mps_med_O[i] = np.median(self.mps[opt])
            self.Ks_med_O[i] = np.median(self.Ks[opt])
            self.Fs_med_O[i] = np.median(self.Fs[opt])
            self.Mss_med_O[i] = np.median(self.Mss[opt])
            self.Rss_med_O[i] = np.median(self.Rss[opt])
            self.Teffs_med_O[i] = np.median(self.Teffs[opt])
            self.dists_med_O[i] = np.median(self.dists[opt])
            self.Prots_med_O[i] = np.median(self.Prots[opt])
            self.vsinis_med_O[i] = np.median(self.vsinis[opt])
            self.Zs_med_O[i] = np.median(self.Zs[opt])
            self.sigmaRV_acts_med_O[i] = np.median(self.sigmaRV_acts[opt])
            self.sigmaRV_planets_med_O[i] = np.median(self.sigmaRV_planets[opt])
            self.Rs_med_O[i] = np.median(self.Rs[opt])
            self.apertures_med_O[i] = np.median(self.apertures[opt])
            self.QEs_med_O[i] = np.median(self.QEs[opt])
            self.fracsigmaK_targets_med_O[i] = np.median(self.fracsigmaK_targets[opt])
            self.sigmaRV_phot_med_O[i] = np.median(self.sigmaRV_phot[opt])
            self.sigmaRV_eff_med_O[i] = np.median(self.sigmaRV_eff[opt])
            self.texps_med_O[i] = np.median(self.texps[opt])
            self.tobss_med_O[i] = np.median(self.tobss[opt])
            self.Nrvs_med_O[i] = np.median(self.Nrvs[opt])

            nir = (self.starnums == self.starnums_med[i]) & \
                  (self.spectrographs == 'I')
            self.ras_med_I[i] = np.median(self.ras[nir])
 	    self.decs_med_I[i] = np.median(self.decs[nir])
            self.Ps_med_I[i] = np.median(self.Ps[nir])
            self.rps_med_I[i] = np.median(self.rps[nir])
            self.mps_med_I[i] = np.median(self.mps[nir])
            self.Ks_med_I[i] = np.median(self.Ks[nir])
            self.Fs_med_I[i] = np.median(self.Fs[nir])
            self.Mss_med_I[i] = np.median(self.Mss[nir])
            self.Rss_med_I[i] = np.median(self.Rss[nir])
            self.Teffs_med_I[i] = np.median(self.Teffs[nir])
            self.dists_med_I[i] = np.median(self.dists[nir])
            self.Prots_med_I[i] = np.median(self.Prots[nir])
            self.vsinis_med_I[i] = np.median(self.vsinis[nir])
            self.Zs_med_I[i] = np.median(self.Zs[nir])
            self.sigmaRV_acts_med_I[i] = np.median(self.sigmaRV_acts[nir])
            self.sigmaRV_planets_med_I[i] = np.median(self.sigmaRV_planets[nir])
            self.Rs_med_I[i] = np.median(self.Rs[nir])
            self.apertures_med_I[i] = np.median(self.apertures[nir])
            self.QEs_med_I[i] = np.median(self.QEs[nir])
            self.fracsigmaK_targets_med_I[i] = np.median(self.fracsigmaK_targets[nir])
            self.sigmaRV_phot_med_I[i] = np.median(self.sigmaRV_phot[nir])
            self.sigmaRV_eff_med_I[i] = np.median(self.sigmaRV_eff[nir])
            self.texps_med_I[i] = np.median(self.texps[nir])
            self.tobss_med_I[i] = np.median(self.tobss[nir])
            self.Nrvs_med_I[i] = np.median(self.Nrvs[nir])


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


if __name__ == '__main__':
    self = RVInformation()
