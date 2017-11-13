import numpy as np
import glob, os
import cPickle as pickle


class RVInformation:

    def __init__(self):

        self.files = np.array(glob.glob('Results/star*/TESS*'))
        self.nfiles = self.files.size
        self._get_data()
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
