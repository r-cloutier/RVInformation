from RVInformationResults import *


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
    np.savetxt('Results/median_results.dat', output, header=hdr, delimiter=',',
               fmt='%.4f')
