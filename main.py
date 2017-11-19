from runTESS import *

planetindex = int(sys.argv[1])
band_strs = [['B','V'], ['Y','J','H'], ['Y','J','H','K']]  # HARPS, NIRPS, & SPIRou
R, aperture_m, noisefloor = [1e5, 1e5, 7e4], [3.6, 3.6, 3.6], [.5, 1., 1.]
nsyst_per_TESSstar = 1
assert nsyst_per_TESSstar > 0

for i in range(len(band_strs)):
    for j in range(nsyst_per_TESSstar):

        testingseed = True if nsyst_per_TESSstar == 1 else False

        if do_i_run_this_job(planetindex, band_strs[i], j):

            estimate_Nrv_TESS(planetindex, band_strs[i], R[i], aperture_m[i],
                              QE=.1, Z=0, sigmaRV_activity=0, sigmaRV_planets=0,
                              sigmaRV_noisefloor=noisefloor[i], testingseed=testingseed,
                              testplanet_sigmaKfrac=0, systnum=j, verbose=False)



def do_i_run_this_job(planetindex, band_strs, systnum):
    fname = 'Results/star%.4d/TESSplanet%.4d_%s_%.4d'%(planetindex,
                                                       planetindex,
                                                       ''.join(band_strs),
                                                       systnum)))
    return not os.path.exists(fname)
