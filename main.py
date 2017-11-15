from runTESS import *

planetindex = int(sys.argv[1])
band_strs = [['B','V'], ['Y','J','H'], ['Y','J','H','K']]  # HARPS, NIRPS, & SPIRou
R, aperture_m = [1e5, 1e5, 7e4], [3.6, 3.6, 3.6]
nsyst_per_TESSstar = 100
assert nsyst_per_TESSstar > 0

for i in range(len(band_strs)):
    for j in range(nsyst_per_TESSstar):

        testingseed = True if nsyst_per_TESSstar == 1 else False

        print estimate_Nrv_TESS(planetindex, band_strs[i], R[i], aperture_m[i],
                                QE=.1, Z=0, sigmaRV_activity=0, sigmaRV_planets=0,
                                sigmaRV_noisefloor=noisefloor[i], testingseed=testingseed,
                                testplanet_sigmaKfrac=0, systnum=j, verbose=True)
