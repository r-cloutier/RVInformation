from runTESS import *


def _TEST_estimate_Nrv_TESS(mags=[15.25, 13.52],
                            band_strs=['B','V'],
                            Teff_round=3300, logg_round=5, vsini=.01, rp=1.1,
                            mp=1.6, K=2.8, R=1e5, aperture_m=3.6, QE=.1, Z=0,
                            sigmaRV_activity=3.33, testplanet_sigmaK=.33):
    startheta = mags, float(Teff_round), float(logg_round), Z, vsini
    planettheta = rp, mp, K
    instrumenttheta = band_strs, R, aperture_m, QE
    Nrv, texp, tobserving, sigmaRV_phot, sigmaRV_eff = \
                        estimate_Nrv(startheta, planettheta, instrumenttheta,
                                     sigmaRV_activity=sigmaRV_activity,
                                     testplanet_sigmaK=testplanet_sigmaK)
    print '\n%35s = %.3f m/s'%('Photon-noise limited RV uncertainty',
                               sigmaRV_phot)
    print '%35s = %.3f m/s'%('Effective RV uncertainty', sigmaRV_eff)
    print '%35s = %i'%('Number of RVs', Nrv)
    print '%35s = %.3f minutes'%('Exposure time', texp)
    print '%35s = %.3f hours'%('Total observing time', tobserving)

    return Nrv, texp, tobserving, sigmaRV_phot, sigmaRV_eff


def Test_Planet(fname):
    '''
    Given the planet parameters in a dat file in ./Nrv_tests, compute Nrv and 
    compare to the measured values.
    e.g. fname = 'Nrv_tests/GJ1132.dat'
    '''
    # Get data
    f = open('%s'%fname, 'r')
    g = f.readlines()
    f.close()
    assert len(g) == 16

    # extract parameters
    mags = list(np.ascontiguousarray(g[1].split(',')).astype(float))
    band_strs = g[2].split(',')
    band_strs[-1] = band_strs[-1].replace('\n','')
    Teff_round, logg_round, vsini = int(g[3]), float(g[4]), float(g[5])
    rp, mp, K = float(g[6]), float(g[7]), float(g[8])
    R, aperture_m, QE, Z = float(g[9]), float(g[10]), float(g[11]), float(g[12])
    sigmaRV_activity = float(g[13])
    fracsigmaK, Nrvtrue = float(g[14]), int(g[15])
    
    # Run calculator on test case
    Nrv, texp, tobserving, sigmaRV_phot, sigmaRV_eff = \
                    _TEST_estimate_Nrv_TESS(mags=mags, band_strs=band_strs,
                                            Teff_round=Teff_round,
                                            logg_round=logg_round,
                                            vsini=vsini, rp=rp, mp=mp, K=K, R=R,
                                            aperture_m=aperture_m, QE=QE, Z=Z,
                                            sigmaRV_activity=sigmaRV_activity,
                                            testplanet_sigmaK=fracsigmaK)
    print '\nNrv should be close to %i'%Nrvtrue


def Test_Planet_FAST(sigmaRV_eff, sigmaK_target):
    '''
    Quick way to make an analytical estimate of the number of RV measurements 
    required to measure K with a particular uncertainty given and effective 
    RV measurement uncertainty.  
    '''
    return int(np.round(2 * (float(sigmaRV_eff) / sigmaK_target)**2))


                     
