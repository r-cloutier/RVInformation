from imports import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from RVInformationResults import *
from scipy.ndimage.filters import gaussian_filter1d

global colH, colN
colH, colN =  '#4444ff', '#a30000'

print '\nRun "plt.xkcd() in the shell'

def compute_Nrv(sigeff, sigK):
    sigeff = np.ascontiguousarray(sigeff).astype(float)
    return 2 * (sigeff / sigK)**2


def MAD(arr):
    return np.median(abs(arr-np.median(arr)))

def _compute_curve_derivative(tobs, Ndet, sigma=5):
    # avoid where the derivative is too close to 0
    g = np.append(True, np.diff(tobs) < 1e7)
    tobs, Ndet = tobs[g], Ndet[g]

    # interpolate to get uniform spacing in time
    fint = interp1d(tobs, Ndet)
    tobs2 = np.arange(tobs.min(), tobs.max(), .1)
    Ndet2 = fint(tobs2)

    # smooth curve
    Ndet_smooth = gaussian_filter1d(Ndet2, sigma=sigma)

    # take derivative
    dt = np.diff(tobs2)[0]
    dNdt = np.gradient(Ndet2, dt)
    return tobs2, dNdt


def plot_Nrvcomparsion(pltt=True, label=False):
    # get data
    starnames = np.genfromtxt('Nrv_tests/FASTtests.dat', dtype='str',
                              usecols=(0))
    _,sigeffs, sigKs, Nrv_true, Ks, rednoiseflag = \
                                    np.genfromtxt('Nrv_tests/FASTtests.dat').T
    Nrv_calc = compute_Nrv(sigeffs, sigKs)
    white = rednoiseflag == 0

    # white 
    plt.xkcd()
    fig = plt.figure(figsize=(4.5,4.5))
    ax = fig.add_subplot(111, aspect='equal')
    ax.plot([1,300],[1,300], 'k--')
    ax.scatter(Nrv_true[white], Nrv_calc[white], c=colH)

    # correlated
    names = ['CoRoT7','K218','Kep21HARPSN','Kep78HARPSN','LHS1140']
    names2 = ['CoRoT-7','K2-18','Kepler-21','Kepler-78','LHS 1140']
    xs, ys = [80,57,82,110,150], [60,85,130,140,200]
    Ks = np.array([3.42, 3.18, 2.12, 1.86, 5.3])
    sigKtargets, Nrvtargets = np.array([.66,.75,.66,.25,1.1]), \
                              np.array([71,75,82,109,144])
    sigeffs = sigKtargets * np.sqrt(Nrvtargets/2.)
    GPNrvs = np.zeros((sigeffs.size, 2))
    for i in range(sigeffs.size):
        sigKs = np.loadtxt('Nrv_tests/GPtest_%s.dat'%names[i])[:,-1]
        sigKs = sigKs[np.isfinite(sigKs)]
        Nrvs = 2 * (sigeffs[i] / sigKs)**2
        GPNrvs[i] = np.median(Nrvs), MAD(Nrvs)
    ax.scatter(Nrvtargets, GPNrvs[:,0], c=colN)

    ax.annotate('EXCELLENT AGREEMENT\nBETWEEN OBSERVED AND\nCALCULATED VALUES',
                xy=(14,14), fontsize=13, arrowprops=dict(arrowstyle='->'),
                xytext=(1.5,100))

    ax.plot([14],[4],'o', c=colH), ax.plot([14],[1.8], 'o', c=colN)
    ax.text(17, 4, "known transiting planets\nwith residual RV noise\ntreated as white", verticalalignment='center', fontsize=9)
    ax.text(17, 1.8, "known transiting planets\nwith correlated residual\nRV noise", verticalalignment='center', fontsize=9)
    
    ax.set_xscale('log'), ax.set_yscale('log')
    ax.set_xlabel('OBSERVED N$_{RV}$', fontsize=14)
    ax.set_ylabel('CALCULATED N$_{RV}$\nBASED ON OUR MODEL', fontsize=14)
    ax.set_xticks(np.logspace(0,2,3))
    ax.set_xticklabels(['%i'%i for i in np.logspace(0,2,3)], fontsize=16)
    ax.set_yticks(np.logspace(0,2,3))
    ax.set_yticklabels(['%i'%i for i in np.logspace(0,2,3)], fontsize=16)
    ax.set_xlim((1,300)), ax.set_ylim((1,300))

    fig.subplots_adjust(bottom=.18, left=.16, right=.99, top=.97)
    if label:
        plt.savefig('/Users/ryancloutier/Talks/ExoplanetsII/Nrvcomparison_poster.png')
    if pltt:
        plt.show()
    plt.close('all')


    
def plot_cumulativetobsGP_50(self, pltt=True, label=False):
    # 5 sigma correction
    corr = (.327/.189)**2 
    g = self.rps_med < 4
    tobsH = self.tobsGPs_med_H / 7.
    tobsN = self.tobsGPs_med_N / 7.

    fig = plt.figure(figsize=(4.2,5.3))
    gs = gridspec.GridSpec(8,1)
    ax1 = plt.subplot(gs[:5,:])
    ax2 = plt.subplot(gs[5:,:])
    
    tobsN = np.sort(tobsN[g] * corr)
    tobsN = np.append(0, np.cumsum(tobsN))[:401]
    Ndet = np.arange(tobsN.size)
    ax1.plot(tobsN, Ndet, '-', c=colN, lw=2, drawstyle='steps')
    print tobsN[Ndet==20] 

    inds = np.arange(g.sum())
    np.random.seed(0)
    for i in range(10):
        np.random.shuffle(inds)
        tobs = self.tobsGPs_med_H[g][inds] * corr/7.
        tobs = np.append(0, np.cumsum(tobs))[:401]
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
        tobs = self.tobsGPs_med_N[g][inds] * corr/7.
        tobs = np.append(0, np.cumsum(tobs))[:401]
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)

    tobsH = np.sort(tobsH[g] * corr)
    tobsH = np.append(0, np.cumsum(tobsH))[:401]
    ax1.plot(tobsH, Ndet, '--', c=colH, lw=2, drawstyle='steps')

    ax2.plot([10,18], np.repeat(1.4,2), '--', c=colH)
    ax2.text(20, 1.4, 'OPTICAL SPECTROGRAPH\n(BVR bands, R=100,000)',
             color=colH, fontsize=10, verticalalignment='center')
    ax2.plot([10,18], np.repeat(2.6,2), '-', c=colN)
    ax2.text(20, 2.6, 'NEAR-IR SPECTROGRAPH\n(YJH bands, R=100,000)',
             color=colN, fontsize=10, verticalalignment='center')
    
    tobs2, dNdt = _compute_curve_derivative(tobsH*7, Ndet)
    ax2.plot(tobs2, 1./dNdt, '--', c=colH, lw=2)

    tobs2, dNdt = _compute_curve_derivative(tobsN*7, Ndet)
    ax2.plot(tobs2, 1./dNdt, '-', c=colN, lw=2)

    ax1.annotate('50 SMALL PLANET \nMASSES AT 5$\sigma$ IN\nAS LITTLE AS ~65\nNIGHTS',
                xy=(65,50), fontsize=14, arrowprops=dict(arrowstyle='->'),
                xytext=(4,210))
    ax1.annotate('SORTED TESS PLANETS\n(I.E. THE MOST\nEFFICIENTLY OBSERVED\nTESS PLANETS)',
                xy=(22,18), fontsize=10, arrowprops=dict(arrowstyle='->'),
                xytext=(1.4,58))
    ax1.annotate('RANDOM\nTESS\nPLANETS\nREQUIRE\nMORE\nOBSERVING\nTIME',
                xy=(220,8), fontsize=10, arrowprops=dict(arrowstyle='->'),
                xytext=(243,60))
    ax2.annotate('CAN DETECT SINGLE\nPLANET MASSES IN\nAS LITTLE AS ~7 HOURS\nAFTER 65 NIGHTS',
                xy=(65,7), fontsize=10, arrowprops=dict(arrowstyle='->'),
                xytext=(1.5,25))
    
    ax1.set_xlim((1,1e3)), ax2.set_xlim((1,1e3))
    ax1.set_ylim((-5,3e2))
    ax1.set_xscale('log'), ax2.set_xscale('log')
    ax2.set_yscale('log'), ax2.set_ylim((30,.8))
    ax1.set_xticklabels('', fontsize=1)
    ax2.set_yticks([1,10]), ax2.set_yticklabels(['1','10'])
    ax2.set_xticks(np.logspace(0,3,4))
    ax2.set_xticklabels(['%i'%i for i in np.logspace(0,3,4)])
    ax2.set_xlabel('CUMULATIVE OBSERVING NIGHTS')
    ax2.set_ylabel('NUMBER OF HOURS\nPER PLANET DETECTION')
    ax1.set_ylabel('CUMULATIVE NUMBER\nOF PLANET DETECTIONS')
    
    fig.subplots_adjust(left=.19, right=.95, top=.98, bottom=.08)
    if label:
        plt.savefig('/Users/ryancloutier/Talks/ExoplanetsII/cumulativetobs50_poster.png')
    if pltt:
        plt.show()
    plt.close('all')
