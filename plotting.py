import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from uncertainties import unumpy as unp
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import matplotlib as mpl
import matplotlib.colors as colors
from CVZ import *
from JWST_SNR import *
from get_tess_data import get_TESS_data
from analysis import *
from rvmodel import get_rv1
from compute_sigK_QPGP import compute_sigmaK_GP
from string import ascii_uppercase

#params = {'text.usetex': True, 
#          'text.latex.preamble': [r'\usepackage{cmbright}',
#                                  r'\usepackage{amsmath}']}
#plt.rcParams.update(params)

global colH, colN, sigmaSB, kB, G, mp
colH, colN, sigmaSB, kB, G = '#4444ff', '#a30000', 5.67e-8, 1.38e-23, 6.67e-11
mp = 1.67e-27

mpl.rc('xtick', labelsize=13)
mpl.rc('ytick', labelsize=13)
mpl.rc('axes', titlepad=3)

global toverhead, Teffedges, Tefflabels, rpedges, rplabels, rplabels2
toverhead = 0.
Teffedges = np.array([25e2, 32e2, 38e2, 76e2, 12e3])
Tefflabels = ['mid-late M', 'early-mid M', 'FGK', 'BA']
rpedges = np.array([0, 1.25, 2, 4, 30])
rplabels = ['Earths','Super-Earths','Neptunes','Giants']
rplabels2 = ['$<$1.25 R$_{\oplus}$','1.25-2 R$_{\oplus}$',
                 '2-4 R$_{\oplus}$','$>$4 R$_{\oplus}$']
    
def add_RVnoise_floor(self, RVnoisefloor_N, RVnoisefloor_H):
    answer = raw_input('Are you sure that you want to add this noise floor? ')
    if answer.lower() in ['y','yes','ye']:
        sigeffN_2 = (self.sigmaRV_eff_med_N**2 + RVnoisefloor_N**2) / \
                    self.sigmaRV_eff_med_N**2
        self.NrvGPs_med_N *= sigeffN_2
        self.tobsGPs_med_N *= sigeffN_2

        sigeffH_2 = (self.sigmaRV_eff_med_H**2 + RVnoisefloor_H**2) / \
                    self.sigmaRV_eff_med_H**2
        self.NrvGPs_med_H *= sigeffH_2
        self.tobsGPs_med_H *= sigeffH_2
        
    

def _truncate_colormap(cmap, minval=0, maxval=1, n=100):
    return colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f},)'.format(n=cmap.name,a=minval,b=maxval),
        cmap(np.linspace(minval, maxval, n)))


def compute_Nrv(sigeff, sigK):
    sigeff = np.ascontiguousarray(sigeff).astype(float)
    return 2 * (sigeff / sigK)**2


def get_planets(tobs, g, N=0, sort=True):
    tobsinds = np.argsort(tobs[g]) if sort else np.arange(tobs[g].size)
    inds = np.arange(tobs.size)[g][tobsinds]
    return inds if N == 0 else inds[:int(N)]


def plot_Nrv_hists(self, nbins=40, pltt=True, label=False):
    bins = np.logspace(0,5,nbins)
    opt, nir = self.spectrographs == 'H', self.spectrographs == 'N'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(self.Nrvs[opt], bins=bins, color='b', alpha=.6, normed=0,
            stacked=1, label='Optical')
    ax.hist(self.Nrvs[nir], bins=bins, color='r', alpha=.6, normed=0,
            stacked=1, label='nIR')
    ax.set_xscale('log')
    ax.set_xlabel('$n_{RV}$')#, ax.set_ylabel('

    if pltt:
	plt.show()
    plt.close('all')


def plot_tobs_hists_v_Teff(self, nbins=40, pltt=True, label=False):
    bins = np.logspace(0,5,nbins)
    opt, nir = self.spectrographs == 'H', self.spectrographs == 'N'
    fig = plt.figure(figsize=(12,3.4))
    for i in range(len(Teffedges)-1): 
        ax = fig.add_subplot(1,len(Teffedges)-1,-i+4)
	gteff = (self.Teffs >= Teffedges[i]) & (self.Teffs <= Teffedges[i+1])
    	opt = (self.spectrographs == 'H') & gteff
	ax.hist(self.tobss[opt], bins=bins, color='b', alpha=.6,
                label='Optical')
	ax.text(.68, .85, 'med(tobs)=%.1f'%(np.median(self.tobss[opt])),
                color='b', transform=ax.transAxes, fontsize=10)
	nir = (self.spectrographs == 'N') & gteff
        ax.hist(self.tobss[nir], bins=bins, color='r', alpha=.6,
                label='nIR')
        ax.text(.68, .8, 'med(tobs)=%.1f'%(np.median(self.tobss[nir])),
                color='r', transform=ax.transAxes, fontsize=10)
    	ax.set_xscale('log'), ax.set_xlabel('$t_{obs}$ [hours]', fontsize=10)
	ax.set_yticklabels('')
	ax.set_title('%s\n%i $\leq T_{eff} <$ %i'%(Tefflabels[i],
                                                   Teffedges[i],
                                                   Teffedges[i+1]),
                     fontsize=10)
        ax.minorticks_on()

    fig.subplots_adjust(left=.04, bottom=.17, top=.87, right=.96)
    if label:
        plt.savefig('plots/tobshist_Teff.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_tobs_hists_v_rp(self, nbins=40, pltt=True, label=False):
    bins = np.logspace(0,5,nbins)
    opt, nir = self.spectrographs == 'H', self.spectrographs == 'N'
    fig = plt.figure(figsize=(12,3.4))
    for i in range(len(rpedges)-1):
        ax = fig.add_subplot(1,len(rpedges)-1,i+1)
        grp = (self.rps >= rpedges[i]) & (self.rps <= rpedges[i+1])
        opt = (self.spectrographs == 'H') & grp
        ax.hist(self.tobss[opt], bins=bins, color='b', alpha=.6,
                label='Optical')
        ax.text(.65, .85, 'med(tobs)=%.1f'%(np.median(self.tobss[opt])),
                color='b', transform=ax.transAxes, fontsize=10)
        nir = (self.spectrographs == 'N') & grp
        ax.hist(self.tobss[nir], bins=bins, color='r', alpha=.6,
                label='nIR')
        ax.text(.65, .8, 'med(tobs)=%.1f'%(np.median(self.tobss[nir])),
                color='r', transform=ax.transAxes, fontsize=10)
        ax.set_xscale('log'), ax.set_xlabel('$t_{obs}$ [hours]', fontsize=10)
        ax.set_yticklabels('')
        ax.set_title('%s\n%s'%(rplabels[i], rplabels2[i]), fontsize=10)

    fig.subplots_adjust(left=.04, bottom=.17, top=.87, right=.96)
    if label:
        plt.savefig('plots/tobshist_rp.png')
    if pltt:
        plt.show()
    plt.close('all')


# bad ones: WASP-19, XO-2, WASP-22
def plot_Nrv_estimates(pltt=True, label=False):
    # Get real planet data
    starnames = np.genfromtxt('Nrv_tests/FASTtests.dat', dtype='str',
			      usecols=(0))
    _,sigeffs, sigKs, Nrv_true, Ks, rednoiseflag = \
                                np.genfromtxt('Nrv_tests/FASTtests.dat').T
    detsig = Ks / sigKs
    c = Ks
    
    # Compute Nrv
    Nrv_calc = compute_Nrv(sigeffs, sigKs)

    for i in range(sigeffs.size):
	print starnames[i], sigeffs[i], sigKs[i], Nrv_true[i], Nrv_calc[i], '\n'

    # plotting
    fig = plt.figure(figsize=(3.8,4.3))
    ax = fig.add_subplot(111, aspect='equal')
    vmin, vmax = 1, 3e2
    g = rednoiseflag == 0
    #np.random.seed(1)
    #np.random.shuffle(g)
    cax = ax.scatter(Nrv_true[g], Nrv_calc[g], edgecolors='none', marker='o',
                     c=c[g], s=60, label='white noise model',
                     norm=LogNorm(vmin=vmin, vmax=vmax),
		     cmap=_truncate_colormap(plt.get_cmap('rainbow'),0,1))
    cbar_axes = fig.add_axes([.08,.1,.84,.04])
    cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')
    cticklabels = ['']*22
    cticklabels[0], cticklabels[9], cticklabels[18] = '1', '10', '100'
    cbar.ax.set_xticklabels(cticklabels) 
    cbar.set_label(r'K [m/s]', labelpad=.1)
   
    # plot bad red nosie sources
    plt_red = 0
    if plt_red:
	g = rednoiseflag == 1
	print 'here %i'%(g.sum())
	ax.scatter(Nrv_true[g], Nrv_calc[g], edgecolors='k', marker='d',
                   c=c[g], s=60, label='correlated noise model',
		   norm=LogNorm(vmin=vmin, vmax=vmax),
                   cmap=_truncate_colormap(plt.get_cmap('rainbow'),0,1)) 

    names = ['CoRoT7','K218','Kep21HARPSN','Kep78HARPSN','LHS1140']
    names2 = ['CoRoT-7','K2-18','Kepler-21','Kepler-78','LHS 1140']
    xs, ys = [80,57,82,110,270], [60,85,130,140,180]
    Ks = np.array([3.42, 3.18, 2.12, 1.86, 4.71])
    sigKtargets, Nrvtargets = np.array([.66,.75,.66,.25,.5]), \
			      np.array([71,75,82,109,293])
    sigeffs = sigKtargets * np.sqrt(Nrvtargets/2.)
    GPNrvs = np.zeros((sigeffs.size, 2))
    for i in range(sigeffs.size):
        sigKs = np.loadtxt('Nrv_tests/GPtest_%s.dat'%names[i])[:,-1]
        sigKs = sigKs[np.isfinite(sigKs)]
        Nrvs = 2 * (sigeffs[i] / sigKs)**2
        GPNrvs[i] = np.median(Nrvs), MAD(Nrvs)
        ax.text(xs[i], ys[i], ascii_uppercase[i], fontsize=8)
        ax.text(30, np.logspace(np.log10(7),np.log10(2),sigeffs.size)[i],
                '%s = %s'%(ascii_uppercase[i], names2[i]), fontsize=8)
    ax.errorbar(Nrvtargets, GPNrvs[:,0], GPNrvs[:,1], fmt='ko', ms=0,
                elinewidth=.8, capsize=0)
    ax.scatter(Nrvtargets, GPNrvs[:,0], edgecolors='k', marker='d',
               c=Ks/sigKtargets, s=60, label='correlated noise model',
               norm=LogNorm(vmin=vmin, vmax=vmax),
	       cmap=_truncate_colormap(plt.get_cmap('rainbow'),0,1))
    
    arr = np.linspace(0, 350, 100)
    ax.plot(arr, arr, 'k--', lw=2)

    ax.set_xscale('log'), ax.set_yscale('log')
    ax.set_xlabel('Observed N$_{RV}$'), ax.set_ylabel('Calculated N$_{RV}$')
    ax.set_xlim((1, 350)), ax.set_ylim((1, 350))
    ax.set_xticks(np.logspace(0,2,3))
    ax.set_xticklabels(['%i'%i for i in np.logspace(0,2,3)])
    ax.set_yticks(np.logspace(0,2,3))
    ax.set_yticklabels(['%i'%i for i in np.logspace(0,2,3)])
    ax.minorticks_on()

    ax.legend(bbox_to_anchor=(.67,.95), fontsize=9,
              handletextpad=.01, labelspacing=.7)
    fig.subplots_adjust(bottom=.27, left=.19, right=.97, top=1)
    if label:
        plt.savefig('plots/Nrvcomparison.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_Nrv_estimates_4referee(pltt=True, label=False):
    # Get real planet data
    starnames = np.genfromtxt('Nrv_tests/FASTtests.dat', dtype='str',
                              usecols=(0))
    _,sigeffs, sigKs, Nrv_true, Ks, rednoiseflag = \
                                np.genfromtxt('Nrv_tests/FASTtests.dat').T
    detsig = Ks / sigKs
    c = Ks

    # Compute Nrv
    Nrv_calc = compute_Nrv(sigeffs, sigKs)
    g = rednoiseflag == 0

    plt.scatter(Nrv_true[g], Nrv_calc[g], c=c[g],
		cmap=_truncate_colormap(plt.get_cmap('rainbow'),0,1))
    #plt.colorbar()
    plt.plot([1,250], [1,250], 'k--')
    plt.xlabel('Observed N_RV'), plt.ylabel('Calculated N_RV')

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
    plt.errorbar(Nrvtargets, GPNrvs[:,0], GPNrvs[:,1], fmt='kd', ms=5,
                 elinewidth=.8, capsize=0)

    Nrv_obs = np.append(Nrv_true[g], Nrvtargets)
    Nrv_cal = np.append(Nrv_calc[g], GPNrvs[:,0])
    print 'rms white = %.3f'%(np.std(Nrv_true[g]-Nrv_calc[g]))
    print 'rms correlated = %.3f'%(np.std(Nrvtargets-GPNrvs[:,0]))
    print 'rms total = %.3f'%(np.std(Nrv_obs-Nrv_cal))
    print 'rms correlated (no LHS1140) = %.3f'%(np.std(Nrvtargets[:-1]-GPNrvs[:-1,0]))
    print 'rms total (no LHS1140) = %.3f'%(np.std(Nrv_obs[:-1]-Nrv_cal[:-1]))


    if pltt:
	plt.show()
    plt.close('all')

    
def plot_tobs_F(self, errorbar=False, pltt=True, label=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    g = self.Teffs_med > 0

    if errorbar:
        ax.errorbar(self.Fs_med[g], self.tobss_med_I[g], self.tobss_emed_N[g],
                    fmt='ro', ms=3, capsize=0, elinewidth=.4, alpha=.6)
        ax.errorbar(self.Fs_med[g], self.tobss_med_H[g], self.tobss_emed_H[g],
                    fmt='bo', ms=3, capsize=0, elinewidth=.4, alpha=.6)
    else:  # point estimates
        ax.scatter(self.Fs_med[g], self.tobss_med_N[g], edgecolor='r',
                   s=20,
                   alpha=1, facecolors='none', marker='d', label='nIR')
        ax.scatter(self.Fs_med[g], self.tobss_med_H[g], edgecolor='b',
                   s=20,
                   alpha=1, facecolors='none', marker='o', label='Optical')
        
    ax.set_xscale('log'), ax.set_yscale('log')
    ax.set_xlim((1e4,6e-2)), ax.set_ylim((1,1e4))
    ax.set_xlabel('Insolation ($S_{\oplus}$)')
    ax.set_ylabel('$t_{obs}$ [hours]')
    ax.legend(bbox_to_anchor=(.25,.96))

    fig.subplots_adjust(bottom=.14, top=.91)
    if label:
        plt.savefig('plots/tobs_F.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_Nrv_mag(self, mag='V', errorbar=False, pltt=True, label=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if mag == 'B':
        xarr = self.Bmags_med
    elif mag == 'V':
        xarr = self.Vmags_med
    elif mag == 'Y':
        xarr = self.Ymags_med
    elif mag == 'J':
        xarr = self.Jmags_med
    elif mag == 'H':
        xarr = self.Hmags_med
    elif mag == 'K':
        xarr = self.Kmags_med
        
    g = self.Teffs_med > 0

    if errorbar:
        ax.errorbar(xarr[g], self.Nrvs_med_N[g], self.Nrvs_emed_N[g],
                    fmt='ro', ms=3, capsize=0, elinewidth=.4, alpha=.5)
        ax.errorbar(xarr[g], self.Nrvs_med_H[g], self.Nrvs_emed_H[g],
                    fmt='bo', ms=3, capsize=0, elinewidth=.4, alpha=.5)
    else:  # point estimates
        ax.plot(xarr[g], self.Nrvs_med_N[g], 'ro', ms=8, alpha=.5)
        ax.plot(xarr[g], self.Nrvs_med_H[g], 'bo', ms=8, alpha=.5)
        
    #ax.set_xlim((1e4,1e-1))
    ax.set_yscale('log'), ax.set_ylim((1,1e3))
    ax.set_xlabel('$%s$'%mag), ax.set_ylabel('$n_{RV}$')
    
    if pltt:
        plt.show()
    plt.close('all')
    

def plot_Nrvratio(self, pltt=True, label=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ratio = unp.uarray(self.Nrvs_med_N, self.Nrvs_emed_N) / \
            unp.uarray(self.Nrvs_med_H, self.Nrvs_emed_H)
    
    ax.errorbar(self.Teffs_med, unp.nominal_values(ratio),
                unp.std_devs(ratio), fmt='k.', ms=12, 
		alpha=.5, capsize=0, elinewidth=.5)
    xlim = ax.get_xlim()
    ax.fill_between(list(xlim), np.ones(2), color='k', alpha=.2)

    # plot planetary systems
    names = ['TRAPPIST-1','Prox Cen','Ross 128','GJ 273']
    Teffs = [2559, 3050, 3192, 3150]
    for i in range(len(Teffs)):
	g = (self.Teffs_med >= Teffs[i]*.75) & \
            (self.Teffs_med <= Teffs[i]*1.25)
	ax.plot(Teffs[i], np.median(unp.nominal_values(ratio)[g]), 'go',
                ms=10)

    ax.set_yscale('log'), ax.set_ylim((1e-2,1e2)), ax.set_xlim(xlim)
    ax.set_xlabel('Teff (K)'), ax.set_ylabel('$n_{RV,I} / n_{RV,O}$')
    ax.minorticks_on()

    plt.subplots_adjust(bottom=.12)
    if label:
        plt.savefig('plots/Nrvcomparison.png')
    if pltt:
        plt.show()
    plt.close('all')

    

def plot_NIRPSvSPIROU(self, pltt=True, label=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    starnums, systnums = np.unique(self.starnums), np.unique(self.systnums)
    N_nirps, N_spirou = np.zeros(0), np.zeros(0)
    Teffs = np.zeros(0)
    for i in range(starnums.size):
        for j in range(systnums.size):
            g = (self.starnums == starnums[i]) & (self.systnums == systnums[j])
            if 'N' in self.spectrographs[g] and 'S' in self.spectrographs[g]:
                N_nirps = np.append(N_nirps,
                                    self.Nrvs[g & (self.spectrographs == 'N')])
                N_spirou = np.append(N_spirou,
                                     self.Nrvs[g & (self.spectrographs == 'S')])
                Teffs = np.append(Teffs, self.Teffs[g][0])

    ax.scatter(Teffs, N_nirps/N_spirou)
    ax.set_yscale('log'), ax.set_ylim((1e-3, 1e3))
    ax.set_xlabel('Teff [K]'), ax.set_ylabel('N_nirps / N_spirou')
    
    if pltt:
        plt.show()
    plt.close('all')

# point estimates
def plot_compare_tobs_H_N(self, pltt=True, label=False, GP=True):

    h = self.spectrographs == 'H'
    n = self.spectrographs == 'N'
    tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
    etobsH = self.tobsGPs_emed_H if GP else self.tobss_emed_H
    tobsN= self.tobsGPs_med_N if GP else self.tobss_med_N
    etobsN= self.tobsGPs_emed_N if GP else self.tobss_emed_N
    ratio = unp.uarray(tobsH / tobsN, etobsH / etobsN)

    fig = plt.figure(figsize=(4,4.3))
    ax = fig.add_subplot(111)

    #ax.errorbar(self.Teffs_med, unp.nominal_values(ratio),
    #            unp.std_devs(ratio),
    #            fmt='ko', ms=1, elinewidth=.1)
    # point estimates
    ax.fill_between([25e2,12e3], 2, 1./2, color='k', alpha=.1)
    inds = np.arange(self.nstars)
    inds = np.delete(inds, np.array([486, 1485, 1645]))
    np.random.shuffle(inds)
    img = ax.scatter(self.Teffs_med[inds], unp.nominal_values(ratio)[inds],
                     s=12, alpha=1, norm=LogNorm(), vmin=.01, vmax=10,
                     cmap=_truncate_colormap(plt.get_cmap('rainbow'),0,1),
                     c=(self.sigmaRV_phot_med_H/self.sigmaRV_phot_med_N)[inds])
    cbar_axes = fig.add_axes([.07, .12, .86, .05])
    cbar = fig.colorbar(img, cax=cbar_axes, orientation='horizontal')
    cbar.set_label('$\sigma_{RV,opt}$ / $\sigma_{RV,nIR}$', labelpad=0)
    ax.text(8.5e3, 3, 'Near-IR (YJH)\nfavored', fontsize=10)
    ax.text(8.6e3, .33, 'Optical (BVR)\nfavored', verticalalignment='top', fontsize=10)

    ax.set_yscale('log')
    ax.set_ylim((5e-3, 100)), ax.set_xlim((25e2, 12e3))
    ax.set_xticks(np.arange(3e3,13e3,3e3))
    ax.set_xticklabels(np.arange(3e3,13e3,3e3,dtype=int), fontsize=12)
    ax.set_ylabel('t$_{obs,opt}$ / t$_{obs,nIR}$', labelpad=0)
    ax.set_xlabel('Effective Temperature [K]')
    ax.minorticks_on()

    fig.subplots_adjust(left=.19, bottom=.29, right=.93, top=.97)
    if label:
        plt.savefig('plots/optnIRratio.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_cumulative_detections_v_tobs(self, pltt=True, label=False,
                                      pltflag=0, tmax=1e3,
                                      xscale='linear', GP=True):
    mpl.rc('xtick', labelsize=12)
    mpl.rc('ytick', labelsize=12)

    fig = plt.figure(figsize=(9.5,4.8))
    ymaxs = 40,110,240,150
    tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
    tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N

    for i in range(rpedges.size-1):
	ax1 = fig.add_subplot(2, rpedges.size-1, i+1)
        ax2 = fig.add_subplot(2, rpedges.size-1, i+1+rpedges.size-1)

	# HARPS total observing time
	if pltflag in [0,2]:
            g = (self.rps_med >= rpedges[i]) & \
                (self.rps_med < rpedges[i+1])
            tobs = np.sort(self.tobss_med_H[g])
	    tobs = np.append(0, np.cumsum(tobs))
	    Ndet = np.arange(tobs.size)
	    g = tobs <= tmax+10
	    ax1.plot(tobs[g], Ndet[g], '--', c=colH, lw=2, drawstyle='steps')
            # plot derivative
	    if i == 0:
                ax1.plot([50,200], np.repeat(35,2), '--', lw=2, c=colH)
		ax1.text(240, 34.5, 'Optical (BVR)', color=colH, fontsize=10, 
			 weight='semibold', verticalalignment='center')
            # plot derivative
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
            
        if pltflag in [0,1]:
            threshold = 1./20
            ax2.plot([0,tmax*1e3], np.repeat(threshold,2), 'k--', lw=.9)
            if i == 0:
                ax2.text(190, threshold*1.05,
                         '%i hours / detection'%(1./threshold),
                         verticalalignment='bottom', fontsize=9)
                
        # NIRPS total observing time
	if pltflag in [0,3]:
            g = (self.rps_med >= rpedges[i]) & (self.rps_med < rpedges[i+1])
            tobs = np.sort(self.tobss_med_N[g])
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
	    g = tobs <= tmax+10
            #_compute_curve_derivative(tobs[g], Ndet[g])
	    ax1.plot(tobs[g], Ndet[g], '-', c=colN, drawstyle='steps')
	    if i == 0:
                ax1.plot([50,200], np.repeat(29,2), '-', c=colN)
                ax1.text(240, 29, 'Near-IR (YJH)', color=colN, fontsize=10, 
                         weight='semibold', verticalalignment='center')
            # plot derivative
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, '-', c=colN)

        if pltflag in [0,1]:
            ax1.set_title('%s\n%s'%(rplabels[i], rplabels2[i]), fontsize=12,
                          y=1.33)
        ax2.set_xlabel('Cumulative observing\ntime [hours]', fontsize=9)
        ax1.set_xscale(xscale), ax1.set_xticklabels('')
        ax1.minorticks_on(), ax2.minorticks_on()
        ax1.set_ylim((0,ymaxs[i]))
        ax2.set_xscale(xscale)
    	ax2.set_yscale('log'), ax2.set_ylim((1e-2,1.3))
    	ax1.set_xlim((0,tmax)), ax2.set_xlim((0,tmax))
        ax2.minorticks_on()
        
        if i == 0:
	    ax1.set_ylabel('Total number of\nplanet detections',
                           fontsize=10)
            ax2.set_ylabel('dN/dt\n[detections / hour]', fontsize=10,
                           labelpad=0)
            ax2.set_yticks(np.logspace(-2,0,3))
            ax2.set_yticklabels(['0.01','0.1','1'])
        else:
            ax2.set_yticklabels('')

        # second derivative axis
        ax3 = ax2.twinx()
        ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
        ax3.set_yscale('log')
        if i < rpedges.size-2:
            ax3.set_yticklabels('')
        else:
            ax3.set_yticks(np.logspace(2,0,3))
            ax3.set_yticklabels(['%i'%i for i in np.logspace(2,0,3)])
            ax3.set_ylabel('dt/dN\n[hours / detection]', fontsize=10,
                           labelpad=0)

        # second x-axis
        ax4 = ax1.twiny()
        ax4.set_xlim(tuple(np.ascontiguousarray(ax1.get_xlim())/7.))
        ax4.set_xticks(np.arange(0,150,50))
        ax4.set_xticklabels(np.arange(0,150,50), va='baseline')
        ax4.set_xlabel('Cumulative observing\ntime [nights]', fontsize=9,
                       labelpad=3)
        ax4.minorticks_on()

        if pltflag not in [0,1]:
            ax1.axis('off')
            ax2.axis('off')
            ax3.axis('off')
            ax4.axis('off')
        
    fig.subplots_adjust(left=.085, bottom=.125, top=.8, right=.925,
                        hspace=0, wspace=.25)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        transparent = False if pltflag in [0,1] else True 
        GPlabel = 'GP' if GP else ''
        plt.savefig('plots/cumulativetobs%s_rp_%s.png'%(GPlabel,label),
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')


def plot_cumulative_detections_v_Nrv(self, pltt=True, label=False,
                                     harps=True, nirps=True, Nmax=1e3,
                                     observebesttargetsfirst=True):
    fig = plt.figure(figsize=(9.5,4.6))

    for i in range(rpedges.size-1):
	ax1 = fig.add_subplot(2, rpedges.size-1, i+1)

	# HARPS total observing time
	if harps:
            g = (self.rps_med >= rpedges[i]) & \
                (self.rps_med < rpedges[i+1])
            Nrvs = np.sort(self.Nrvs_med_H[g]) if observebesttargetsfirst \
                   else self.Nrvs_med_H[g]
	    Nrvs = np.append(0, np.cumsum(Nrvs))
	    Ndet = np.arange(Nrvs.size)
	    g = Nrvs <= Nmax
	    ax1.plot(Nrvs[g], Ndet[g], '--', c=colH, lw=2, drawstyle='steps')
            # plot derivative
	    if i == 0:
                ax1.plot([410,560], np.repeat(3.2,2), '--', lw=2, c=colH)
		ax1.text(6e2, 3.2, 'Optical', color=colH, fontsize=10, 
			 weight='semibold', verticalalignment='center')
            # plot derivative
            ax2 = fig.add_subplot(2, rpedges.size-1, i+1+rpedges.size-1)
            Nrvs2, dNdt = _compute_curve_derivative(Nrvs[g], Ndet[g])
            ax2.plot(Nrvs2, dNdt, '--', c=colH, lw=2)
            threshold = 1./20
            ax2.plot([0,Nmax*1e3], np.repeat(threshold,2), 'k--', lw=.9)
            if i == 0:
                ax2.text(Nmax*.2, threshold*1.03,
                         '%i measurements /\n     detection'%(1./threshold),
                         verticalalignment='bottom', fontsize=9)
                
        # NIRPS total observing time
	if nirps:
            g = (self.rps_med >= rpedges[i]) & (self.rps_med < rpedges[i+1])
            Nrvs = np.sort(self.Nrvs_med_N[g]) if observebesttargetsfirst \
                   else self.Nrvs_med_N[g]
            Nrvs = np.append(0, np.cumsum(Nrvs))
            Ndet = np.arange(Nrvs.size)
	    g = Nrvs <= Nmax
	    ax1.plot(Nrvs[g], Ndet[g], '-', c=colN, drawstyle='steps')
	    if i == 0:
                ax1.plot([410,550], np.repeat(1.2,2), '-', c=colN)
                ax1.text(6e2, 1.2, 'Near-IR', color=colN, fontsize=10, 
                         weight='semibold', verticalalignment='center')
            # plot derivative
            Nrvs2, dNdt = _compute_curve_derivative(Nrvs[g], Ndet[g])
            ax2.plot(Nrvs2, dNdt, '-', c=colN)

        ax1.set_title('%s\n%s'%(rplabels[i], rplabels2[i]), fontsize=10)
        ax2.set_xlabel('Cumulative n$_{RV}$', fontsize=11)
        ax1.set_xticklabels('')
        ax1.minorticks_on(), ax2.minorticks_on()
    	ax2.set_yscale('log'), ax2.set_ylim((1e-2,1.3))
    	ax1.set_xlim((0,Nmax)), ax2.set_xlim((0,Nmax))

        if i == 0:
	    ax1.set_ylabel('Total number of\nplanet detections',
                           fontsize=10)
            ax2.set_ylabel('dN/dn$_{RV}$\n[detections / measurement]',
                           fontsize=10, labelpad=0)
            ax2.set_yticks(np.logspace(-2,0,3))
            ax2.set_yticklabels(['0.01','0.1','1'])
        else:
            ax2.set_yticklabels('')

        # second derivative axis
        ax3 = ax2.twinx()
        ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
        ax3.set_yscale('log')
        if i < rpedges.size-2:
            ax3.set_yticklabels('')
        else:
            ax3.set_yticks(np.logspace(2,0,3))
            ax3.set_yticklabels(['%i'%i for i in np.logspace(2,0,3)])
            ax3.set_ylabel('dn$_{RV}$/dN\n[measurements / detection]',
                           fontsize=10, labelpad=0)
        
    fig.subplots_adjust(left=.085, bottom=.125, top=.92, right=.92,
                        hspace=0, wspace=.22)
    if label:
        plt.savefig('plots/cumulativeNrv_rp.png')    
    if pltt:
	plt.show()
    plt.close('all')

    

def plot_cumulative_detections_v_tobs_Teff(self, pltt=True, label=False,
                                           tmax=1e3,
                                           observebesttargetsfirst=True,
                                           GP=True, pltflag=0):
    '''pltflag: 0 (full), 1 (bkgd), 2 (optical), 3 (nearIR)'''
    mpl.rc('xtick', labelsize=12)
    mpl.rc('ytick', labelsize=12)

    fig = plt.figure(figsize=(9.5,4.8))
    ymaxs = 60,200,300,40
    tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
    tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
    
    for i in range(Teffedges.size-1):
	ax1 = fig.add_subplot(2, Teffedges.size-1, -1*i+4)
        ax2 = fig.add_subplot(2, Teffedges.size-1, -1*i+4+Teffedges.size-1)
            
	# HARPS total observing time
	if pltflag in [0,2]:
            g = (self.Teffs_med >= Teffedges[i]) & \
                (self.Teffs_med < Teffedges[i+1])
            tobs = np.sort(tobsH[g]) if observebesttargetsfirst \
                   else tobsH[g]
	    tobs = np.append(0, np.cumsum(tobs))
	    Ndet = np.arange(tobs.size)
	    g = tobs <= tmax+100
	    ax1.plot(tobs[g], Ndet[g], '--', c=colH, lw=2, drawstyle='steps')
            # plot derivative
	    if i == 3:
                ax1.plot([410,560], np.repeat(29,2), '--', c=colH, lw=2)
		ax1.text(6e2, 28, 'Optical\n(BVR)', color=colH, fontsize=9,
			 weight='semibold', verticalalignment='center')
            # plot derivative
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
        threshold = 1./20
        if pltflag in [0,1]:
            ax2.plot([0,tmax], np.repeat(threshold,2), 'k--', lw=.9)
            if i == 2:
                ax2.text(100, threshold*1.03,
                         '%i hours / detection'%(1./threshold),
                         verticalalignment='bottom', fontsize=10)
                    
        # NIRPS total observing time
	if pltflag in [0,3]:
            g = (self.Teffs_med >= Teffedges[i]) & \
                (self.Teffs_med < Teffedges[i+1])
            tobs = np.sort(tobsN[g]) if observebesttargetsfirst \
                   else tobsN[g]
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
	    g = tobs <= tmax+100
            #_compute_curve_derivative(tobs[g], Ndet[g])
	    ax1.plot(tobs[g], Ndet[g], '-', c=colN, drawstyle='steps')
	    if i == 3:
                ax1.plot([410,550], np.repeat(20,2), 'r-')
                ax1.text(6e2, 19, 'Near-IR\n(YJH)', color=colN, fontsize=9, 
                         weight='semibold', verticalalignment='center')
            # plot derivative
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, '-', c=colN)
            
        if pltflag in [0,1]:
            ax1.set_title('%s\n%i $> $T$_{eff} \geq$ %i'%(Tefflabels[i],
                                     		          Teffedges[i+1],
                                     		          Teffedges[i]),
                          fontsize=12, y=1.33)
        ax2.set_xlabel('Cumulative observing\ntime [hours]', fontsize=9)
        ax1.set_xticklabels('')
        ax1.set_ylim((0,ymaxs[i]))
        ax1.minorticks_on(), ax2.minorticks_on()
    	ax2.set_yscale('log'), ax2.set_ylim((1e-2,1.3))
    	ax1.set_xlim((0,tmax)), ax2.set_xlim((0,tmax))
        ax2.set_xticks(np.arange(0,1001,500))
        ax2.set_xticklabels(np.arange(0,1001,500), fontsize=12)
        
        if i == 3:
	    ax1.set_ylabel('Total number of\nplanet detections',
                           fontsize=10)
            ax2.set_ylabel('dN/dt\n[detections / hour]', fontsize=10,
                           labelpad=0)
            ax2.set_yticks(np.logspace(-2,0,3))
            ax2.set_yticklabels(['0.01','0.1','1'])
        else:
            ax2.set_yticklabels('')

        # second derivative axis
        ax3 = ax2.twinx()
        ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
        ax3.set_yscale('log')
        if i == 0:
            ax3.set_yticks(np.logspace(2,0,3))
            ax3.set_yticklabels(['%i'%i for i in np.logspace(2,0,3)])
            ax3.set_ylabel('dt/dN\n[hours / detection]', fontsize=10,
                           labelpad=0)
        else:
            ax3.set_yticklabels('')

        # second x-axis
        ax4 = ax1.twiny()
        ax4.set_xlim(tuple(np.ascontiguousarray(ax1.get_xlim())/7.))
        ax4.set_xticks(np.arange(0,150,50))
        ax4.set_xticklabels(np.arange(0,150,50), va='baseline')
        ax4.set_xlabel('Cumulative observing\ntime [nights]', fontsize=9,
                       labelpad=3)
        ax4.minorticks_on()

        if pltflag not in [0,1]:
            ax1.axis('off')
            ax2.axis('off')
            ax3.axis('off')
            ax4.axis('off')
        
    fig.subplots_adjust(left=.085, bottom=.125, top=.8, right=.925,
                        hspace=0, wspace=.25)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        transparent = False if pltflag in [0,1] else True 
        GPlabel = 'GP' if GP else ''
        plt.savefig('plots/cumulativetobs%s_Teff_%s.png'%(GPlabel,label),
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')


def plot_cumulative_detections_v_Nrv_Teff(self, pltt=True, label=False,
                                          harps=True, nirps=True, Nmax=1e3,
                                          observebesttargetsfirst=True):
    fig = plt.figure(figsize=(9.5,4.6))

    for i in range(Teffedges.size-1):
	ax1 = fig.add_subplot(2, Teffedges.size-1, -1*i+4)

	# HARPS total observing time
	if harps:
            g = (self.Teffs_med >= Teffedges[i]) & \
                (self.Teffs_med < Teffedges[i+1])
            Nrvs = np.sort(self.Nrvs_med_H[g]) if observebesttargetsfirst \
                   else self.Nrvs_med_H[g]
	    Nrvs = np.append(0, np.cumsum(Nrvs))
	    Ndet = np.arange(Nrvs.size)
	    g = Nrvs <= Nmax
	    ax1.plot(Nrvs[g], Ndet[g], '--', c=colH, lw=2, drawstyle='steps')
            # plot derivative
	    if i == 3:
                ax1.plot([410,560], np.repeat(18,2), '--', c=colH, lw=2)
		ax1.text(6e2, 18, 'Optical', color=colH, fontsize=10,
			 weight='semibold', verticalalignment='center')
            # plot derivative
            ax2 = fig.add_subplot(2, Teffedges.size-1, -1*i+4+Teffedges.size-1)
            Nrvs2, dNdt = _compute_curve_derivative(Nrvs[g], Ndet[g])
            ax2.plot(Nrvs2, dNdt, '--', c=colH, lw=2)
            threshold = 1./20
            ax2.plot([0,Nmax], np.repeat(threshold,2), 'k--', lw=.9)
            if i == 3:
                ax2.text(Nmax*.2, threshold*1.03,
                         '%i measurements /\n     detection'%(1./threshold),
                         verticalalignment='bottom', fontsize=9)
                
        # NIRPS total observing time
	if nirps:
            g = (self.Teffs_med >= Teffedges[i]) & \
                (self.Teffs_med < Teffedges[i+1])
            Nrvs = np.sort(self.Nrvs_med_N[g]) if observebesttargetsfirst \
                   else self.Nrvs_med_N[g]
            Nrvs = np.append(0, np.cumsum(Nrvs))
            Ndet = np.arange(Nrvs.size)
	    g = Nrvs <= Nmax+100
	    ax1.plot(Nrvs[g], Ndet[g], '-', c=colN, drawstyle='steps')
	    if i == 3:
                ax1.plot([410,550], np.repeat(14,2), 'r-')
                ax1.text(6e2, 14, 'Near-IR', color=colN, fontsize=10, 
                         weight='semibold', verticalalignment='center')
            # plot derivative
            Nrvs2, dNdt = _compute_curve_derivative(Nrvs[g], Ndet[g])
            ax2.plot(Nrvs2, dNdt, '-', c=colN)
            
        ax1.set_title('%s\n%i $> T_{eff} \geq$ %i'%(Tefflabels[i],
                                     		    Teffedges[i+1],
                                     		    Teffedges[i]),
                      fontsize=10)
        ax2.set_xlabel('Cumulative n$_{RV}$', fontsize=11)
        ax1.set_xticklabels('')
        ax1.minorticks_on(), ax2.minorticks_on()
    	ax2.set_yscale('log'), ax2.set_ylim((1e-2,1.3))
    	ax1.set_xlim((0,Nmax)), ax2.set_xlim((0,Nmax))

        if i == 3:
	    ax1.set_ylabel('Total number of\nplanet detections',
                           fontsize=10)
            ax2.set_ylabel('dN/dn$_{RV}$\n[detections / measurement]',
                           fontsize=10, labelpad=0)
            ax2.set_yticks(np.logspace(-2,0,3))
            ax2.set_yticklabels(['0.01','0.1','1'])
        else:
            ax2.set_yticklabels('')

        # second derivative axis
        ax3 = ax2.twinx()
        ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
        ax3.set_yscale('log')
        if i == 0:
            ax3.set_yticks(np.logspace(2,0,3))
            ax3.set_yticklabels(['%i'%i for i in np.logspace(2,0,3)])
            ax3.set_ylabel('dn$_{RV}$/dN\n[measuremens / detection]',
                           fontsize=10, labelpad=0)
        else:
            ax3.set_yticklabels('')
            
    fig.subplots_adjust(left=.085, bottom=.125, top=.92, right=.92,
                        hspace=0, wspace=.23)
    if label:
        plt.savefig('plots/cumulativeNrv_Teff.png')
    if pltt:
	plt.show()
    plt.close('all')


def plot_cumulative_detections_v_tobs_50_random(self, pltt=True, label=False,
                                                Nrand=10, seed=0, pltflag=0,
                                                sigma=5., tmax=1e5, GP=True,
                                                xscale='log'):
    '''0=all, 1=bkgd, 2=optical curves, 3=nIR curves''' 
    fig = plt.figure(figsize=(5.1,6.3))# 5.2
    gs = gridspec.GridSpec(10,1)
    ax1 = plt.subplot(gs[:7,0])
    ax2 = plt.subplot(gs[7:,0])
    ax4 = fig.add_axes([.22, .67, .2, .18])
    #ax1 = fig.add_subplot(211)
    #ax2 = fig.add_subplot(212)
    np.random.seed(int(seed))
    corr = (.327/.189)**2 if sigma == 5 else 1
    g = self.rps_med < 4
    tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
    tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
    
    # HARPS total observing time
    if pltflag in [0,2]:
        tobs = np.sort(tobsH[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))[:51]
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax2.text(.51, .79, 'Optical (BVR)', color=colH, fontsize=13,
                 weight='semibold', transform=ax2.transAxes)
        ax2.plot([1e2,250], np.repeat(1,2), '--', c=colH)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
    if pltflag in [0,1]:
        threshold = 1./20  # detections / hour
        ax2.plot([1,tmax], np.repeat(threshold,2), 'k--')
        ax2.text(4, threshold*1.06,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='bottom', fontsize=10)
    if pltflag in [0,2]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = tobsH[g][inds]  * corr
            tobs = np.append(0, np.cumsum(tobs))[:51]
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)

    # NIRPS total observing time
    if pltflag in [0,3]:
        tobs = np.sort(tobsN[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))[:51]
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
        ax2.text(.51, .63, 'Near-IR (YJH)', color=colN, fontsize=13, 
		 weight='semibold', transform=ax2.transAxes)
        ax2.plot([1e2,250], np.repeat(.5,2), '-', c=colN)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '-', c=colN)
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = tobsN[g][inds]  * corr
            tobs = np.append(0, np.cumsum(tobs))[:51]
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    # plot rp histogram
    if pltflag in [0,1]:
        nbins = 8
        best = np.unique(np.append(get_planets(tobsH, g, N=50),
                                   get_planets(tobsN, g, N=50))) 
        ax4.hist(self.rps_med[best], bins=np.linspace(1,4,nbins), color='k',
                 alpha=.5)
        ax4.hist(self.rps_med[best], bins=np.linspace(1,4,nbins), color='k',
                 lw=1.5, histtype='step')
        ax4.set_xticks(np.arange(1,5))
        ax4.set_xticklabels(np.arange(1,5), fontsize=9), ax4.set_yticklabels('')
        ax4.minorticks_on()
        ax4.set_xlabel('r$_p$ [R$_{\oplus}$]', fontsize=10, labelpad=0)

    ax1.set_xscale(xscale), ax1.set_xlim((1,tmax))
    ax1.set_ylim((0,50))
    ax1.set_xticklabels('')
    ax2.set_yscale('log'), ax2.set_ylim((3e-2,2))
    ax2.set_xscale(xscale), ax2.set_xlim((1,tmax))
    ax1.minorticks_on()
    if pltflag in [0,1]:
        ax1.set_ylabel('Total number of\nplanet detections', fontsize=13)
        ax1.set_yticks(np.arange(0,51,10))
        ax1.set_yticklabels(np.arange(0,51,10), fontsize=13)
        ax2.set_yticks([.1,1])
        ax2.set_yticklabels(['0.1','1'], fontsize=13)
        ax2.set_xlabel('Cumulative observing time [hours]', fontsize=13)
        ax2.set_ylabel('dN / dt\n[detections / hour]', fontsize=12, labelpad=0)
        ax1.set_title('TESS level one science requirement', fontsize=13,
                      y=1.14)
    else:
        ax1.set_yticklabels('')
        ax2.set_xticklabels('')
        ax2.set_yticklabels('')
        
    # second derivative axis
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    if pltflag in [0,1]:
        ax3.set_ylabel('dt / dN\n[hours / detection]', fontsize=12, labelpad=4)
        ax3.set_yticks([10,1])
        ax3.set_yticklabels(['10','1'], fontsize=13)
    else:
        ax3.set_yticklabels('')

    # second x-axis
    ax5 = ax1.twiny()
    ax5.set_xlim(tuple(np.ascontiguousarray(ax1.get_xlim())/7.))
    if pltflag in [0,1]:
        ax5.set_xticks(10**np.arange(0,5))
        ax5.set_xticklabels(['10$^{%i}$'%i for i in np.arange(0,5)],
                            fontsize=12, va='baseline')
        ax5.set_xlabel('Cumulative observing time [nights]', fontsize=11)
        ax5.set_xlim((1/7.,tmax/7.)), ax5.set_xscale('log')
        ax5.minorticks_on()
    else:
        ax5.set_yticklabels('')

    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
        ax4.axis('off')
        ax5.axis('off')
    fig.subplots_adjust(left=.19, bottom=.09, top=.89, right=.84, hspace=0)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        elif pltflag == 4:
            label = 'thres'
        transparent = False if pltflag in [0,1] else True 
        GPlabel = 'GP' if GP else ''
        plt.savefig('plots/cumulativetobs%s_50random4_%s.png'%(GPlabel,label),
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')


def plot_cumulative_detections_v_Nrv_50_random(self, pltt=True, label=False,
                                               Nrand=10, seed=0, pltflag=0,
                                               sigma=3.):
    '''0=all, 1=bkgd, 2=optical curves, 3=nIR curves, 4=threshold line''' 
    fig = plt.figure(figsize=(5.5,5.4))# 5.2
    gs = gridspec.GridSpec(10,1)
    ax1 = plt.subplot(gs[:7,0])
    ax2 = plt.subplot(gs[7:,0])
    #ax1 = fig.add_subplot(211)
    #ax2 = fig.add_subplot(212)
    np.random.seed(int(seed))
    corr = (.327/.189)**2 if sigma == 5 else 1
    g = self.rps_med < 4
    
    # HARPS total observing time
    if pltflag in [0,2]:
        Nrvs = np.sort(self.Nrvs_med_H[g] * corr)
        Nrvs = np.append(0, np.cumsum(Nrvs))[:51]
        Ndet = np.arange(Nrvs.size)
        ax1.plot(Nrvs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax2.text(.7, .79, 'Optical', color=colH, fontsize=13,
                 weight='semibold', transform=ax2.transAxes)
        # plot derivative
        Nrvs2, dNdt = _compute_curve_derivative(Nrvs, Ndet)
        ax2.plot(Nrvs2, dNdt, '--', c=colH, lw=2)
    if pltflag in [0,1]:
        threshold = 1./20  # detections / hour
        ax2.plot([1,1e5], np.repeat(threshold,2), 'k--')
        ax2.text(4, threshold*1.06,
                 '%i measurements / detection'%(1./threshold),
                 verticalalignment='bottom', fontsize=9)
    if pltflag in [0,2]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            Nrvs = self.Nrvs_med_H[g][inds]  * corr
            Nrvs = np.append(0, np.cumsum(Nrvs))[:51]
            Ndet = np.arange(Nrvs.size)
            ax1.plot(Nrvs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)

    # NIRPS total observing time
    if pltflag in [0,3]:
        Nrvs = np.sort(self.Nrvs_med_N[g] * corr)
        Nrvs = np.append(0, np.cumsum(Nrvs))[:51]
	Ndet = np.arange(Nrvs.size)
	ax1.plot(Nrvs, Ndet, '-', c=colN, drawstyle='steps')
	ax2.text(.7, .63, 'Near-IR', color=colN, fontsize=13, 
		 weight='semibold', transform=ax2.transAxes)
        Nrvs2, dNdt = _compute_curve_derivative(Nrvs, Ndet)
        ax2.plot(Nrvs2, dNdt, '-', c=colN)
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            Nrvs = self.Nrvs_med_N[g][inds]  * corr
            Nrvs = np.append(0, np.cumsum(Nrvs))[:51]
            Ndet = np.arange(Nrvs.size)
            ax1.plot(Nrvs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    ax1.set_xscale('log'), ax1.set_xlim((1,1e5))
    ax1.set_ylim((0,50))
    ax1.set_xticklabels('')
    ax2.set_yscale('log'), ax2.set_ylim((3e-2,2))
    ax2.set_xscale('log'), ax2.set_xlim((1,1e5))
    ax1.minorticks_on()
    if pltflag in [0,1]:
        ax1.set_ylabel('Total number of\nplanet detections', fontsize=13)
        ax1.set_yticks(np.arange(0,51,10))
        ax1.set_yticklabels(np.arange(0,51,10), fontsize=13)
        ax2.set_yticks([.1,1])
        ax2.set_yticklabels(['0.1','1'], fontsize=13)
        ax2.set_xlabel('Cumulative n$_{RV}$', fontsize=13)
        ax2.set_ylabel('dN / dn$_{RV}$\n[detections / measurement]',
                       fontsize=12, labelpad=0)
        ax1.set_title('r$_p <$ 4 R$_{\oplus}$', fontsize=13)
    else:
        ax1.set_yticklabels('')
        ax2.set_xticklabels('')
        ax2.set_yticklabels('')
        
    # second derivative axis
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    if pltflag in [0,1]:
        ax3.set_ylabel('dn$_{RV}$ / dN\n[measurements / detection]',
                       fontsize=12, labelpad=4)
        ax3.set_yticks([10,1])
        ax3.set_yticklabels(['10','1'], fontsize=13)
    else:
        ax3.set_yticklabels('')
    
    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
    fig.subplots_adjust(left=.15, bottom=.1, top=.96, right=.86, hspace=0)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        elif pltflag == 4:
            label = 'thres'
        transparent = False if pltflag in [0,1] else True 
        plt.savefig('plots/cumulativeNrv_50random_%s.png'%label,
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')


    
def plot_cumulative_detections_v_tobs_MR(self, pltt=True, label=False,
                                         tmax=1e5, Nrand=10, seed=0, pltflag=0,
                                         GP=True, xscale='log', sigma=5.):
    '''0=all, 1=bkgd, 2=optical curves, 3=nIR curves, 4=threshold line'''
    fig = plt.figure(figsize=(5.5,5.4))# 5.2
    gs = gridspec.GridSpec(10,1)
    ax1 = plt.subplot(gs[:7,0])
    ax2 = plt.subplot(gs[7:,0])
    ax4 = fig.add_axes([.24, .6, .32, .32])
    np.random.seed(int(seed))
    corr = (.327/.189)**2
    g = (self.rps_med >= 1.5) & (self.rps_med <= 2)

    # NIRPS total observing time
    if pltflag in [0,3]:
        tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
        tobs = np.sort(tobsN[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
	ax2.text(.7, .57, 'Near-IR', color=colN, fontsize=13, 
		 weight='semibold', transform=ax2.transAxes)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '-', c=colN)
        ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
            tobs = tobsN[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    # HARPS total observing time
    if pltflag in [0,2]:
        tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
        tobs = np.sort(tobsH[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax2.text(.7, .73, 'Optical', color=colH, fontsize=13,
                 weight='semibold', transform=ax2.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
        # plot linear subplot
        ax4.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
    if pltflag in [0,1]:
        threshold = 1./20  # detections / hour
        ax2.plot([1,tmax], np.repeat(threshold,2), 'k--')
        ax2.text(1e3, threshold*1.06,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='bottom', fontsize=9)
    if pltflag in [0,2]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
            tobs = tobsH[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            ax4.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)
            
    ax1.set_xscale(xscale)
    ax1.set_xlim((1,tmax))
    if xscale=='log':
        ax1.set_ylim((0,400))
    else:
        ax1.set_ylim((0,60))
    ax1.set_xticklabels('')
    ax2.set_yscale('log'), ax2.set_ylim((1e-2,2))
    ax2.set_xscale(xscale), ax2.set_xlim((1,tmax))
    ax1.minorticks_on()
    ax4.set_xlim((0,1e3)), ax4.set_ylim((0,70))
    ax4.set_xticks(np.linspace(0,1e3,5))
    ax4.set_xticklabels(['0','','500','','1000'], fontsize=12)
    ax4.set_yticks(np.linspace(0,50,6))
    ax4.set_yticklabels(np.linspace(0,50,6,dtype=int), fontsize=12)
    ax4.minorticks_on()
    if pltflag in [0,1]:
        ax1.set_ylabel('Total number of\nplanet detections', fontsize=13)
        #ax1.set_yticks(np.arange(0,51,10))
        #ax1.set_yticklabels(np.arange(0,51,10), fontsize=13)
        ax2.set_xlabel('Cumulative observing time [hours]', fontsize=13)
        ax2.set_ylabel('dN / dt\n[detections / hour]', fontsize=12, labelpad=0)
        ax1.set_title('1.5$\leq$ r$_p$/R$_{\oplus} \leq$2', fontsize=13)
    else:
        ax1.set_yticklabels('')
        ax2.set_xticklabels('')
        ax2.set_yticklabels('')
        
    # second derivative axis
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    if pltflag in [0,1]:
        ax3.set_ylabel('dt / dN\n[hours / detection]', fontsize=12, labelpad=0)
        ax3.set_yticks([100,10,1])
        ax3.set_yticklabels(['100','10','1'], fontsize=13)
    else:
        ax3.set_yticklabels('')
    
    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
        ax4.axis('off')
    fig.subplots_adjust(left=.17, bottom=.1, top=.96, right=.84, hspace=0)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        transparent = False if pltflag in [0,1] else True
        GPlabel = 'GP' if GP else ''
        plt.savefig('plots/cumulativetobs%s_MR_%s.png'%(GPlabel,label),
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')

    
def plot_cumulative_detections_v_tobs_Fulton(self, pltt=True, label=False,
                                             tmax=4e5, Nrand=10, seed=0,
                                             pltflag=0, GP=True, xscale='log',
                                             sigma=5.):
    '''0=all, 1=bkgd, 2=optical curves, 3=nIR curves, 4=threshold line'''
    fig = plt.figure(figsize=(5.1,6.3))# 4.8
    gs = gridspec.GridSpec(10,1)
    ax1 = plt.subplot(gs[:7,0])
    ax2 = plt.subplot(gs[7:,0])
    ax4 = fig.add_axes([.27, .57, .3, .25])
    np.random.seed(int(seed))
    corr = (.327/.189)**2 if sigma == 5 else 1.
    rp_upper = 10**(-0.09*np.log10(self.Ps_med)+.44)
    rp_lower = 10**(-0.09*np.log10(self.Ps_med)+.29)
    g1 = (self.rps_med <= 2) & \
         (self.rps_med >= rp_lower) & (self.rps_med <= rp_upper)
    g2 = (self.rps_med > 2) & \
         (self.rps_med >= rp_lower) & (self.rps_med <= rp_upper)
    g = g1 + g2
    print '%i planets'%g.sum()

    # mix each class of planet so we capture both in the best cumulative tobs
    inds1 = get_planets(self.tobsGPs_med_N, g1)
    inds2 = get_planets(self.tobsGPs_med_N, g2)
    sortN = np.zeros(inds1.size + inds2.size)
    lim = inds1.size if inds1.size < inds2.size else inds2.size
    sortN[np.arange(0,lim*2,2)] = inds1  # may have to switch to inds2 here
    sortN[np.arange(1,lim*2,2)] = inds2[:lim]
    leftover = abs(inds2.size-inds1.size)
    assert not np.any(sortN[-leftover:] > 0)
    sortN[-leftover:] = inds2[-leftover:]
    sortN = sortN.astype(int)

    inds1 = get_planets(self.tobsGPs_med_H, g1)
    inds2 = get_planets(self.tobsGPs_med_H, g2)
    sortH = np.zeros(inds1.size + inds2.size)
    lim = inds1.size if inds1.size < inds2.size else inds2.size
    sortH[np.arange(0,lim*2,2)] = inds1  # may have to switch to inds2 here
    sortH[np.arange(1,lim*2,2)] = inds2[:lim]
    leftover = abs(inds2.size-inds1.size)
    assert not np.any(sortH[-leftover:] > 0)
    sortH[-leftover:] = inds2[-leftover:]
    sortH = sortH.astype(int)
    
    # NIRPS total observing time
    if pltflag in [0,3]:
        tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
        #tobs1 = np.sort(tobsN[g] * corr)
        tobs = tobsN[sortN]*corr
        tobs = np.append(0, np.cumsum(tobs))
        Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
	ax2.text(.55, .66, 'Near-IR (YJH)', color=colN, fontsize=13, 
		 weight='semibold', transform=ax2.transAxes)
        ax2.plot([3e2,9e2], np.repeat(.4,2), '-', c=colN)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '-', c=colN)
        ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
            tobs = tobsN[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    # HARPS total observing time
    if pltflag in [0,2]:
        tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
        #tobs = np.sort(tobsH[g] * corr)
        tobs = tobsH[sortH]*corr
        tobs = np.append(0, np.cumsum(tobs))
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax2.text(.55, .82, 'Optical (BVR)', color=colH, fontsize=13,
                 weight='semibold', transform=ax2.transAxes)
        ax2.plot([3e2,9e2], np.repeat(1,2), '--', c=colH)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
        # plot linear subplot
        ax4.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
    if pltflag in [0,1]:
        threshold = 1./20  # detections / hour
        ax2.plot([1,tmax], np.repeat(threshold,2), 'k--')
        ax2.text(1.5, threshold*.86,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='top', fontsize=10)
    if pltflag in [0,2]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
            tobs = tobsH[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            ax4.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)
            
    ax1.set_xscale(xscale)
    ax1.set_xlim((1,tmax))
    if xscale=='log':
        ax1.set_ylim((0,g.sum()))
    else:
        ax1.set_ylim((0,60))
    ax1.set_xticklabels('')
    ax2.set_yscale('log'), ax2.set_ylim((1e-2,2))
    ax2.set_xscale(xscale), ax2.set_xlim((1,tmax))
    ax1.minorticks_on()
    ax4.set_xlim((0,900)), ax4.set_ylim((0,60))
    ax4.set_xticks(np.arange(0,901,300))
    ax4.set_xticklabels(np.arange(0,901,300), fontsize=11)
    ax4.set_yticks(np.arange(0,65,20))
    ax4.set_yticklabels(np.arange(0,65,20), fontsize=11)
    ax4.minorticks_on()
    if pltflag in [0,1]:
        ax1.set_ylabel('Total number of\nplanet detections', fontsize=13)
        #ax1.set_yticks(np.arange(0,51,10))
        #ax1.set_yticklabels(np.arange(0,51,10), fontsize=13)
        ax2.set_xlabel('Cumulative observing time [hours]', fontsize=13)
        ax2.set_ylabel('dN / dt\n[detections / hour]', fontsize=12, labelpad=0)
        ax1.set_title('TOIs spanning the radius valley', fontsize=13, y=1.14)
    else:
        ax1.set_yticklabels('')
        ax2.set_xticklabels('')
        ax2.set_yticklabels('')

    # second derivative axis
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    if pltflag in [0,1]:
        ax3.set_ylabel('dt / dN\n[hours / detection]', fontsize=12, labelpad=0)
        ax3.set_yticks([100,10,1])
        ax3.set_yticklabels(['100','10','1'], fontsize=13)
    else:
        ax3.set_yticklabels('')

    # second x-axis
    ax5 = ax4.twiny()
    ax5.set_xlim(tuple(np.ascontiguousarray(ax4.get_xlim())/7.))
    if pltflag in [0,1]:
        ax5.set_xticks(np.arange(0,130,30))
        ax5.set_xticklabels(np.arange(0,130,30), fontsize=11,
                            verticalalignment='baseline')
        #ax1.text(.25, .94, '[nights]', fontsize=10, transform=ax1.transAxes)
        #ax1.text(.26, .37, '[hours]', fontsize=10, transform=ax1.transAxes)
        ax4.set_xlabel('[hours]', fontsize=10)
        ax5.set_xlabel('[nights]', fontsize=10)
        ax5.set_xlim((0,9e2/7))
        ax5.minorticks_on()
    else:
        ax5.set_yticklabels('')

    # second x-axis
    ax6 = ax1.twiny()
    ax6.set_xlim(tuple(np.ascontiguousarray(ax1.get_xlim())/7.))
    if pltflag in [0,1]:
        ax6.set_xticks(10**np.arange(0,5))
        ax6.set_xticklabels(['10$^{%i}$'%i for i in np.arange(0,5)],
                            fontsize=12, va='baseline')
        ax6.set_xlabel('Cumulative observing time [nights]', fontsize=11)
        ax6.set_xlim((1/7.,tmax/7.)), ax6.set_xscale('log')
        ax6.minorticks_on()
    else:
        ax6.set_yticklabels('')

    #ax1.minorticks_on()
    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
        ax4.axis('off')
        ax5.axis('off')
        ax6.axis('off')
    fig.subplots_adjust(left=.19, bottom=.09, top=.89, right=.84, hspace=0)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        transparent = False if pltflag in [0,1] else True
        GPlabel = 'GP' if GP else ''
        plt.savefig('plots/cumulativetobs%s_radvalley_%s.png'%(GPlabel,label),
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')

    return sortH, sortN  # return indices sorted by tobsN for histogram if desired


def plot_cumulative_detections_v_tobs_transmission(self, pltt=True, label=False,
                                                   tmax=7e5, Nrand=10, seed=0,
                                                   pltflag=0, GP=True,
                                                   sigmarho=3., sigma=5.):
    '''0=full, 1=background, 2=opt curves, 3=nir curves'''
    fig = plt.figure(figsize=(5.1,6.3))# 5.2
    gs = gridspec.GridSpec(10,1)
    ax1 = plt.subplot(gs[:7,0])
    ax2 = plt.subplot(gs[7:,0])
    ax4 = fig.add_axes([.28, .58, .28, .24])
    np.random.seed(int(seed))
    _, snr = compute_transmission_SNR(self)
    g = snr >= 10
    #_,corr = _get_tobs_scaling(sigmarho, self.starnums_med[g], self.mps_med[g])
    #print '%i planets'%np.where(np.isfinite(corr))[0].size
    corr = (.327/.189)**2 if sigma == 5 else 1
    
    # NIRPS total observing time
    if pltflag in [0,3]:
        tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
        tobs = np.sort(tobsN[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
	ax2.text(1500, .4, 'Near-IR (YJH)', color=colN, fontsize=13, 
		 weight='semibold', verticalalignment='center')
        ax2.plot([3.5e2,1e3], np.repeat(.4,2), '-', c=colN)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '-', c=colN)
        ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
            tobs = tobsN[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    # HARPS total observing time
    if pltflag in [0,2]:
        tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
        tobs = np.sort(tobsH[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax2.text(1500, 1, 'Optical (BVR)', color=colH, fontsize=13,
                 weight='semibold', verticalalignment='center')
        ax2.plot([3.5e2,1e3], np.repeat(1,2), '--', c=colH)        
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
        # plot linear subplot
        ax4.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
    if pltflag in [0,1]:
        threshold = 1./20  # detections / hour
        ax2.plot([1,tmax], np.repeat(threshold,2), 'k--')
        ax2.text(5, threshold*.88,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='top', fontsize=10)
    if pltflag in [0,2]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
            tobs = tobsH[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            ax4.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)
            
    ax1.set_xscale('log'),
    ax1.set_xlim((1,tmax)), ax1.set_ylim((0,g.sum()))
    ax1.set_xticklabels('')
    ax1.minorticks_on()
    ax2.set_yscale('log'), ax2.set_ylim((1e-2,2))
    ax2.set_xscale('log'),
    ax2.set_xlim((1,tmax))
    if pltflag in [0,1]:
        ax1.set_ylabel('Total number of\nplanet detections', fontsize=13)
        ax2.set_xlabel('Cumulative observing time [hours]', fontsize=13)
        ax2.set_ylabel('dN / dt\n[detections / hour]', fontsize=12, labelpad=0)
        ax1.set_title('TOIs for JWST transmission follow-up', fontsize=13, y=1.14)
    else:
        ax1.set_yticklabels('')
        ax2.set_xticklabels('')
        ax2.set_yticklabels('')

    # second derivative axis
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    if pltflag in [0,1]:
        ax3.set_ylabel('dt / dN\n[hours / detection]', fontsize=12, labelpad=1)
        ax3.set_yticks([100,10,1])
        ax3.set_yticklabels(['100','10','1'], fontsize=13)
    else:
        ax3.set_yticklabels('')

    # add inset
    ax4.set_xlim((0,28e2)), ax4.set_ylim((0,240))
    ax4.set_xticks(np.arange(0,3000,1e3))
    ax4.set_xticklabels(np.arange(0,3000,1e3,dtype=int), fontsize=11)
    ax4.set_yticks(np.arange(0,210,100))
    ax4.set_yticklabels(np.arange(0,210,100,dtype=int), fontsize=11)
    ax4.minorticks_on()

    # second x-axis
    ax5 = ax4.twiny()
    ax5.set_xlim(tuple(np.ascontiguousarray(ax4.get_xlim())/7.))
    if pltflag in [0,1]:
        ax5.set_xticks(np.arange(0,401,100))
        ax5.set_xticklabels(np.arange(0,401,100), fontsize=11,
                            verticalalignment='baseline')
        ax4.set_xlabel('[hours]', fontsize=10, labelpad=0)
        ax5.set_xlabel('[nights]', fontsize=10)
        ax5.set_xlim((0,28e2/7))
        ax5.minorticks_on()
    else:
        ax5.set_yticklabels('')

    # second x-axis
    ax6 = ax1.twiny()
    ax6.set_xlim(tuple(np.ascontiguousarray(ax1.get_xlim())/7.))
    if pltflag in [0,1]:
        ax6.set_xticks(10**np.arange(0,5))
        ax6.set_xticklabels(['10$^{%i}$'%i for i in np.arange(0,5)],
                            fontsize=12, va='baseline')
        ax6.set_xlabel('Cumulative observing time [nights]', fontsize=11)
        ax6.set_xlim((1/7.,tmax/7.)), ax6.set_xscale('log')
        ax6.minorticks_on()
    else:
        ax6.set_yticklabels('')
    
    ax1.minorticks_on()
    ax2.minorticks_on()    
    ax3.minorticks_on()    
    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
        ax4.axis('off')
        ax5.axis('off')
        ax6.axis('off')
    fig.subplots_adjust(left=.19, bottom=.09, top=.89, right=.84, hspace=0)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        transparent = False if pltflag in [0,1] else True
        GPlabel = 'GP' if GP else ''
        plt.savefig('plots/cumulativetobs%s_transmission_%s.png'%(GPlabel,
                                                                  label),
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')


def plot_transmission_scatter(self, pltt=True, label=False, sigmarho=3):
    g = (self.Jmags_med >= 6.1) & (self.Jmags_med <= 10.73) & \
        (self.transmission_ppm > 30)
    _,corr = _get_tobs_scaling(sigmarho, self.starnums_med, self.mps_med)
    jsN = get_planets(self.tobsGPs_med_N*corr, g, 120)
    jsH = get_planets(self.tobsGPs_med_H*corr, g, 120)

    fig = plt.figure(figsize=(5.7,4.7))
    gs = gridspec.GridSpec(8,8)
    ax1 = plt.subplot(gs[1:,:-1])
    ax2 = plt.subplot(gs[0,:-1]) 
    ax3 = plt.subplot(gs[1:,-1]) 
    img = ax1.scatter(self.rps_med[jsN], self.Tps_med[jsN], marker='o', s=20,
                      c=self.transmission_ppm[jsN],
                      norm=LogNorm(vmin=30, vmax=1e3),
		      cmap=_truncate_colormap(plt.get_cmap('rainbow'),0,1))
    cbar_axes = fig.add_axes([.08, .08, .84, .032])
    cbar = fig.colorbar(img, cax=cbar_axes, orientation='horizontal')
    cbar.set_label('Expected transmission depth [ppm]', fontsize=12)

    ax1.scatter(self.rps_med[jsH], self.Tps_med[jsH], marker='d', s=25,
                c=self.transmission_ppm[jsH],
                norm=LogNorm(vmin=30, vmax=1e3),
		cmap=_truncate_colormap(plt.get_cmap('rainbow'),0,1))

    ax1.set_xscale('log'), ax1.set_xlim((1,23)), ax1.set_ylim((200,4e3))

    
    if pltt:
        plt.show()
    plt.close('all')
    
    
    
def plot_cumulative_detections_v_tobs_CVZ(self, pltt=True, label=False,
                                          Nrand=10, seed=0, pltflag=0,
                                          tmax=2e5, sigmarho=3., xscale='log'):
    '''0=all, 1=bkgd, 2=optical curves, 3=nIR curves''' 
    fig = plt.figure(figsize=(5.5,5.4))# 5.2
    gs = gridspec.GridSpec(10,1)
    ax1 = plt.subplot(gs[:7,0])
    ax2 = plt.subplot(gs[7:,0])
    #ax1 = fig.add_subplot(211)
    #ax2 = fig.add_subplot(212)
    np.random.seed(int(seed))
    inCVZ = is_star_in_CVZ(self.ras_med, self.decs_med)
    g = inCVZ == 1

    # get scaling to 3sigma density detection
    _,scaling = _get_tobs_scaling(sigmarho, self.starnums_med, self.mps_med)
    
    # HARPS total observing time
    if pltflag in [0,2]:
        tobs = np.sort(self.tobss_med_H[g] * scaling[g])
        tobs = np.append(0, np.cumsum(tobs))
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax2.text(.12, .28, 'Optical', color=colH, fontsize=11,
                 weight='semibold', transform=ax2.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
    if pltflag in [0,1]:
        threshold = 1./20  # detections / hour
        ax2.plot([1,tmax], np.repeat(threshold,2), 'k--')
        ax2.text(4e2, threshold*1.07,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='bottom', fontsize=10)
    if pltflag in [0,2]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = self.tobss_med_H[g][inds] * scaling[g][inds]
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)

    # NIRPS total observing time
    if pltflag in [0,3]:
        tobs = np.sort(self.tobss_med_N[g] * scaling[g])
        tobs = np.append(0, np.cumsum(tobs))[:51]
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
	ax2.text(.12, .17, 'Near-IR', color=colN, fontsize=11, 
		 weight='semibold', transform=ax2.transAxes)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '-', c=colN)
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = self.tobss_med_N[g][inds] * scaling[g][inds]
            tobs = np.append(0, np.cumsum(tobs))[:51]
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    ax1.set_xscale(xscale), ax1.set_xlim((1,tmax))
    ax1.set_ylim((-1,inCVZ.sum()))    
    ax1.set_xticklabels('')
    ax2.set_yscale('log'), ax2.set_ylim((8e-4,.8))
    ax2.set_xscale(xscale), ax2.set_xlim((1,tmax))
    ax1.minorticks_on(), ax2.minorticks_on()
    if pltflag in [0,1]:
        ax1.set_ylabel('Total number of\nplanet detections', fontsize=13)
        ax1.set_yticks(np.arange(0,31,10))
        ax1.set_yticklabels(np.arange(0,31,10))
        ax2.set_xlabel('Cumulative observing time [hours]', fontsize=13)
        ax2.set_ylabel('dN / dt\n[detections / hour]', fontsize=12, labelpad=0)
        ax1.set_title('TOIs in the JWST CVZ', fontsize=13)
    else:
        ax1.set_yticklabels('')
        ax2.set_xticklabels('')
        ax2.set_yticklabels('')
        
    # second derivative axis
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    if pltflag in [0,1]:
        ax3.set_ylabel('dt / dN\n[hours / detection]', fontsize=12, labelpad=4)
        #ax3.set_yticks(np.logspace(3,0,4))
        #ax3.set_yticklabels(['%i'%i for i in np.logspace(2,0,3)])
    else:
        ax3.set_yticklabels('')
    
    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
    fig.subplots_adjust(left=.18, bottom=.1, top=.96, right=.84, hspace=0)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        elif pltflag == 4:
            label = 'optrand'
        elif pltflag == 5:
            label = 'nirrand'
        elif pltflag == 6:
            label = 'labels'
        elif pltflag == 7:
            label = 'threshold'
        transparent = False if pltflag in [0,1] else True 
        plt.savefig('plots/cumulativetobs_CVZ_%s.png'%label,
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')


def _get_tobs_scaling(sigmarho, starnums, mps, sigP=5e-5,
                      fracsigMs=.1, fracsigRs=.1):
    '''
    Compute the scaling of tobs from a 3sig mass detection to a sigmarho 
    density detection.
    '''
    g = np.argsort(starnums)
    # Get rp uncertainty
    _,_,rp,P,_,K,Rs,_,_,_,_,_,_,_,logsigV,_,_ = get_TESS_data()
    depth    = compute_depth(rp, Rs)
    Gamma    = compute_Gamma()
    T        = compute_transit_duration(rp, P, K, Rs)
    Q        = compute_Q(Gamma, T, depth, logsigV)
    sigdepth = compute_sigdepth(depth, Q)
    depth = unp.uarray(depth, sigdepth)
    Rs2 = rvs.m2Rearth(rvs.Rsun2m(unp.uarray(Rs, fracsigRs*Rs)))
    rp2 = unp.sqrt(depth)*Rs2
    sigrp = unp.std_devs(rp2)

    # Get mass uncertainty
    fracsigrho = 1./sigmarho
    sigmp = mps * np.sqrt(fracsigrho**2 - (3*sigrp[g]/rp[g])**2)
    if not np.any(np.isfinite(sigmp)):
        raise ValueError('Cannot measure the bulk density this precisely ' + \
                         'for any TESS planet.')

    # Get K uncertainty
    fracsigK = np.sqrt((sigmp/mps)**2 - (sigP/(3*P[g]))**2 - (2*fracsigMs/3)**2)

    # return scaling of tobs from 3sigma mass detection
    return fracsigK, (.327 / fracsigK)**2
    
    

def plot_cumulative_detections_v_tobs_HZ(self, pltt=True, label=False,
                                         tmax=2e4, Nrand=10, seed=0, pltflag=0,
                                         GP=True, sigma=3.):
    '''0=full, 1=background, 2=opt curves, 3=nir curves'''
    fig = plt.figure(figsize=(5.1,6.3))# 5.2
    gs = gridspec.GridSpec(10,1)
    ax1 = plt.subplot(gs[:7,0])
    ax2 = plt.subplot(gs[7:,0])
    #ax4 = fig.add_axes([.25, .62, .25, .2])
    np.random.seed(int(seed))
    corr = (.327/.189)**2 if sigma == 5 else 1.
    g = (self.HZflags_med == 1) & (self.rps_med <= 2)

    # NIRPS total observing time
    if pltflag in [0,3]:
        tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
        tobs = np.sort(tobsN[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
	ax2.text(2e2, 3e-1, 'Near-IR (YJH)', color=colN, fontsize=13, 
		 weight='semibold')
        ax2.plot([70,160], np.repeat(4e-1,2), '-', c=colN)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '-', c=colN)
        #ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobsN = self.tobsGPs_med_N if GP else self.tobss_med_N
            tobs = tobsN[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    # HARPS total observing time
    if pltflag in [0,2]:
        tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
        tobs = np.sort(tobsH[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax2.text(2e2, 7e-1, 'Optical (BVR)', color=colH, fontsize=13,
                 weight='semibold')
        ax2.plot([70,170], np.repeat(9e-1,2), '--', c=colH)        
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
        # plot linear subplot
        #ax4.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
    if pltflag in [0,1]:
        threshold = 1./20  # detections / hour
        ax2.plot([1,tmax], np.repeat(threshold,2), 'k--')
        ax2.text(130, threshold*1.06,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='bottom', fontsize=10)
    if pltflag in [0,2]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobsH = self.tobsGPs_med_H if GP else self.tobss_med_H
            tobs = tobsH[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #ax4.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)
            
    ax1.set_xscale('log'),
    ax1.set_xlim((1,tmax)), ax1.set_ylim((0,g.sum()))
    ax1.set_xticklabels('')
    ax2.set_yscale('log'), ax2.set_ylim((1e-2,2))
    ax2.set_xscale('log'),
    ax2.set_xlim((1,tmax))
    if pltflag in [0,1]:
        ax1.set_ylabel('Total number of\nplanet detections', fontsize=13)
        #ax1.set_yticks(np.arange(0,51,10))
        #ax1.set_yticklabels(np.arange(0,51,10), fontsize=13)
        ax2.set_xlabel('Cumulative observing time [hours]', fontsize=13)
        ax2.set_ylabel('dN / dt\n[detections / hour]', fontsize=12, labelpad=0)
        ax1.set_title('Temperate Earths & super-Earths', fontsize=13, y=1.14)
    else:
        ax1.set_yticklabels('')
        ax2.set_xticklabels('')
        ax2.set_yticklabels('')

    ax1.minorticks_on()
    ax2.minorticks_on()
    #ax4.set_xlim((0,70)), ax4.set_ylim((0,6))
    #ax4.set_xticks(np.arange(0,70,20))
    #ax4.set_xticklabels(np.arange(0,70,20), fontsize=11)
    #ax4.set_yticks(np.arange(0,7,2))
    #ax4.set_yticklabels(np.arange(0,7,2,dtype=int), fontsize=11)
    #ax4.minorticks_on()

    # second derivative axis
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    if pltflag in [0,1]:
        ax3.set_ylabel('dt / dN\n[hours / detection]', fontsize=12, labelpad=1)
        ax3.set_yticks([100,10,1])
        ax3.set_yticklabels(['100','10','1'], fontsize=13)
    else:
        ax3.set_yticklabels('')

    # second x-axis
    #ax5 = ax4.twiny()
    #ax5.set_xlim(tuple(np.ascontiguousarray(ax4.get_xlim())/7.))
    #if pltflag in [0,1]:
        #ax5.set_xticks(np.arange(0,11,5))
        #ax5.set_xticklabels(np.arange(0,11,5), fontsize=11,
        #                    verticalalignment='baseline')
        #ax4.set_xlabel('[hours]', fontsize=10)
        #ax5.set_xlabel('[nights]', fontsize=10)
        #ax5.set_xlim((0,70/7.))
        #ax5.minorticks_on()
    #else:
        #ax5.set_yticklabels('')

    # second x-axis
    ax6 = ax1.twiny()
    ax6.set_xlim(tuple(np.ascontiguousarray(ax1.get_xlim())/7.))
    if pltflag in [0,1]:
        ax6.set_xticks(10**np.arange(0,5))
        ax6.set_xticklabels(['10$^{%i}$'%i for i in np.arange(0,5)],
                            fontsize=12, va='baseline')
        ax6.set_xlabel('Cumulative observing time [nights]', fontsize=11)
        ax6.set_xlim((1/7.,tmax/7.)), ax6.set_xscale('log')
        ax6.minorticks_on()
    else:
        ax6.set_yticklabels('')

    ax1.minorticks_on()
    ax2.minorticks_on()    
    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
        #ax4.axis('off')
        #ax5.axis('off')
        ax6.axis('off')
    fig.subplots_adjust(left=.19, bottom=.09, top=.89, right=.84, hspace=0)
    if label:
        if pltflag == 0:
            label = 'full'
        elif pltflag == 1:
            label = 'bkgd'
        elif pltflag == 2:
            label = 'opt'
        elif pltflag == 3:
            label = 'nir'
        transparent = False if pltflag in [0,1] else True
        GPlabel = 'GP' if GP else ''
        plt.savefig('plots/cumulativetobs%s_HZ_%s.png'%(GPlabel, label),
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')
    

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


def plot_identifying_best_50(self, s=20, pltt=True, label=False, pltflag=0,
                             GP=True, sigma=5.):
    '''
    Plot variables to try and guess which parameters best tell which small 
    planets we can observe most efficiently.
    e.g. self_xarr = self.Vmags_med
    0: full, 1:bkgd, 2:best, 3:regions
    '''
    tobsH = self.tobsGPs_med_H/7 if GP else self.tobss_med_H/7
    tobsN = self.tobsGPs_med_N/7 if GP else self.tobss_med_N/7
    scale = (.327/.189)**2 if sigma == 5 else 1.
    tobsN *= scale
    tobsH *= scale
    transit_quantity = lambda P, rp: rp/P**(1./3)
    TQs = transit_quantity(self.Ps_med, self.rps_med)
    
    fig = plt.figure(figsize=(4,6.2))  # 5.7,5.1
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    colmap = _truncate_colormap(plt.get_cmap('rainbow'),0,1)
    # set colorbar
    vmin, vmax = 1, 1e3
    if pltflag in [0,1]:
        img = ax1.scatter(self.Vmags_med, TQs,
                          c=tobsH, cmap=plt.get_cmap(colmap), s=0,
                          norm=colors.LogNorm(vmin=vmin, vmax=vmax))
        # add transluscent points
        ax1.scatter(self.Vmags_med, TQs, c=tobsH,
                    facecolors='none', cmap=plt.get_cmap(colmap), alpha=1, s=s,
                    norm=colors.LogNorm(vmin=vmin,vmax=vmax))
        cbar_axes = fig.add_axes([.08, .08, .84, .032])
        cbar = fig.colorbar(img, cax=cbar_axes, orientation='horizontal')
        cbar.set_label('Total observing time per TOI [nights]', fontsize=12)
    
    # Get 50 best
    if pltflag in [0,2]:
        sort = get_planets(tobsH, self.rps_med<4, N=50)
	print self.Vmags_med[sort].size
        ax1.scatter(self.Vmags_med[sort], TQs[sort],
                    facecolor='k', marker='^', edgecolor='k', s=s+10)
    
    # fill 'good' region
    if pltflag in [0,3]:
        m, b = get_best_fraction(self.Vmags_med, TQs, tobsH, self.rps_med)
        Vmax = 10.7
        print 'V slope and intercept: ', m, b
        print 'Vmax: ', Vmax
        line = lambda x: m*x + b
        ax1.fill_between([3.5,Vmax], [line(3.5),line(Vmax)], 30, color='k',
                         alpha=.15)
        ax1.plot(np.repeat(Vmax,2), [line(Vmax),30], 'k--', lw=1.7)
        ax1.plot([3.5,Vmax], [line(3.5),line(Vmax)], 'k--', lw=1.7)
    
    ax1.set_xlabel('V', labelpad=1, fontsize=14, style='italic')
    ax1.set_ylabel('$\Omega \equiv$ r$_p$ / P$^{1/3}$',
                   fontsize=12)
    ax1.set_xlim((3.5,20)), ax1.set_ylim((0,4))
    ax1.set_yticks(np.arange(0,5)), ax1.set_yticklabels(np.arange(0,5))
    ax1.minorticks_on()

    # Plot Jmag
    if pltflag in [0,1]:
        # set colorbar
        ax2.scatter(self.Jmags_med, TQs, c=tobsN, cmap=plt.get_cmap(colmap),
                    s=0, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
        # add transluscent points
        ax2.scatter(self.Jmags_med, TQs, c=tobsN,
                    facecolors='none', cmap=plt.get_cmap(colmap), alpha=1, s=s,
                    norm=colors.LogNorm(vmin=vmin,vmax=vmax))

    # Get 50 best
    if pltflag in [0,2]:
        sort = get_planets(tobsN, self.rps_med<4, N=50)
	print self.Jmags_med[sort].size
        ax2.scatter(self.Jmags_med[sort], TQs[sort], facecolor='k',
                    edgecolor='k', marker='^', s=s+10)

    # fill 'good' region
    if pltflag in [0,3]:
        m, b = get_best_fraction(self.Jmags_med, TQs, tobsN, self.rps_med)
        Jmax = 11.7
        print 'J slope and intercept: ', m, b
        print 'Jmax: ', Jmax
        line = lambda x: m*x + b
        ax2.fill_between([3.5,Jmax], [line(3.5),line(Jmax)], 30, color='k',
                         alpha=.15)
        ax2.plot(np.repeat(Jmax,2), [line(Jmax),30], 'k--', lw=1.7)
        ax2.plot([3.5,Jmax], [line(3.5),line(Jmax)], 'k--', lw=1.7)

    ax2.set_xlabel('J', labelpad=1, fontsize=14, style='italic')
    ax2.set_ylabel('$\Omega \equiv$ r$_p$ / P$^{1/3}$',
                   fontsize=12)
    ax2.set_xlim((3.5,15)), ax2.set_ylim((0,4))
    ax2.set_yticks(np.arange(0,5)), ax2.set_yticklabels(np.arange(0,5))
    ax2.minorticks_on()

    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
    fig.subplots_adjust(top=.98, bottom=.19, left=.14, right=.94, hspace=.21)
    if label:
        if pltflag == 0:
            label = 'full'
            transparent = False
        elif pltflag == 1:
            label = 'bkgd'
            transparent = False
        elif pltflag == 2:
            label = 'best'
            transparent = True
        elif pltflag == 3:
            label = 'region'
            transparent = True
        plt.savefig('plots/identify50_%s.png'%label, transparent=transparent)
    if pltt:
        plt.show()
    plt.close('all')


def get_best_fraction(mags, TQs, tobss, rps):
    #sort = TQs <= 5
    sort = get_planets(tobss, rps<4, N=50)
    mags, TQs, tobss = mags[sort], TQs[sort], tobss[sort]
    w = np.ones(mags.size)
    w[abs(TQs-np.mean(TQs)) > 3*np.std(TQs)] = 0.  # clip outliers
    m, b =  np.polyfit(mags, TQs, 1, w=w/w.sum())
    line = lambda x: m*x + b
    g = TQs > line(mags)
    while g.sum() < TQs.size - .05*TQs.size:
        b -= 1e-2
        g = TQs > line(mags)
    return m, b


def MAD(xarr):
    return np.median(abs(xarr - np.median(xarr)))


def _get_transmission_depth(self):
    #S = self.Fs_med*1367
    #Teq = (.25*S(1.-albedo)/sigmaSB)**(.25)
    smas = rvs.semimajoraxis(self.Ps_med, self.Mss_med, self.mps_med)
    self.Tps_med = self.Teffs_med * np.sqrt(rvs.Rsun2m(self.Rss_med)/ \
                                        (2*rvs.AU2m(smas)))
    self.mus_med = np.repeat(2., self.nstars)
    self.mus_med[self.rps_med <= 2] = 29.
    self.transmission_ppm = rvs.transmission_spectroscopy_depth(self.Rss_med,
                                                                self.mps_med,
                                                                self.rps_med,
                                                                self.Tps_med,
                                                                self.mus_med)


def plot_individual_Nrv_histogram(self, starnum, nbins=15,
                                  label=False, pltt=True):
    assert 0 <= starnum <= 1983
    fig = plt.figure(figsize=(5,6))
    ax = fig.add_subplot(111)
    scale = (.327/.189)**2
    Nrvs = self.NrvGPs*scale
    Nrvs[Nrvs < 10] = 10.
    tobs = self.texps*Nrvs / 60.

    # plot optical
    g = (self.starnums == starnum) & (self.spectrographs == 'H')
    assert g.sum() == 100
    ax.hist(Nrvs[g], bins=np.logspace(0,np.log10(Nrvs[g].max()),nbins),
            histtype='step', color=colH, lw=3, label='Optical')
    ax.axvline(np.median(Nrvs[g]), ls='--', lw=.9, color=colH)
    ax.set_xlabel('N$_{RV}$ (5$\sigma$)'), ax.set_xscale('log')
    ax.set_xlim((1,3e3))

    # plot nearIR
    gN = (self.starnums == starnum) & (self.spectrographs == 'N')
    assert g.sum() == 100
    ax.hist(Nrvs[gN], bins=np.logspace(0,np.log10(Nrvs[gN].max()),nbins),
            histtype='step', color=colN, lw=3, label='Near-IR')
    ax.axvline(np.median(Nrvs[gN]), ls='--', lw=.9, color=colN)
    ax.legend()
    ax.set_title('TOI = %.4d'%starnum, weight='semibold')
    
    # first compute transmission depth
    sma = rvs.AU2m(rvs.semimajoraxis(self.Ps[g][0], self.Mss[g][0],
                                     self.mps[g][0]))
    Teq = self.Teffs[g][0] * np.sqrt(rvs.Rsun2m(self.Rss[g][0]) / (2*sma))
    mu = 30. if self.rps[g][0] <=2 else 2.
    transmission_ppm = rvs.transmission_spectroscopy_depth(self.Rss[g][0],
                                                           self.mps[g][0],
                                                           self.rps[g][0],
                                                           Teq, mu)

    if np.median(tobs[g]) < np.median(tobs[gN]):
        weightH, weightN = 'semibold', 'normal'
    else:
        weightH, weightN =  'normal', 'semibold'

    # Stellar/planet properties
    gmed = self.starnums_med == starnum
    ax.text(0, 1.45, 'J = %.2f'%self.Jmags_med[gmed][0],
            transform=ax.transAxes)    
    ax.text(0, 1.4, 'M$_s$ = %.2f M$_{\odot}$'%self.Mss[g][0],
            transform=ax.transAxes)
    ax.text(0, 1.35, 'Teff = %i K'%self.Teffs[g][0], transform=ax.transAxes)
    ax.text(0, 1.3, '$\sigma_{activity}$ = %.1f m/s'%self.sigmaRV_acts_med_H[gmed], transform=ax.transAxes)
    ax.text(0, 1.25, 't$_{exp,opt}$ = %.1f min'%self.texps_med_H[gmed],
            transform=ax.transAxes, color=colH, weight=weightH)
    ax.text(0,1.2,'med($\sigma_{RV,opt}$) = %.1f m/s'%self.sigmaRV_phot_med_H[gmed],
            transform=ax.transAxes, color=colH, weight=weightH)
    ax.text(0, 1.15, 't$_{exp,nIR}$ = %.1f min'%self.texps_med_N[gmed],
            transform=ax.transAxes, color=colN, weight=weightN)
    ax.text(0,1.1,'med($\sigma_{RV,nIR}$) = %.1f m/s'%self.sigmaRV_phot_med_N[gmed], transform=ax.transAxes, color=colN, weight=weightN)

    ax.text(.5, 1.45, 'P = %.2f days'%self.Ps[g][0], transform=ax.transAxes)
    ax.text(.5, 1.4, 'r$_p$ = %.2f R$_{\oplus}$'%self.rps[g][0],
            transform=ax.transAxes)
    ax.text(.5, 1.35, 'K = %.2f m/s ($\sigma_K$ = %.2f m/s)'%(self.Ks[g][0],
                                                              self.Ks[g][0]/5.),
            transform=ax.transAxes)
    ax.text(.5, 1.3, '$\delta F/F$ (transmission) = %i ppm'%transmission_ppm,
            transform=ax.transAxes)
    ax.text(.5, 1.25, 'med(N$_{RV,opt}$) = %.1f'%np.median(Nrvs[g]),
            transform=ax.transAxes, color=colH, weight=weightH)
    ax.text(.5, 1.2, 'med(t$_{obs,opt}$) = %.1f hours'%np.median(tobs[g]),
            transform=ax.transAxes, color=colH, weight=weightH)
    ax.text(.5, 1.15, 'med(N$_{RV,nIR}$) = %.1f'%np.median(Nrvs[gN]),
            transform=ax.transAxes, color=colN, weight=weightN)
    ax.text(.5, 1.1, 'med(t$_{obs,nIR}$) = %.1f hours'%np.median(tobs[gN]),
            transform=ax.transAxes, color=colN, weight=weightN)
    
    fig.subplots_adjust(bottom=.1, top=.7)
    if label:
        plt.savefig('plots/TOIhistograms/TOIhist_%.4d.png'%starnum)
    if pltt:
        plt.show()
    plt.close('all')



def plot_WP2_Nrvs_histograms(self, nbins=15, label=False, pltt=True):
    inds = np.array([  12,   47,   92,  138,  167,  200,  219,  356,  387,  423,  431,
                       446,  470,  487,  515,  526,  544,  708,  815,  817,  832,  840,
                       844,  857,  859,  876,  877,  933,  935,  944,  947,  969,  992,
                       1038, 1042, 1055, 1058, 1071, 1072, 1075, 1120, 1161, 1164, 1165,
                       1166, 1188, 1210, 1213, 1221, 1223, 1254, 1259, 1292, 1295, 1307,
                       1321, 1337, 1392, 1405, 1412, 1420, 1438, 1451, 1475, 1481, 1498,
                       1513, 1522, 1526, 1560, 1565, 1607, 1666, 1799, 1806, 1834, 1860,
                       1873, 1929, 1934, 1944, 1956, 1967, 1976])

    # plot Nrvs
    fig = plt.figure(0, figsize=(6,5))
    ax = fig.add_subplot(111)
    g = np.in1d(self.starnums_med, inds)
    scale = (.327/.189)**2
    ax.hist(self.NrvGPs_med_N[g]*scale,
            bins=np.logspace(0,np.log10(self.NrvGPs_med_N[g].max()*scale),nbins), histtype='step', lw=3, color=colN, label='Near-IR')
    ax.axvline(np.median(self.NrvGPs_med_N[g]*scale), ls='--', lw=.9,
               color=colN)
    ax.text(.65, .85, 'median = %.2f'%np.median(self.NrvGPs_med_N[g]*scale),
            transform=ax.transAxes)
    ax.set_xlabel('median(N$_{RV}$) (5$\sigma$)'), ax.set_xscale('log')
    ax.set_xlim((1,3e3)), ax.set_title('%i WP2 TOIs'%inds.size)

    if label:
        plt.savefig('plots/WP2Nrvs.png')
    if pltt:
        plt.show()
    plt.close('all')



def plot_WP2_texps_histograms(self, nbins=15, label=False, pltt=True):
    inds = np.array([  12,   47,   92,  138,  167,  200,  219,  356,  387,  423,  431,
                       446,  470,  487,  515,  526,  544,  708,  815,  817,  832,  840,
                       844,  857,  859,  876,  877,  933,  935,  944,  947,  969,  992,
                       1038, 1042, 1055, 1058, 1071, 1072, 1075, 1120, 1161, 1164, 1165,
                       1166, 1188, 1210, 1213, 1221, 1223, 1254, 1259, 1292, 1295, 1307,
                       1321, 1337, 1392, 1405, 1412, 1420, 1438, 1451, 1475, 1481, 1498,
                       1513, 1522, 1526, 1560, 1565, 1607, 1666, 1799, 1806, 1834, 1860,
                       1873, 1929, 1934, 1944, 1956, 1967, 1976])

    # plot exposure times
    fig = plt.figure(1, figsize=(6,5))
    ax = fig.add_subplot(111)
    g = np.in1d(self.starnums_med, inds)
    scale = (.327/.189)**2    
    ax.hist(self.texps_med_N[g], bins=np.linspace(10,60,nbins),
            histtype='step', lw=3, color=colN, label='Near-IR')
    ax.axvline(np.median(self.texps_med_N[g]), ls='--', lw=.9, color=colN)
    ax.text(.65, .85, 'median = %.2f'%np.median(self.texps_med_N[g]),
            transform=ax.transAxes)
    ax.set_xlabel('median(t$_{exp}$) [minutes]'), ax.set_xlim((0,60)),
    ax.set_title('%i WP2 TOIs'%inds.size)

    if label:
        plt.savefig('plots/WP2texps.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_WP2_tobs_histograms(self, nbins=15, label=False, pltt=True):
    inds = np.array([  12,   47,   92,  138,  167,  200,  219,  356,  387,  423,  431,
                       446,  470,  487,  515,  526,  544,  708,  815,  817,  832,  840,
                       844,  857,  859,  876,  877,  933,  935,  944,  947,  969,  992,
                       1038, 1042, 1055, 1058, 1071, 1072, 1075, 1120, 1161, 1164, 1165,
                       1166, 1188, 1210, 1213, 1221, 1223, 1254, 1259, 1292, 1295, 1307,
                       1321, 1337, 1392, 1405, 1412, 1420, 1438, 1451, 1475, 1481, 1498,
                       1513, 1522, 1526, 1560, 1565, 1607, 1666, 1799, 1806, 1834, 1860,
                       1873, 1929, 1934, 1944, 1956, 1967, 1976])
    
    # plot observing times
    fig = plt.figure(2, figsize=(6,5))
    ax = fig.add_subplot(111)
    g = np.in1d(self.starnums_med, inds)
    scale = (.327/.189)**2
    ax.hist(self.tobsGPs_med_N[g]*scale,
            bins=np.logspace(0,np.log10(self.tobsGPs_med_N[g].max()*scale),nbins), histtype='step', lw=3, color=colN, label='Near-IR')
    ax.axvline(np.median(self.tobsGPs_med_N[g]*scale), ls='--', lw=.9,
               color=colN)
    ax.text(.65, .85, 'median = %.2f'%np.median(self.tobsGPs_med_N[g]*scale),
            transform=ax.transAxes)
    ax.set_xlabel('median(t$_{obs}$) [hours]'), ax.set_xscale('log')
    ax.set_xlim((1,3e3)), ax.set_title('%i WP2 TOIs'%inds.size)

    if label:
        plt.savefig('plots/WP2tobs.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_sigRV_phot(self, pltt=True, label=False):
    fig = plt.figure(figsize=(4.5,3.7))
    ax = fig.add_subplot(111)

    bins = np.logspace(-.5,np.log10(30),12)
    ax.hist(np.array([self.sigmaRV_phot_med_H, self.sigmaRV_phot_med_N]).T, bins,
            histtype='bar', color=[colH,colN])
    ax.set_xscale('log'), ax.set_xlabel('$\sigma_{RV}$ [m/s]')
    ax.set_ylabel('Number of TOIs'), ax.minorticks_on()
    ax.set_xticks([.3,1,3,10,30])
    ax.set_xticklabels(['0.3','1','3','10','30'])
    
    ax.text(.08, .88, 'Optical (BVR)', weight='semibold', color=colH,
            transform=ax.transAxes, fontsize=12)
    ax.text(.08, .8, 'Near-IR (YJH)', weight='semibold', color=colN,
            transform=ax.transAxes, fontsize=12)

    fig.subplots_adjust(bottom=.14, top=.96, right=.96, left=.16)
    if label:
        plt.savefig('plots/sigRVhist.png')
    if pltt:
        plt.show()
    plt.close('all')
