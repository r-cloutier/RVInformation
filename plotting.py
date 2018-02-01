import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from uncertainties import unumpy as unp
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import matplotlib as mpl
import matplotlib.colors as colors
from CVZ import *
from get_tess_data import get_TESS_data
from analysis import *
from rvmodel import get_rv1
from compute_sigK_QPGP import compute_sigmaK_GP

#params = {'text.usetex': True, 
#          'text.latex.preamble': [r'\usepackage{cmbright}',
#                                  r'\usepackage{amsmath}']}
#plt.rcParams.update(params)

global colH, colN
colH, colN = '#4444ff', '#a30000'
    
mpl.rc('xtick', labelsize=13)
mpl.rc('ytick', labelsize=13)
mpl.rc('axes', titlepad=3)

global toverhead, Teffedges, Tefflabels, rpedges, rplabels, rplabels2
toverhead = 5.
Teffedges = np.array([25e2, 32e2, 38e2, 76e2, 12e3])
Tefflabels = ['mid-late M', 'early-mid M', 'FGK', 'BA']
rpedges = np.array([0, 1.25, 2, 4, 30])
rplabels = ['Earths','Super-Earths','Sub-Neptunes','Giants']
rplabels2 = ['$<$1.25 R$_{\oplus}$','1.25-2 R$_{\oplus}$',
                 '2-4 R$_{\oplus}$','$>$4 R$_{\oplus}$']

def _truncate_colormap(cmap, minval=0, maxval=1, n=100):
    return colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f},)'.format(n=cmap.name,a=minval,b=maxval),
        cmap(np.linspace(minval, maxval, n)))

def compute_Nrv(sigeff, sigK):
    sigeff = np.ascontiguousarray(sigeff).astype(float)
    return 2 * (sigeff / sigK)**2


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
    _,sigeffs, sigKs, Nrv_true, Ks, rednoiseflag = \
                                np.genfromtxt('Nrv_tests/FASTtests.dat').T
    detsig = Ks / sigKs
    c = detsig
    
    # Compute Nrv
    Nrv_calc = compute_Nrv(sigeffs, sigKs)

    # plotting
    fig = plt.figure(figsize=(3.8,4.5))
    ax = fig.add_subplot(111, aspect='equal')
    g = rednoiseflag == 0
    cax = ax.scatter(Nrv_true[g], Nrv_calc[g], edgecolors='none', marker='o',
                     c=c[g], s=60, label='white activity model',
                     norm=LogNorm(vmin=1, vmax=3e2),
		     cmap=_truncate_colormap(plt.get_cmap('Blues'),.3,1))
    cbar_axes = fig.add_axes([.08,.13,.84,.04])
    cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')
    cticklabels = ['']*22
    cticklabels[0], cticklabels[9], cticklabels[18] = '1', '10', '100'
    cbar.ax.set_xticklabels(cticklabels) 
    cbar.set_label(r'K / $\sigma_K$', labelpad=.1)
    
    #g = rednoiseflag == 1
    #ax.scatter(Nrv_true[g], Nrv_calc[g], edgecolors='none', marker='d',
    #           c=c[g], s=70, label='red activity model',
    #           norm=LogNorm(vmin=1, vmax=3e2),
    # 	       cmap=_truncate_colormap(plt.get_cmap('Blues'),.3,1))
    names = ['K218','LHS1140','Kep78HARPSN','Kep21HARPSN','CoRoT7']
    Ks = np.array([3.18, 5.34, 1.86, 2.12, 3.42])
    sigKtargets, Nrvtargets = np.array([.75, 1.1, .25, .66, .66]), \
                              np.array([75, 144, 109, 82, 71])
    sigeffs = sigKtargets * np.sqrt(Nrvtargets/2.)
    GPNrvs = np.zeros((sigeffs.size, 2))
    for i in range(sigeffs.size):
        sigKs = np.loadtxt('Nrv_tests/GPtest_%s.dat'%names[i])[:,-1]
        sigKs = sigKs[np.isfinite(sigKs)]
        Nrvs = 2 * (sigeffs[i] / sigKs)**2
        GPNrvs[i] = np.median(Nrvs), MAD(Nrvs)
        ax.text(Nrvtargets[i], GPNrvs[i,0]*.97, names[i], fontsize=6)
    ax.errorbar(Nrvtargets, GPNrvs[:,0], GPNrvs[:,1], fmt='ko', ms=0,
                elinewidth=.8, capsize=0)
    ax.scatter(Nrvtargets, GPNrvs[:,0], edgecolors='none', marker='d',
               c=Ks/sigKtargets, s=60, label='red activity model',
               norm=LogNorm(vmin=1, vmax=3e2),
	       cmap=_truncate_colormap(plt.get_cmap('Blues'),.3,1))
    
    arr = np.linspace(0, 3e2, 100)
    ax.plot(arr, arr, 'k--', lw=2)
    #for i in np.arange(.1,.5,.1):
    #	ax.fill_between(arr, (1-i)*arr, (1+i)*arr, alpha=.5)
    	#ax.plot(arr, (1+i)*arr, 'k--', lw=.6)
    	#ax.plot(arr, (1-i)*arr, 'k--', lw=.6)

    ax.set_xscale('log'), ax.set_yscale('log')
    ax.set_xlabel('Observed n$_{RV}$'), ax.set_ylabel('Calculated n$_{RV}$')
    ax.set_xlim((1, 3e2)), ax.set_ylim((1, 3e2))
    ax.set_xticks(np.logspace(0,2,3))
    ax.set_xticklabels(['%i'%i for i in np.logspace(0,2,3)])
    ax.set_yticks(np.logspace(0,2,3))
    ax.set_yticklabels(['%i'%i for i in np.logspace(0,2,3)])
    ax.minorticks_on()

    ax.legend(bbox_to_anchor=(.62,.95), fontsize=9,
              handletextpad=.01, labelspacing=.7)
    fig.subplots_adjust(bottom=.27, left=.19, right=.97, top=1)
    if label:
        plt.savefig('plots/Nrvcomparison.png')
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
def plot_compare_tobs_H_N(self, pltt=True, label=False):

    h = self.spectrographs == 'H'
    n = self.spectrographs == 'N'
    ratio = unp.uarray(self.tobss_med_H / self.tobss_med_N,
                       self.tobss_emed_H / self.tobss_emed_N)

    fig = plt.figure(figsize=(4.8,3.6))
    ax = fig.add_subplot(111)

    #ax.errorbar(self.Teffs_med, unp.nominal_values(ratio),
    #            unp.std_devs(ratio),
    #            fmt='ko', ms=1, elinewidth=.1)
    # point estimates
    img = ax.scatter(self.Teffs_med, unp.nominal_values(ratio), s=2, alpha=1,
                     c='k')#self.Bmags_med-self.Vmags_med)
    #cbar_axes = fig.add_axes([.9, .1, .05, .8])
    #cbar = fig.colorbar(img, cax=cbar_axes)
    #ax.fill_between(ax.get_xlim(), 8.1, 1./8.1, color='k',
    #                alpha=.1)
    ax.fill_between(ax.get_xlim(), 3, 1./3, color='k',
                    alpha=.1)

    ax.set_yscale('log')
    ax.set_ylim((1e-2, 1e2)), ax.set_xlim((25e2, 12e3))
    ax.set_ylabel('median(t$_{obs,opt}$) / median(t$_{obs,nIR}$)')
    ax.set_xlabel('Effective Temperature [K]')
    ax.minorticks_on()

    fig.subplots_adjust(left=.17, bottom=.14, right=.93, top=.91)
    if label:
        plt.savefig('plots/optnIRratio.png')
    if pltt:
        plt.show()
    plt.close('all')
        

def plot_cumulative_detections_v_tobs(self, pltt=True, label=False,
                                      harps=True, nirps=True, tmax=1e3,
                                      xscale='linear'):
    #tarr = np.arange(0,tmax,1)
    fig = plt.figure(figsize=(9.5,4.6))

    for i in range(rpedges.size-1):
	ax1 = fig.add_subplot(2, rpedges.size-1, i+1)

	# HARPS total observing time
	if harps:
            g = (self.rps_med >= rpedges[i]) & \
                (self.rps_med < rpedges[i+1])
            tobs = np.sort(self.tobss_med_H[g])
	    tobs = np.append(0, np.cumsum(tobs))
	    Ndet = np.arange(tobs.size)
	    g = tobs <= tmax+10
	    ax1.plot(tobs[g], Ndet[g], '--', c=colH, lw=2, drawstyle='steps')
            # plot derivative
	    if i == 0:
                ax1.plot([410,560], np.repeat(3.2,2), '--', lw=2, c=colH)
		ax1.text(6e2, 3.2, 'Optical', color=colH, fontsize=10, 
			 weight='semibold', verticalalignment='center')
            # plot derivative
            ax2 = fig.add_subplot(2, rpedges.size-1, i+1+rpedges.size-1)
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
            threshold = 1./20
            ax2.plot([0,tmax*1e3], np.repeat(threshold,2), 'k--', lw=.9)
            if i == 0:
                ax2.text(tmax*.2, threshold*1.03,
                         '%i hours / detection'%(1./threshold),
                         verticalalignment='bottom', fontsize=9)
                
        # NIRPS total observing time
	if nirps:
            g = (self.rps_med >= rpedges[i]) & (self.rps_med < rpedges[i+1])
            tobs = np.sort(self.tobss_med_N[g])
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
	    g = tobs <= tmax+10
            #_compute_curve_derivative(tobs[g], Ndet[g])
	    ax1.plot(tobs[g], Ndet[g], '-', c=colN, drawstyle='steps')
	    if i == 0:
                ax1.plot([410,550], np.repeat(1.2,2), '-', c=colN)
                ax1.text(6e2, 1.2, 'Near-IR', color=colN, fontsize=10, 
                         weight='semibold', verticalalignment='center')
            # plot derivative
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, '-', c=colN)

        ax1.set_title('%s\n%s'%(rplabels[i], rplabels2[i]), fontsize=10)
        ax2.set_xlabel('Cumulative observing\ntime [hours]', fontsize=9)
        ax1.set_xscale(xscale), ax1.set_xticklabels('')
        ax1.minorticks_on(), ax2.minorticks_on()
        ax2.set_xscale(xscale)
    	ax2.set_yscale('log'), ax2.set_ylim((1e-2,1.3))
    	ax1.set_xlim((0,tmax)), ax2.set_xlim((0,tmax))

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
        
    fig.subplots_adjust(left=.085, bottom=.125, top=.92, right=.925,
                        hspace=0, wspace=.22)
    if label:
        plt.savefig('plots/cumulativetobs_rp.png')
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
                                           harps=True, nirps=True, tmax=1e3,
                                           observebesttargetsfirst=True):
    fig = plt.figure(figsize=(9.5,4.6))

    for i in range(Teffedges.size-1):
	ax1 = fig.add_subplot(2, Teffedges.size-1, -1*i+4)

	# HARPS total observing time
	if harps:
            g = (self.Teffs_med >= Teffedges[i]) & \
                (self.Teffs_med < Teffedges[i+1])
            tobs = np.sort(self.tobss_med_H[g]) if observebesttargetsfirst \
                   else self.tobss_med_H[g]
	    tobs = np.append(0, np.cumsum(tobs))
	    Ndet = np.arange(tobs.size)
	    g = tobs <= tmax+100
	    ax1.plot(tobs[g], Ndet[g], '--', c=colH, lw=2, drawstyle='steps')
            # plot derivative
	    if i == 3:
                ax1.plot([410,560], np.repeat(25,2), '--', c=colH, lw=2)
		ax1.text(6e2, 25, 'Optical', color=colH, fontsize=10,
			 weight='semibold', verticalalignment='center')
            # plot derivative
            ax2 = fig.add_subplot(2, Teffedges.size-1, -1*i+4+Teffedges.size-1)
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
            threshold = 1./20
            ax2.plot([0,tmax], np.repeat(threshold,2), 'k--', lw=.9)
            if i == 3:
                ax2.text(tmax*.25, threshold*1.03,
                         '%i hours / detection'%(1./threshold),
                         verticalalignment='bottom', fontsize=9)
                
        # NIRPS total observing time
	if nirps:
            g = (self.Teffs_med >= Teffedges[i]) & \
                (self.Teffs_med < Teffedges[i+1])
            tobs = np.sort(self.tobss_med_N[g]) if observebesttargetsfirst \
                   else self.tobss_med_N[g]
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
	    g = tobs <= tmax+100
            #_compute_curve_derivative(tobs[g], Ndet[g])
	    ax1.plot(tobs[g], Ndet[g], '-', c=colN, drawstyle='steps')
	    if i == 3:
                ax1.plot([410,550], np.repeat(20,2), 'r-')
                ax1.text(6e2, 20, 'Near-IR', color=colN, fontsize=10, 
                         weight='semibold', verticalalignment='center')
            # plot derivative
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, '-', c=colN)
            
        ax1.set_title('%s\n%i $> T_{eff} \geq$ %i'%(Tefflabels[i],
                                     		      Teffedges[i+1],
                                     		      Teffedges[i]),
                      fontsize=10)
        ax2.set_xlabel('Cumulative observing\ntime [hours]', fontsize=9)
        ax1.set_xticklabels('')
        ax1.minorticks_on(), ax2.minorticks_on()
    	ax2.set_yscale('log'), ax2.set_ylim((1e-2,1.3))
    	ax1.set_xlim((0,tmax)), ax2.set_xlim((0,tmax))

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
            
    fig.subplots_adjust(left=.085, bottom=.125, top=.92, right=.925,
                        hspace=0, wspace=.23)
    if label:
        plt.savefig('plots/cumulativetobs_Teff.png')
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
                                                sigma=3., tmax=1e5,
                                                xscale='log'):
    '''0=all, 1=bkgd, 2=optical curves, 3=nIR curves''' 
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
        tobs = np.sort(self.tobss_med_H[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))[:51]
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax2.text(.7, .79, 'Optical', color=colH, fontsize=13,
                 weight='semibold', transform=ax2.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
    if pltflag in [0,1]:
        threshold = 1./20  # detections / hour
        ax2.plot([1,tmax], np.repeat(threshold,2), 'k--')
        ax2.text(4, threshold*1.06,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='bottom', fontsize=9)
    if pltflag in [0,2]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = self.tobss_med_H[g][inds]  * corr
            tobs = np.append(0, np.cumsum(tobs))[:51]
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)

    # NIRPS total observing time
    if pltflag in [0,3]:
        tobs = np.sort(self.tobss_med_N[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))[:51]
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
	ax2.text(.7, .63, 'Near-IR', color=colN, fontsize=13, 
		 weight='semibold', transform=ax2.transAxes)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '-', c=colN)
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = self.tobss_med_N[g][inds]  * corr
            tobs = np.append(0, np.cumsum(tobs))[:51]
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

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
        ax3.set_ylabel('dt / dN\n[hours / detection]', fontsize=12, labelpad=4)
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
        plt.savefig('plots/cumulativetobs_50random_%s.png'%label,
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
                                         tmax=1e6, Nrand=10, seed=0, pltflag=0,
                                         xscale='log', sigma=5.):
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
        tobs = np.sort(self.tobss_med_N[g] * corr)
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
            tobs = self.tobss_med_N[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    # HARPS total observing time
    if pltflag in [0,2]:
        tobs = np.sort(self.tobss_med_H[g] * corr)
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
            tobs = self.tobss_med_H[g][inds] * corr
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
    ax4.set_xlim((0,1e3)), ax4.set_ylim((0,50))
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
        plt.savefig('plots/cumulativetobs_MR_%s.png'%label,
                    transparent=transparent)
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
    scaling = _get_tobs_scaling(sigmarho, self.starnums_med, self.mps_med)
    
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
                 verticalalignment='bottom', fontsize=9)
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
    return (.327 / fracsigK)**2
    
    

def plot_cumulative_detections_v_tobs_HZ(self, pltt=True, label=False,
                                         tmax=2e5, Nrand=10, seed=0, pltflag=0,
                                         sigma=5.):
    '''0=full, 1=background, 2=opt curves, 3=nir curves'''
    fig = plt.figure(figsize=(5.5,5.4))# 5.2
    gs = gridspec.GridSpec(10,1)
    ax1 = plt.subplot(gs[:7,0])
    ax2 = plt.subplot(gs[7:,0])
    #ax4 = fig.add_axes([.24, .6, .32, .32])
    np.random.seed(int(seed))
    corr = (.327/.189)**2 if sigma == 5 else 1
    g = (self.HZflags_med == 1) & (self.rps_med < 2)

    # NIRPS total observing time
    if pltflag in [0,3]:
        tobs = np.sort(self.tobss_med_N[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
	ax1.text(.12, .65, 'Near-IR', color=colN, fontsize=13, 
		 weight='semibold', transform=ax1.transAxes)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '-', c=colN)
        #ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps')
    if pltflag in [0,3]:
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = self.tobss_med_N[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #ax4.plot(tobs, Ndet, '-', c=colN, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '-', c=colN, lw=.3)

    # HARPS total observing time
    if pltflag in [0,2]:
        tobs = np.sort(self.tobss_med_H[g] * corr)
        tobs = np.append(0, np.cumsum(tobs))
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
        ax1.text(.12, .73, 'Optical', color=colH, fontsize=13,
                 weight='semibold', transform=ax1.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, '--', c=colH, lw=2)
        # plot linear subplot
        #ax4.plot(tobs, Ndet, '--', c=colH, lw=2, drawstyle='steps')
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
            tobs = self.tobss_med_H[g][inds] * corr
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
            ax1.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #ax4.plot(tobs, Ndet, '--', c=colH, drawstyle='steps', lw=.3)
            #tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            #ax2.plot(tobs2, dNdt, '--', c=colH, lw=.3)
            
    ax1.set_xscale('log'),
    ax1.set_xlim((1,tmax)), ax1.set_ylim((0,g.sum()))
    ax1.set_xticklabels('')
    ax2.set_yscale('log'), ax2.set_ylim((8e-4,.2))
    ax2.set_xscale('log'),
    ax2.set_xlim((1,tmax))
    ax1.minorticks_on()
    '''ax4.set_xlim((0,1e3)), ax4.set_ylim((0,6))
    ax4.set_xticks(np.linspace(0,1e3,5))
    ax4.set_xticklabels(['0','','500','','1000'], fontsize=12)
    ax4.set_yticks(np.arange(0,7,2))
    ax4.set_yticklabels(np.arange(0,7,2,dtype=int), fontsize=12)
    ax4.minorticks_on()'''
    if pltflag in [0,1]:
        ax1.set_ylabel('Total number of\nplanet detections', fontsize=13)
        #ax1.set_yticks(np.arange(0,51,10))
        #ax1.set_yticklabels(np.arange(0,51,10), fontsize=13)
        ax2.set_xlabel('Cumulative observing time [hours]', fontsize=13)
        ax2.set_ylabel('dN / dt\n[detections / hour]', fontsize=12, labelpad=0)
        ax1.set_title('Potentially habitable TOIs', fontsize=12)
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
        #ax3.set_yticks([100,10,1])
        #ax3.set_yticklabels(['100','10','1'], fontsize=13)
    else:
        ax3.set_yticklabels('')
    
    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')
        #ax4.axis('off')
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
        plt.savefig('plots/cumulativetobs_HZ_%s.png'%label,
                    transparent=transparent)
    if pltt:
	plt.show()
    plt.close('all')
    

def _compute_curve_derivative(tobs, Ndet, sigma=5):
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


def plot_identifying_best_50(self, s=10, pltt=True, label=False, pltflag=0):
    '''
    Plot variables to try and guess which parameters best tell which small 
    planets we can observe most efficiently.
    e.g. self_xarr = self.Vmags_med
    0: full, 1:bkgd, 2:best, 3:regions
    '''
    g = (self.rps_med < 4)
    tobss_med = np.append(self.tobss_med_H[g], self.tobss_med_N[g])
    
    fig = plt.figure(figsize=(4,6.2))  # 5.7,5.1
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    colmap = _truncate_colormap(plt.get_cmap('Blues'),.3,1)
    '''img = ax1.scatter(self.Vmags_med[g], self_yarr[g], c=tobss_med[g],
    facecolors='none',
    cmap=plt.get_cmap('hot_r'), alpha=.4, s=s)#,
    #norm=colors.LogNorm(vmin=1, vmax=tobss_med[g].max()))'''
    # set colorbar
    if pltflag in [0,1]:
        img = ax1.scatter(self.Vmags_med[g], self.Ks_med[g],
                          c=self.sigmaRV_eff_med_H[g], cmap=plt.get_cmap(colmap),
                          s=0,
                          norm=colors.LogNorm(vmin=1, vmax=tobss_med.max()))
        # add transluscent points
        ax1.scatter(self.Vmags_med[g], self.Ks_med[g], c=self.tobss_med_H[g],
                    facecolors='none', cmap=plt.get_cmap(colmap), alpha=1, s=s,
                    norm=colors.LogNorm(vmin=1,vmax=tobss_med.max()))
        cbar_axes = fig.add_axes([.08, .08, .84, .032])
        cbar = fig.colorbar(img, cax=cbar_axes, orientation='horizontal')
        cbar.set_label('Total observing time per TOI [hours]', fontsize=12)
    
    # Get 50 best
    if pltflag in [0,2]:
        sort = np.argsort(self.tobss_med_H[g])[:50]
        x, y, c = self.Vmags_med[g], self.Ks_med[g], self.tobss_med_H[g]
        ax1.scatter(x[sort], y[sort], facecolor='none', edgecolor='k', s=s+2)
    
    # fill 'good' region
    if pltflag in [0,3]:
        m, b = get_best_fraction(self.Vmags_med[g], self.Ks_med[g],
                                 self.tobss_med_H[g])
        print m, b
        line = lambda x: m*x + b
        ax1.fill_between([3.5,11.35], [line(3.5),line(11.35)], 30, color='k',
                         alpha=.15)
        ax1.plot(np.repeat(11.35,2), [line(11.35),30], 'k--', lw=1.7)
        ax1.plot([3.5,11.35], [line(3.5),line(11.35)], 'k--', lw=1.7)
    
    ax1.set_xlabel('V', labelpad=1, fontsize=14, style='italic')
    ax1.set_ylabel('RV semi-amplitude [m s$^{-1}$]', fontsize=12)
    ax1.set_xlim((3.5,20)), ax1.set_ylim((0,12))
    ax1.set_yticks(np.arange(0,13,3)), ax1.set_yticklabels(np.arange(0,13,3))
    ax1.minorticks_on()

    # Plot Jmag
    if pltflag in [0,1]:
        # set colorbar
        ax2.scatter(self.Jmags_med[g], self.Ks_med[g],
                    c=self.sigmaRV_eff_med_N[g], cmap=plt.get_cmap(colmap), s=0, 
                    norm=colors.LogNorm(vmin=1, vmax=tobss_med.max()))
        # add transluscent points
        ax2.scatter(self.Jmags_med[g], self.Ks_med[g], c=self.tobss_med_N[g],
                    facecolors='none', cmap=plt.get_cmap(colmap), alpha=1, s=s,
                    norm=colors.LogNorm(vmin=1,vmax=tobss_med.max()))
    
    # Get 50 best
    if pltflag in [0,2]:
        sort = np.argsort(self.tobss_med_N[g])[:50]
        x, y = self.Jmags_med[g][sort], self.Ks_med[g][sort]
        ax2.scatter(x, y, facecolor='none', edgecolor='k', s=s+2)

    # fill 'good' region
    if pltflag in [0,3]:
        m, b = get_best_fraction(self.Jmags_med[g], self.Ks_med[g],
                                 self.tobss_med_N[g])
        print m, b
        line = lambda x: m*x + b
        ax2.fill_between([3.5,10.], [line(3.5),line(10.)], 30, color='k',
                         alpha=.15)
        ax2.plot(np.repeat(10,2), [line(10),30], 'k--', lw=1.7)
        ax2.plot([3.5,10], [line(3.5),line(10)], 'k--', lw=1.7)
    
    ax2.set_xlabel('J', labelpad=1, fontsize=14, style='italic')
    ax2.set_ylabel('RV semi-amplitude [m s$^{-1}$]', fontsize=12)
    ax2.set_xlim((3.5,15)), ax2.set_ylim((0,12))
    ax2.set_yticks(np.arange(0,13,3)), ax2.set_yticklabels(np.arange(0,13,3))
    ax2.minorticks_on()

    if pltflag not in [0,1]:
        ax1.axis('off')
        ax2.axis('off')
    fig.subplots_adjust(top=.98, bottom=.19, left=.15, right=.94, hspace=.21)
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


def get_best_fraction(mags, Ks, tobss):
    sort = np.argsort(tobss)[:50]
    mags, Ks, tobss = mags[sort], Ks[sort], tobss[sort]
    sort = Ks <= 12
    mags, Ks, tobss = mags[sort], Ks[sort], tobss[sort]
    m, b =  np.polyfit(mags, Ks, 1)
    line = lambda x: m*x + b
    g = Ks / line(mags) > 1
    while g.sum() < Ks.size-5:
        b -= 1e-2
        g = Ks > line(mags)
    return m, b


def MAD(xarr):
    return np.median(abs(xarr - np.median(xarr)))
