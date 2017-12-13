import numpy as np
import pylab as plt
from uncertainties import unumpy as unp
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import matplotlib as mpl
import matplotlib.colors as colors

mpl.rc('xtick', labelsize=12)
mpl.rc('ytick', labelsize=12)
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


def plot_Nrv_estimates(pltt=True, label=False):
    # Get real planet data
    _,sigeffs, sigKs, Nrv_true, multflag = \
                                np.genfromtxt('Nrv_tests/FASTtests.dat').T

    # Compute Nrv
    Nrv_calc = compute_Nrv(sigeffs, sigKs)

    # plotting
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, aspect='equal')
    g = np.isfinite(multflag)#multflag == 1
    cax = ax.scatter(Nrv_true[g], Nrv_calc[g], edgecolors='none',
                     c=np.log10(sigeffs[g]), marker='o', s=100)
    cbar_axes = fig.add_axes([.1,.9,.87,.04])
    cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')
    cbar.set_label('log sigK')
    
    Nrvs = np.append(Nrv_true, Nrv_calc)
    arr = np.linspace(0, 3e2, 100)
    ax.plot(arr, arr, 'k--', lw=2)
    #ax.fill_between(arr, 3e2, arr, alpha=.3)
    ax.fill_between(arr, .9*arr, 1.1*arr, alpha=.5)
    ax.fill_between(arr, .8*arr, 1.2*arr, alpha=.5)
    ax.plot(arr, 1.1*arr, 'k--', lw=.9)
    ax.plot(arr, .9*arr, 'k--', lw=.9)

    ax.set_xscale('log'), ax.set_yscale('log')
    ax.set_xlabel('Actual n$_{RV}$'), ax.set_ylabel('Calculated n$_{RV}$')
    ax.set_xlim((1, 3e2)), ax.set_ylim((1, 3e2))
    ax.minorticks_on()

    plt.subplots_adjust(bottom=.12)
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


def plot_compare_tobs_H_N(self, pltt=True, label=False):

    h = self.spectrographs == 'H'
    n = self.spectrographs == 'N'
    ratio_NH = self.sigmaRV_acts[n] / self.sigmaRV_acts[h][:n.sum()]
    #ratio_NH = unp.uarray(self.tobss_med_N/self.tobss_med_H, self.tobss_emed_N/self.tobss_emed_H)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #ax.errorbar(self.Teffs_med, unp.nominal_values(ratio_NH), unp.std_devs(ratio_NH),
    #            fmt='ko', ms=0, elinewidth=.01)
    img = ax.scatter(self.Teffs[n], ratio_NH, s=2, alpha=.2, c=np.log10(self.vsinis[n]))
    #    cbar_axes = fig.add_axes([.9, .1, .05, .8])
    #    cbar = fig.colorbar(img, cax=cbar_axes)
    #ax.hist2d(self.Teffs[n], ratio_NH, bins=(50,np.logspace(-2,2,50)))

    ax.set_yscale('log')
    ax.set_ylim((1e-5, 1e5)), ax.set_xlim((25e2, 12e3))
    ax.set_ylabel('t$_{obs,nIR}$ / t$_{obs,opt}$'), ax.set_xlabel('Effective Temperature [K]')

    if pltt:
        plt.show()
    plt.close('all')
        
    
def plot_cumulative_detections_v_tobs(self, pltt=True, label=False,
                                      harps=True, nirps=True, tmax=1e3,
                                      observebesttargetsfirst=True):
    #tarr = np.arange(0,tmax,1)
    fig = plt.figure(figsize=(9.5,4.6))

    for i in range(rpedges.size-1):
	ax1 = fig.add_subplot(2, rpedges.size-1, i+1)

	# HARPS total observing time
	if harps:
            g = (self.rps_med >= rpedges[i]) & \
                (self.rps_med < rpedges[i+1])
            tobs = np.sort(self.tobss_med_H[g]) if observebesttargetsfirst \
                   else self.tobss_med_H[g]
	    tobs = np.append(0, np.cumsum(tobs))
	    Ndet = np.arange(tobs.size)
	    g = tobs <= tmax
	    ax1.plot(tobs[g], Ndet[g], 'b--', lw=2, drawstyle='steps')
            # plot derivative
	    if i == 0:
                ax1.plot([410,560], np.repeat(3.2,2), 'b--', lw=2)
		ax1.text(6e2, 3.2, 'Optical', color='b', fontsize=10, 
			 weight='normal', verticalalignment='center')
            # plot derivative
            ax2 = fig.add_subplot(2, rpedges.size-1, i+1+rpedges.size-1)
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, 'b--', lw=2)
            threshold = 1./20
            ax2.plot([0,tmax], np.repeat(threshold,2), 'k--', lw=.9)
            if i == 0:
                ax2.text(tmax*.25, threshold*1.03,
                         '%i hours / detection'%(1./threshold),
                         verticalalignment='bottom', fontsize=9)
                
        # NIRPS total observing time
	if nirps:
            g = (self.rps_med >= rpedges[i]) & (self.rps_med < rpedges[i+1])
            tobs = np.sort(self.tobss_med_N[g]) if observebesttargetsfirst \
                   else self.tobss_med_N[g]
            tobs = np.append(0, np.cumsum(tobs))
            Ndet = np.arange(tobs.size)
	    g = tobs <= tmax
            #_compute_curve_derivative(tobs[g], Ndet[g])
	    ax1.plot(tobs[g], Ndet[g], 'r-', drawstyle='steps')
	    if i == 0:
                ax1.plot([410,550], np.repeat(1.2,2), 'r-')
                ax1.text(6e2, 1.2, 'Near-IR', color='r', fontsize=10, 
                         weight='normal', verticalalignment='center')
            # plot derivative
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, 'r-')
                
        ax1.set_title('%s\n%s'%(rplabels[i], rplabels2[i]), fontsize=10)
        ax2.set_xlabel('Cumulative observing\ntime [hours]', fontsize=9)
        ax1.set_xticklabels('')
        ax1.minorticks_on(), ax2.minorticks_on()
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
        
    fig.subplots_adjust(left=.08, bottom=.125, top=.92, right=.929,
                        hspace=0, wspace=.22)
    if label:
        plt.savefig('plots/cumulativetobs_rp.png')
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
	    g = tobs <= tmax
	    ax1.plot(tobs[g], Ndet[g], 'b--', lw=2, drawstyle='steps')
            # plot derivative
	    if i == 3:
                ax1.plot([410,560], np.repeat(25,2), 'b--', lw=2)
		ax1.text(6e2, 25, 'Optical', color='b', fontsize=10, 
			 weight='normal', verticalalignment='center')
            # plot derivative
            ax2 = fig.add_subplot(2, Teffedges.size-1, -1*i+4+Teffedges.size-1)
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, 'b--', lw=2)
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
	    g = tobs <= tmax
            #_compute_curve_derivative(tobs[g], Ndet[g])
	    ax1.plot(tobs[g], Ndet[g], 'r-', drawstyle='steps')
	    if i == 3:
                ax1.plot([410,550], np.repeat(20,2), 'r-')
                ax1.text(6e2, 20, 'Near-IR', color='r', fontsize=10, 
                         weight='normal', verticalalignment='center')
            # plot derivative
            tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
            ax2.plot(tobs2, dNdt, 'r-')
            
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
            
    fig.subplots_adjust(left=.08, bottom=.125, top=.92, right=.929,
                        hspace=0, wspace=.22)
    if label:
        plt.savefig('plots/cumulativetobs_Teff.png')
    if pltt:
	plt.show()
    plt.close('all')
    
    
def plot_cumulative_detections_v_tobs_50(self, pltt=True, label=False,
                                         harps=True, nirps=True):
    fig = plt.figure(figsize=(5,5))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    # HARPS total observing time
    if harps:
        g = self.rps_med <= 2
        tobs = np.sort(self.tobss_med_H[g])
        tobs = np.append(0, np.cumsum(tobs))[:51]
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, 'b-', drawstyle='steps')
	ax1.text(.6, .3, 'Optical', color='b', fontsize=14,
		 weight='semibold', transform=ax1.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, 'b-')
        threshold = 1./20  # detections / hour
        ax2.plot(ax1.get_xlim(), np.repeat(threshold,2), 'k--')
        ax2.text(ax1.get_xlim()[1]*.08, threshold*.9,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='top', fontsize=10)
            
    # NIRPS total observing time
    if nirps:
        g = self.rps_med <= 2
        tobs = np.sort(self.tobss_med_N[g])
        tobs = np.append(0, np.cumsum(tobs))[:51]
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, 'r-', drawstyle='steps')
	ax1.text(.6, .25, 'nIR', color='r', fontsize=12, 
		 weight='semibold', transform=ax1.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, 'r-')    


    ax2.set_xlabel('Cumulative observing time [hours]', fontsize=11)
    ax2.set_ylabel('dN/dt\n[detections / hours]', fontsize=11, labelpad=0)
    ax1.set_ylabel('Total number of\nplanet detections', fontsize=11)	
    ax1.minorticks_on(), ax2.minorticks_on()
    ax1.set_xticklabels('')
    ax2.set_yscale('log')
    ax2.set_ylim((1e-2,1e0))
    ax2.set_xlim(ax2.set_xlim()), ax1.set_xlim(ax2.set_xlim())

    
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    ax3.set_yticks(np.logspace(2,0,3))
    ax3.set_yticklabels(['%i'%i for i in np.logspace(2,0,3)])
    ax3.set_ylabel('dt/dN\n[hours / detection]', fontsize=11,
                   labelpad=0)

    fig.subplots_adjust(left=.16, bottom=.09, top=.98, right=.86,
                        hspace=0)
    if label:
        plt.savefig('plots/cumulativetobs_50.png')
    if pltt:
	plt.show()
    plt.close('all')


def plot_cumulative_detections_v_tobs_50_random(self, pltt=True, label=False,
                                                harps=True, nirps=True,
                                                Nrand=10, seed=0):
    fig = plt.figure(figsize=(5.4,5))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    np.random.seed(int(seed))
    
    # HARPS total observing time
    if harps:
        g = self.rps_med < 4
        tobs = np.sort(self.tobss_med_H[g])
        tobs = np.append(0, np.cumsum(tobs))[:51]
        Ndet = np.arange(tobs.size)
        ax1.plot(tobs, Ndet, 'b--', lw=2, drawstyle='steps')
	ax1.text(.1, .8, 'Optical', color='b', fontsize=14,
		 weight='normal', transform=ax1.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, 'b--', lw=2)
        threshold = 1./20  # detections / hour
        ax2.plot([1,1e5], np.repeat(threshold,2), 'k--')
        ax2.text(2, threshold*1.05,
                 '%i hours / detection'%(1./threshold),
                 verticalalignment='bottom', fontsize=10)
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = self.tobss_med_H[g][inds]
            tobs = np.append(0, np.cumsum(tobs))[:51]
            ax1.plot(tobs, Ndet, 'b--', drawstyle='steps', lw=.3)
            tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            ax2.plot(tobs2, dNdt, 'b--', lw=.3)
            
    # NIRPS total observing time
    if nirps:
        g = self.rps_med < 4
        tobs = np.sort(self.tobss_med_N[g])
        tobs = np.append(0, np.cumsum(tobs))[:51]
	Ndet = np.arange(tobs.size)
	ax1.plot(tobs, Ndet, 'r-', drawstyle='steps')
	ax1.text(.1, .73, 'Near-IR', color='r', fontsize=12, 
		 weight='normal', transform=ax1.transAxes)
        tobs2, dNdt = _compute_curve_derivative(tobs, Ndet)
        ax2.plot(tobs2, dNdt, 'r-')
        inds = np.arange(g.sum())
        for i in range(Nrand):
            np.random.shuffle(inds)
            tobs = self.tobss_med_N[g][inds]
            tobs = np.append(0, np.cumsum(tobs))[:51]
            ax1.plot(tobs, Ndet, 'r-', drawstyle='steps', lw=.3)
            tobs2, dNdt = _compute_curve_derivative(tobs, Ndet, 50)
            ax2.plot(tobs2, dNdt, 'r-', lw=.3)

    ax1.set_xscale('log'), ax1.set_xlim((1,1e5)), ax1.set_ylim((-1,50))
    ax1.set_ylabel('Total number of\nplanet detections', fontsize=11)	
    ax1.minorticks_on()
    ax1.set_xticklabels('')
    ax2.set_xlabel('Cumulative observing time [hours]', fontsize=11)
    ax2.set_ylabel('dN/dt\n[detections / hour]', fontsize=11, labelpad=0)
    ax1.set_title('r$_p < 4$ R$_{\oplus}$')
    ax2.set_yscale('log'), ax2.set_ylim((1e-4,2))
    ax2.set_xscale('log'), ax2.set_xlim((1,1e5))
    #ax2.set_yticks(np.logspace(-4,0,5))
    #ax2.set_yticklabels(['0.01','0.1','1'])

    # second derivative axis
    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    #ax3.set_yticks(np.logspace(4,0,3))
    #ax3.set_yticklabels(['%i'%i for i in np.logspace(2,0,3)])
    ax3.set_ylabel('dt/dN\n[hours / detection]', fontsize=10,
                   labelpad=0)

    fig.subplots_adjust(left=.14, bottom=.14, top=.95, right=.87, hspace=0)
    if label:
        plt.savefig('plots/cumulativetobs_50random.png')
    if pltt:
	plt.show()
    plt.close('all')

    
def plot_cumulative_detections_v_tobs_MR(self, pltt=True, label=False,
                                         harps=True, nirps=True, tmax=1e3,
                                         observebesttargetsfirst=True):
    fig = plt.figure(figsize=(7,5.5))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    # HARPS total observing time
    if harps:
        g = self.rps_med <= 4
        Nrvs = self.Nrvs_med_H[g] * (.327 / .188)**2
        tobss = Nrvs * (self.texps_med_H[g] + toverhead) / 60
        tobs = np.sort(tobss) if observebesttargetsfirst else tobss
        tobs = np.append(0, np.cumsum(tobs))
        Ndet = np.arange(tobs.size)
	g = tobs <= tmax
	#_compute_curve_derivative(tobs[g], Ndet[g])
        ax1.plot(tobs[g], Ndet[g], 'b-', drawstyle='steps')
	ax1.text(.6, .3, 'Optical', color='b', fontsize=14,
		 weight='normal', transform=ax1.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
        ax2.plot(tobs2, dNdt, 'b-')
        threshold = 1./20  # detections / hour
        ax2.plot([0,tmax], np.repeat(threshold,2), 'k--')
        ax2.text(tmax*.08, threshold*.8, '%i hours / detection'%(1./threshold),
                 verticalalignment='top', fontsize=10)
 
    # NIRPS total observing time
    if nirps:
        g = self.rps_med <= 4
        Nrvs = self.Nrvs_med_N[g] * (.327 / .188)**2
        tobss = Nrvs * (self.texps_med_N[g] + toverhead) / 60
        tobs = np.sort(tobss) if observebesttargetsfirst else tobss
        tobs = np.append(0, np.cumsum(tobs))
	Ndet = np.arange(tobs.size)
	g = tobs <= tmax
	#_compute_curve_derivative(tobs[g], Ndet[g])
        ax1.plot(tobs[g], Ndet[g], 'r-', drawstyle='steps')
	ax1.text(.6, .25, 'nIR', color='r', fontsize=12, 
		 weight='normal', transform=ax1.transAxes)
        # plot derivative
        tobs2, dNdt = _compute_curve_derivative(tobs[g], Ndet[g])
        ax2.plot(tobs2, dNdt, 'r-')

    ax2.set_xlabel('Cumulative observing time [hours]', fontsize=11)
    ax2.set_ylabel('dN/dt\n[detections / hours]', fontsize=11)
    ax1.set_ylabel('Total number of\nplanet detections', fontsize=11)	
    ax1.minorticks_on(), ax2.minorticks_on()
    ax1.set_xticklabels('')
    ax2.set_yscale('log')
    ax2.set_ylim((1e-2,1e0))
    ax1.set_xlim((0,tmax)), ax2.set_xlim((0,tmax))

    ax3 = ax2.twinx()
    ax3.set_ylim(tuple(1./np.ascontiguousarray(ax2.get_ylim())))
    ax3.set_yscale('log')
    ax3.set_yticks(np.logspace(2,0,3))
    ax3.set_yticklabels(['%i'%i for i in np.logspace(2,0,3)],
                        fontsize=11)
    ax3.set_ylabel('dt/dN\n[hours / detection]', fontsize=11)

    #fig.subplots_adjust(left=.08, bottom=.17, top=.86, right=.97, wspace=.24)
    if label:
        plt.savefig('plots/cumulativetobs_MR.png')
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


def plot_identifying_best_50(self, self_xarr, self_yarr,
                             s=3, pltt=True, label=False):
    '''
    Plot variables to try and guess which parameters best tell which small 
    planets we can observe most efficiently.
    e.g. self_xarr = self.Vmags_med
    '''
    g = (self.rps_med < 4)# & (self.spectrographs == 'N')
    #tobss_med = np.mean([self.tobss_med_H, self.tobss_med_N], axis=0)
    tobss_med = self.tobss_med_N
    assert self_xarr.size == tobss_med.size
    assert self_yarr.size == tobss_med.size
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colmap = _truncate_colormap(plt.get_cmap('hot_r'),.1,1)
    img = ax.scatter(self_xarr[g], self_yarr[g], c=np.log10(tobss_med[g]),
                     cmap=plt.get_cmap(colmap), alpha=.5, s=s)
    cbar_axes = fig.add_axes([.9, .1, .05, .8])
    cbar = fig.colorbar(img, cax=cbar_axes)
    
    # Get 50*100 best
    sort = np.argsort(tobss_med[g])[:50]
    x, y = self_xarr[g][sort], self_yarr[g][sort]
    ax.scatter(x, y, facecolor='none', edgecolor='k', s=s)

    #ax.set_yscale('log')

    if pltt:
        plt.show()
    plt.close('all')
