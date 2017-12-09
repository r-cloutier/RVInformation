import numpy as np
import pylab as plt
from uncertainties import unumpy as unp

global toverhead, Teffedges, Tefflabels, rpedges, rplabels, rplabels2
toverhead = 5.
Teffedges = np.array([25e2, 32e2, 38e2, 76e2, 12e3])
Tefflabels = ['mid-late M', 'early-mid M', 'FGK', 'BA']
rpedges = np.array([0, 1.25, 2, 4, 30])
rplabels = ['Earths','Super-Earths','Sub-Neptunes','Giants']
rplabels2 = ['$<1.25 R_{\oplus}$','$1.25-2 R_{\oplus}$',
                 '$2-4 R_{\oplus}$','$>4 R_{\oplus}$']

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


def plot_cumulative_detections_v_tobs(self, pltt=True, label=False,
                                      harps=True, nirps=True, tmax=3e2):
    #tarr = np.arange(0,tmax,1)
    fig = plt.figure(figsize=(12,3.3))

    for i in range(Teffedges.size-1):
	ax = fig.add_subplot(1, Teffedges.size-1, i+1)

	# HARPS total observing time
	if harps:
            g = (self.rps_med >= rpedges[i]) & \
                (self.rps_med < rpedges[i+1])
	    tobs = np.append(0, np.cumsum(np.sort(self.tobss_med_H[g])))
	    Ndet = np.arange(tobs.size)
	    g = tobs <= tmax
	    ax.plot(tobs[g], Ndet[g], 'b-', drawstyle='steps')
	    if i == 0:
		ax.text(.5, .3, 'HARPS-like', color='b', fontsize=10, 
			weight='normal', transform=ax.transAxes)
 
        # NIRPS total observing time
	if nirps:
            g = (self.rps_med >= rpedges[i]) & (self.rps_med < rpedges[i+1])
            tobs = np.append(0, np.cumsum(np.sort(self.tobss_med_N[g])))
            Ndet = np.arange(tobs.size)
	    g = tobs <= tmax
            ax.plot(tobs[g], Ndet[g], 'r-', drawstyle='steps')
	    if i == 0:
                ax.text(.5, .2, 'NIRPS-like', color='r', fontsize=11, 
                        weight='normal', transform=ax.transAxes)

        ax.set_title('%s\n%s'%(rplabels[i], rplabels2[i]), fontsize=10)
        ax.set_xlabel('Total observing time [hours]', fontsize=9)
	if i == 0:
	    ax.set_ylabel('Cumulative number of\nplanet detections')	
	ax.minorticks_on()
    	#ax.set_xscale('log')
    	ax.set_xlim((0,tmax))

    fig.subplots_adjust(left=.08, bottom=.17, top=.86, right=.97, wspace=.24)
    if label:
        plt.savefig('plots/cumulativetobs_rp.png')
    if pltt:
	plt.show()
    plt.close('all')

    
def plot_cumulative_detections_v_tobs_MR(self, pltt=True, label=False,
                                         harps=True, nirps=True, tmax=1e3):
    fig = plt.figure(figsize=(7,5.5))
    ax = fig.add_subplot(111)
    
    # HARPS total observing time
    if harps:
        g = self.rps_med <= 2
        Nrvs = self.Nrvs_med_H[g] * (.327 / .188)**2
        tobss = Nrvs * (self.texps_med_H[g] + toverhead) / 60
	tobs = np.append(0, np.cumsum(np.sort(tobss)))
	Ndet = np.arange(tobs.size)
	g = tobs <= tmax
	ax.plot(tobs[g], Ndet[g], 'b-', drawstyle='steps')
	ax.text(.6, .3, 'Optical', color='b', fontsize=14,
		weight='normal', transform=ax.transAxes)
 
    # NIRPS total observing time
    if nirps:
        g = self.rps_med <= 2
        Nrvs = self.Nrvs_med_N[g] * (.327 / .188)**2
        tobss = Nrvs * (self.texps_med_N[g] + toverhead) / 60
        tobs = np.append(0, np.cumsum(np.sort(tobss)))
	Ndet = np.arange(tobs.size)
	g = tobs <= tmax
	ax.plot(tobs[g], Ndet[g], 'r-', drawstyle='steps')
	ax.text(.6, .25, 'nIR', color='r', fontsize=12, 
		weight='normal', transform=ax.transAxes)

    ax.set_xlabel('Total observing time [hours]')
    ax.set_ylabel('Cumulative number of\nplanet detections')	
    ax.minorticks_on()
    ax.set_xlim((0,tmax))

    #fig.subplots_adjust(left=.08, bottom=.17, top=.86, right=.97, wspace=.24)
    if label:
        plt.savefig('plots/cumulativetobs_MR.png')
    if pltt:
	plt.show()
    plt.close('all')

    #return tobs, Ndet

