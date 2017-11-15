from imports import *


def compute_Nrv(sigeff, sigK):
    sigeff = np.ascontiguousarray(sigeff).astype(float)
    return 2 * (sigeff / sigK)**2


def plot_Nrv_estimates(pltt=True, label=False):
    # Get real planet data
    _,sigeffs, sigKs, Nrv_true, multflag = np.genfromtxt('Nrv_tests/FASTtests.dat').T

    # Compute Nrv
    Nrv_calc = compute_Nrv(sigeffs, sigKs)

    # plotting
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, aspect='equal')
    g = np.isfinite(multflag)#multflag == 1
    cax = ax.scatter(Nrv_true[g], Nrv_calc[g], edgecolors='none', c=np.log10(sigeffs[g]),
                     marker='o', s=100)
    cbar_axes = fig.add_axes([.1,.9,.87,.04])
    cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')#, pad=.11)
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

    
def plot_Nrv_F(self, errorbar=False, pltt=True, label=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    g = self.Teffs_med > 0

    if errorbar:
        ax.errorbar(self.Fs_med[g], self.Nrvs_med_I[g], self.Nrvs_emed_I[g],
                    fmt='r.', ms=12, capsize=0, elinewidth=.5)
        ax.errorbar(self.Fs_med[g], self.Nrvs_med_O[g], self.Nrvs_emed_O[g],
                    fmt='b.', ms=12, capsize=0, elinewidth=.5)
    else:
        ax.plot(self.Fs_med[g], self.Nrvs_med_I[g], 'ro', ms=10, alpha=.5)
        ax.plot(self.Fs_med[g], self.Nrvs_med_O[g], 'bo', ms=10, alpha=.5)
        
    ax.set_xscale('log')#, ax.set_yscale('log')
    ax.set_xlim((1e4,1e-1)), ax.set_ylim((1,2e2))
    ax.set_xlabel('Insolation ($S_{\oplus}$)'), ax.set_ylabel('$n_{RV}$')
    
    if pltt:
        plt.show()
    plt.close('all')


def plot_Nrv_mag(self, mag='V', errorbar=True, pltt=True, label=False):
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
        ax.errorbar(xarr[g], self.Nrvs_med_I[g], self.Nrvs_emed_I[g],
                    fmt='r.', ms=12, capsize=0, elinewidth=.5)
        ax.errorbar(xarr[g], self.Nrvs_med_O[g], self.Nrvs_emed_O[g],
                    fmt='b.', ms=12, capsize=0, elinewidth=.5)
    else:
        ax.plot(xarr[g], self.Nrvs_med_I[g], 'ro', ms=10, alpha=.5)
        ax.plot(xarr[g], self.Nrvs_med_O[g], 'bo', ms=10, alpha=.5)
        
    #ax.set_xlim((1e4,1e-1))
    ax.set_ylim((1,2e2))
    ax.set_xlabel('$%s$'%mag), ax.set_ylabel('$n_{RV}$')
    
    if pltt:
        plt.show()
    plt.close('all')
    

def plot_Nrvratio(self, pltt=True, label=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ratio = unumpy.uarray(self.Nrvs_med_I, self.Nrvs_emed_I) / \
            unumpy.uarray(self.Nrvs_med_O, self.Nrvs_emed_O)
    
    ax.errorbar(self.Teffs_med, unumpy.nominal_values(ratio),
                unumpy.std_devs(ratio),
                fmt='k.', ms=12, capsize=0, elinewidth=.5)
    xlim = ax.get_xlim()
    ax.fill_between(list(xlim), np.ones(2), color='k', alpha=.2)

    ax.set_yscale('log'), ax.set_ylim((1e-2,1e2)), ax.set_xlim(xlim)
    ax.set_xlabel('Teff (K)'), ax.set_ylabel('$n_{RV,I} / n_{RV,O}$')
    ax.minorticks_on()
    
    if pltt:
        plt.show()
    plt.close('all')
