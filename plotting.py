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

    
def plot_Nrv_F(self, pltt=True, label=False):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    g = self.spectrographs == 'I'
    xarr = np.median(self.Fs[g])
    
    ax.scatter(self.Fs[g], self.Nrvs[g])
    ax.set_xlim((1e3, 1e-1)), ax.set_ylim((0,200))
    ax.set_xscale('log')

    if pltt:
        plt.show()
    plt.close('all')


def median4star(starnums, attribute):
    '''
    Get the median attribute across all iterations for a single star.
    '''
    assert starnums.shape == attribute.shape
    unique_starnums = np.unique(starnums)
    nstars = unique_starnums.size
    attribute_out = np.zeros(nstars)
    for i in range(nstars):
        g = starnums == unique_starnums[i]
        attribute_out[i] = np.median(attribute[g])
    return unique_starnums, attribute_out
