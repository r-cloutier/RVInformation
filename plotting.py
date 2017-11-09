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
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')
    g = multflag == 1
    cax = ax.scatter(Nrv_true[g], Nrv_calc[g], edgecolors='none', c=np.log10(sigKs[g]), 
		     marker='o', s=100)
    g = multflag == 0
    cax = ax.scatter(Nrv_true[g], Nrv_calc[g], edgecolors='none', c=np.log10(sigKs[g]),
                     marker='s', s=100)
    #cbar_axes = fig.add_axes([.1,.9,.87,.04])
    #cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')#, pad=.11)
    #cbar.set_label('log sigK')
    
    Nrvs = np.append(Nrv_true, Nrv_calc)
    arr = np.linspace(0, 3e2, 100)
    ax.plot(arr, arr, 'k--', lw=2)
    ax.fill_between(arr, 3e2, arr, alpha=.3)

    ax.set_xscale('log'), ax.set_yscale('log')
    ax.set_xlabel('True n$_{RV}$'), ax.set_ylabel('Calculated n$_{RV}$')
    ax.set_xlim((0, 3e2)), ax.set_ylim((0, 3e2))
    ax.minorticks_on()
    
    if pltt:
        plt.show()
    plt.close('all')
