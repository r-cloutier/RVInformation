'''
The value of sigK from the Fisher information when using a QP GP activity model 
depends on the parameter values {Krv,a,l,G,P,s}. From the RV analysis of real 
systems these parameter values are not known exactly and have an associated 
uncertainty which leads to an uncertainty in sigK and hence on nRV. Here we 
sample parameter values from their marginalized posteriors to compute the 
distribution in Fisher sigK.
'''
from imports import *
from compute_sigK_QPGP import compute_sigmaK_GP
from rvmodel import get_rv1

global fnames
fnames = ['K218','LHS1140','Kep78HARPSN','Kep21HARPSN','CoRoT7']


def sigK_K218(N=1e3, P=32.93963, T0=2457264.39157):
    '''From Cloutier et al 2017'''
    # get time-series
    bjd, rv, erv = np.loadtxt('Nrv_tests/K218.dat').T
    kepc = get_rv1((8.962, 2457264.55, 0, 4.63, 0, 0), bjd)
    rv -= kepc

    # sample parameter posteriors approximating the PDFs as Gaussian
    N = int(N)
    Krvs = _random_normal_draws(3.18, .71, N, positive=True)
    As = _random_normal_draws(.67, 1.8, N, positive=True)
    ls = _random_normal_draws(64.49, 15, N, positive=True)
    Gs = _random_normal_draws(1.2, .45, N, positive=True)
    Ps = _random_normal_draws(38.6, .9, N, positive=True)
    ss = _random_normal_draws(.25, .4, N, positive=True)

    # Monte-Carlo computation of sigKs
    sigKs = np.zeros(N)
    for i in range(N):
        theta = P, T0, Krvs[i], As[i], ls[i], Gs[i], Ps[i], ss[i]
        sigKs[i] = compute_sigmaK_GP(theta, bjd, rv, erv)

    # save output
    _save_draws('K218', Krvs, As, ls, Gs, Ps, ss, sigKs)



def sigK_LHS1140(N=1e3, P=24.73712, T0=2456915.6997):
    '''From Dittmann et al 2017'''
    # get time-series
    bjd, rv, erv = np.loadtxt('Nrv_tests/LHS1140.dat').T
    
    # sample parameter posteriors approximating the PDFs as Gaussian
    N = int(N)
    Krvs = _random_normal_draws(5.3, 1.1, N, positive=True)
    As = _random_normal_draws(9, 5, N, positive=True)
    ls = _random_normal_draws(393/np.sqrt(2), 30/np.sqrt(2), N, positive=True)
    Gs = _random_normal_draws(2, .2, N, positive=True)
    Ps = _random_normal_draws(134, 10, N, positive=True) # median of posterior
    ss = _random_normal_draws(3, 1, N, positive=True)

    # Monte-Carlo computation of sigKs
    sigKs = np.zeros(N)
    for i in range(N):
        theta = P, T0, Krvs[i], As[i], ls[i], Gs[i], Ps[i], ss[i]
        sigKs[i] = compute_sigmaK_GP(theta, bjd, rv, erv)

    # save output
    _save_draws('LHS1140', Krvs, As, ls, Gs, Ps, ss, sigKs)



def sigK_Kep78(N=1e3, P=.35500744, T0=2454953.95995):
    '''From Grunblatt et al 2015'''
    # get time-series
    bjd, rv, erv = np.loadtxt('Nrv_tests/Kepler78_HARPSN.dat').T
    
    # sample parameter posteriors approximating the PDFs as Gaussian
    N = int(N)
    Krvs = _random_normal_draws(1.86, .25, N, positive=True)
    As = _random_normal_draws(5.6, 1.7, N, positive=True)
    ls = _random_normal_draws(26.1/np.sqrt(2), 15/np.sqrt(2), N, positive=True)
    Gs = _random_normal_draws(1/np.sqrt(2*.28**2), 1/np.sqrt(2*.28**2)/5.6,
                              N, positive=True)
    Ps = _random_normal_draws(13.26, .12, N, positive=True)
    ss = _random_normal_draws(1.1, .5, N, positive=True)

    # Monte-Carlo computation of sigKs
    sigKs = np.zeros(N)
    for i in range(N):
        theta = P, T0, Krvs[i], As[i], ls[i], Gs[i], Ps[i], ss[i]
        sigKs[i] = compute_sigmaK_GP(theta, bjd, rv, erv)
        sigKs[i] *= np.sqrt(109./193)  # correct for lack of HIRES RVs to compare to measured sigK = 0.25 m/s

    # save output
    _save_draws('Kep78HARPSN', Krvs, As, ls, Gs, Ps, ss, sigKs)

    

def sigK_Kep21(N=1e3, P=2.78578, T0=2456798.7188):
    '''From Lopez-Morales et al 2016'''
    # get time-series
    bjd, rv, erv = np.loadtxt('Nrv_tests/Kepler21_HARPSN.dat').T
    
    # sample parameter posteriors approximating the PDFs as Gaussian
    N = int(N)
    Krvs = _random_normal_draws(2.12, .66, N, positive=True)
    As = _random_normal_draws(6.7, 1.4, N, positive=True)
    ls = _random_normal_draws(24.04/np.sqrt(2), .09/np.sqrt(2),
                              N, positive=True)
    Gs = _random_normal_draws(1/.42, .12/.42, N, positive=True)
    Ps = _random_normal_draws(12.6, .02, N, positive=True)
    rmsinst, rmsgran = unp.uarray(.9, .1), unp.uarray(1.76, .04)
    s = unp.sqrt(rmsinst**2 + rmsgran**2)
    ss = _random_normal_draws(float(unp.nominal_values(s)), 
			      float(unp.std_devs(s)), N, 
			      positive=True)

    # Monte-Carlo computation of sigKs
    sigKs = np.zeros(N)
    for i in range(N):
        theta = P, T0, Krvs[i], As[i], ls[i], Gs[i], Ps[i], ss[i]
        sigKs[i] = compute_sigmaK_GP(theta, bjd, rv, erv)

    # save output
    _save_draws('Kep21HARPSN', Krvs, As, ls, Gs, Ps, ss, sigKs)


    
def sigK_CoRoT7(N=1e3, P=.85359165, T0=2454398.0769):
    '''From Haywood et al 2015'''
    # get time-series
    bjd, rv, erv = np.loadtxt('Nrv_tests/CoRoT7.dat').T
    kepc = get_rv1((3.7, 2455953.54, 0, 6.01, 0, 0), bjd)
    rv -= kepc

    # sample parameter posteriors approximating the PDFs as Gaussian
    N = int(N)
    Krvs = _random_normal_draws(3.42, .66, N, positive=True)
    As = _random_normal_draws(7, 2, N, positive=True)
    ls = _random_normal_draws(20.6, 2.5, N, positive=True)
    Gs = _random_normal_draws(1, .1, N, positive=True)
    Ps = _random_normal_draws(23.81, .03, N, positive=True)
    rmsRVrot, rmsRVconv, rmsRVtot = .46, 1.82, 3.92  # FF' contributions
    ss = _random_normal_draws(np.sqrt(rmsRVtot**2 - rmsRVrot**2 - rmsRVconv**2), 
			      0, N, positive=True)

    # Monte-Carlo computation of sigKs
    sigKs = np.zeros(N)
    for i in range(N):
        theta = P, T0, Krvs[i], As[i], ls[i], Gs[i], Ps[i], ss[i]
        sigKs[i] = compute_sigmaK_GP(theta, bjd, rv, erv)

    # save output
    _save_draws('CoRoT7', Krvs, As, ls, Gs, Ps, ss, sigKs)

    
def _random_normal_draws(mu, sig, N, positive=False):
    draws = mu + np.random.randn(N) * sig
    if positive:
        draws2 = draws + 0
        while np.any(draws < 0) or draws.size < N:
            draws2 = np.append(draws2, mu + np.random.randn(N) * sig)
            draws = draws2[draws2 >= 0]
        draws = draws[:N]
    return draws


def _save_draws(fname, Krvs, As, ls, Gs, Ps, ss, sigKs):
    outp = np.array([Krvs, As, ls, Gs, Ps, ss, sigKs]).T
    np.savetxt('Nrv_tests/GPtest_%s.dat'%fname, outp, delimiter='\t',
               fmt='%.5f')


def plot_correlations(fname, sigKtarget=0, pltt=True, label=False):
    d = np.loadtxt('Nrv_tests/GPtest_%s.dat'%fname)
    labels = ['$K_{RV}$','$a$','$\lambda$','$\Gamma$','$P_{GP}$','$\sigma$']
    fig = plt.figure(figsize=(14,3))
    for i in range(6):
        ax = fig.add_subplot(1,6,i+1)
        ax.scatter(d[:,i], d[:,-1], s=10)
        ax.set_xlabel(labels[i])
        ax.set_yscale('log')
        if i == 0:
            ax.set_ylabel('$\sigma_K$ (m/s)')
        else:
            ax.set_yticklabels('')
        if sigKtarget != 0:
            ax.axhline(sigKtarget, lw=.7)    
    fig.subplots_adjust(bottom=.19, left=.07, right=.99, top=.96, wspace=.05)
    if label:
        plt.savefig('plots/GPtestcorrelations_%s.png'%fname)
    if pltt:
        plt.show()
    plt.close('all')


def plot_sigK_hist(fname, sigKtarget=0, pltt=True, label=False):
    sigKs = np.loadtxt('Nrv_tests/GPtest_%s.dat'%fname)[:,-1]
    sigKs = sigKs[np.isfinite(sigKs)]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(sigKs, bins=np.logspace(np.log10(sigKs.min()),
                                    np.log10(sigKs.max()),30), alpha=.9)
    ax.set_xscale('log')
    if sigKtarget != 0:
        ax.axvline(sigKtarget, lw=1.1,
                   label='Measured $\sigma_K$\n= %.2f m/s'%sigKtarget)
    v = np.percentile(sigKs, (16,84))
    ylim = ax.get_ylim()
    ax.fill_between(v, ylim[0], ylim[1], alpha=.7, color='k')
    ax.set_xlabel('$\sigma_K$ (m/s)')
    ax.set_ylim(ylim)
    ax.legend(bbox_to_anchor=(.35, .8))
    fig.subplots_adjust(bottom=.13)
    if label:
        plt.savefig('plots/GPtestsigKhist_%s.png'%fname)
    if pltt:
        plt.show()
    plt.close('all')
