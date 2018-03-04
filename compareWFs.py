from imports import *
from compute_sigK_QPGP import compute_sigmaK_GP


def set_thetas():
    a, l, G, Pgp, s = 2.8,  59.1, 1.2, 38.6, 0.25
    P, K = 32.93963, 3.18
    sigRV = 3.63
    return (a, l, G, Pgp, s), (P, K), sigRV


def get_true_WF():
    bjd = np.loadtxt('/Users/ryancloutier/Research/K2_18/input_data/2M1130+0735_05052017_w7s5c5_NAIRA_v1.rdb', skiprows=2)[:-4,0]
    return bjd


def get_uniform_WF():
    bjd = np.loadtxt('/Users/ryancloutier/Research/K2_18/input_data/2M1130+0735_05052017_w7s5c5_NAIRA_v1.rdb', skiprows=2)[:-4,0]
    return np.linspace(0, bjd.max()-bjd.min(), bjd.size)


def compare_WFs(N):
    '''Compute a set of sigKs using a GP with the K2-18 HARPS WF and a uniform 
    WF to test the sensitivity of sigK on the assumed WF.'''
    GPtheta, keptheta, sig_phot = set_thetas()
    a, l, G, Pgp, s = GPtheta
    P, K = keptheta
    sigK_target = K / 3.
    
    N = int(N)
    sigKs_uni, sigKs_harps = np.zeros(N), np.zeros(N)
    for i in range(N):
        print float(i) / N
        gp = george.GP(a*(george.kernels.ExpSquaredKernel(l) + \
                          george.kernels.ExpSine2Kernel(G,Pgp)))
        t1 = get_uniform_WF()
        t2 = get_true_WF()
        assert t1.size == t2.size
        erv = np.repeat(sig_phot, t1.size)
        gp.compute(t1, np.sqrt(erv**2 + s**2))
        rv_act1 = gp.sample(t1)
        gp.compute(t2, np.sqrt(erv**2 + s**2))
        rv_act2 = gp.sample(t2)

        # get planet model
        rv_kep1 = -K*np.sin(2*np.pi*foldAt(t1, P))
        rv_kep2 = -K*np.sin(2*np.pi*foldAt(t2, P))

        # get total rv signal with noise
        rv1 = rv_act1 + rv_kep1 + np.random.randn(t1.size) * sig_phot
        rv2 = rv_act2 + rv_kep2 + np.random.randn(t1.size) * sig_phot

        # compute sigK
        theta = P, 0, K, a, l, G, Pgp, s
        sigKs_uni[i] = compute_sigmaK_GP(theta, t1, rv1, erv)
        sigKs_harps[i] = compute_sigmaK_GP(theta, t2, rv2, erv)

    return sigKs_uni, sigKs_harps


def plot_compareWF(sigKs_uni, sigKs_harps):
    plt.scatter(sigKs_uni, sigKs_harps)
    sigKs = np.append(sigKs_uni, sigKs_harps)
    g = np.isfinite(sigKs)
    plt.plot([sigKs[g].min(),sigKs[g].max()], [sigKs[g].min(),sigKs[g].max()],
             'k--')
    plt.xlabel('$\sigma_{K,uniform}$')
    plt.ylabel('$\sigma_{K,HARPS}$')
    g = np.isfinite(sigKs_harps) & np.isfinite(sigKs_uni)
    plt.title('residual RMS = %.3f'%(np.std(sigKs_uni[g]-sigKs_harps[g])))
    plt.show()
