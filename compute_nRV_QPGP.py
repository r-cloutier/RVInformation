from imports import *
import sympy
from rvmodel import get_rv1

sympy.init_printing()


def get_timeseries():
    tarr = np.random.rand(100) * 100
    P, T0, K, sig = 10., 0., 3., 1.
    rvarr = get_rv1((P,T0,0,K,0,0), tarr) + np.random.randn(tarr.size)*sig
    ervarr = np.repeat(sig, tarr.size)
    return tarr, rvarr, ervarr


def compute_sigmaK_GP(theta, t_arr, rv_arr, erv_arr):
    '''P, T0, Krv, a, l, G, Pgp = theta'''
    B = _compute_Fisher_information_GP(theta, t_arr, rv_arr, erv_arr)
    return np.sqrt(np.diag(np.linalg.inv(B))[0])


def _covariance_matrix(theta, t, sig=np.zeros(0)):
    a, l, G, P = theta
    t = t.reshape(t.size, 1)
    K = np.zeros ((t.size, t.size))
    K += -.5/l**2 * (np.tile(t[:,0], t.shape) - \
                     np.tile(t[:,0], t.shape).T)**2
    K += -G**2 * np.sin(np.pi/P * abs(np.tile(t[:,0], t.shape) - \
                                      np.tile(t[:,0], t.shape).T))**2
    K = a*a*np.exp(K)
    if sig.size == t.size:
        return K + sig**2 * np.eye(t.size)
    else:
        return K


def _compute_Fisher_information_GP(theta, t_arr, rv_arr, erv_arr):
    '''Compute the Fisher information matrix with very complex terms for a 
    circular keplerian RV model and a QP GP red noise model.
    theta = list of parameter values (P, T0, Krv, a, l, G, Pgp)
    '''
    # get orbital phase array :)
    assert len(theta) == 7
    P, T0, Krv, a, l, G, Pgp = theta
    sort = np.argsort(t_arr)
    t_arr, rv_arr, erv_arr = t_arr[sort], rv_arr[sort], erv_arr[sort]
    phase_arr = foldAt(t_arr, P, T0)
    thetavals = theta[2:]
    
    # define variables
    Krv, a, l, G, P = sympy.symbols('Krv a lambda Gamma P')
    symbol_vals = Krv, a, l, G, P
    assert len(symbol_vals) == len(thetavals)
    
    # define time-series symbols
    t, phi, rv, sig = sympy.symbols('t phi RV sigma')
    y = rv - (-Krv*sympy.sin(phi))
    
    # compute QP covariance function
    k = a*a*sympy.exp(-.5*t**2 / l**2 - G*G*sympy.sin(np.pi*abs(t)/P)**2)
    symbol_arrs = t, phi, rv, y, k

    # get arrays
    Kinv = np.linalg.inv(_covariance_matrix(theta[3:], t_arr, erv_arr))
    thetaarrs = t_arr, phase_arr, rv_arr, Kinv

    # compute Fisher matrix
    Nparams = len(thetavals)
    B = np.zeros((Nparams, Nparams))
    for i in range(Nparams):
        for j in range(Nparams):
            if i > j:
                B[i,j] = B[j,i]
            else:
                B[i,j] = _compute_Fisher_entry(symbol_vals[i], symbol_vals[j],
                                               symbol_vals, symbol_arrs,
                                               thetavals, thetaarrs)
    return B


def _compute_Fisher_entry(symbol_i, symbol_j, symbol_values, symbol_arrays,
                          thetavals, thetaarrs):
    '''
    Compute the Fisher information entry for one pair of model 
    parameters. I.e. the partial of the lnlikelihood wrt to each 
    symbol.
    '''
    # compute partial expressions
    Krv_sym, a_sym, l_sym, G_sym, P_sym = symbol_values
    t_sym, phi_sym, rv_sym, y_sym, K_sym = symbol_arrays
    dy_didj = sympy.lambdify([Krv_sym, phi_sym, rv_sym],
                             sympy.diff(y_sym, symbol_i, symbol_j), 'numpy')
    dy_di = sympy.lambdify([Krv_sym, phi_sym, rv_sym],
                           sympy.diff(y_sym, symbol_i), 'numpy')
    dy_dj = sympy.lambdify([Krv_sym, phi_sym, rv_sym],
                           sympy.diff(y_sym, symbol_j), 'numpy')
    dK_didj = sympy.lambdify([a_sym, l_sym, G_sym, P_sym, t_sym],
                             sympy.diff(K_sym, symbol_i, symbol_j), 'numpy')
    dK_di = sympy.lambdify([a_sym, l_sym, G_sym, P_sym, t_sym],
                           sympy.diff(K_sym, symbol_i), 'numpy')
    dK_dj = sympy.lambdify([a_sym, l_sym, G_sym, P_sym, t_sym],
                           sympy.diff(K_sym, symbol_j), 'numpy')
    
    # evaluate partials at input values
    K_val, a_val, l_val, G_val, P_val = thetavals
    t_arr, phase_arr, rv_arr, Kinv =  thetaarrs
    N = t_arr.size
    deltat_arr = np.tile(t_arr, (N,1)) - np.tile(t_arr, (N,1)).T    
    y_arr = _intovector(rv_arr - _keplerian(K_val, phase_arr), N)
    dy_didj = _intovector(dy_didj(K_val, phase_arr, rv_arr), N)
    dy_di = _intovector(dy_di(K_val, phase_arr, rv_arr), N)
    dy_dj = _intovector(dy_dj(K_val, phase_arr, rv_arr), N)
    dK_didj = _intomatrix(dK_didj(a_val, l_val, G_val, P_val, deltat_arr), N)
    dK_di = _intomatrix(dK_di(a_val, l_val, G_val, P_val, deltat_arr), N)
    dK_dj = _intomatrix(dK_dj(a_val, l_val, G_val, P_val, deltat_arr), N)

    # get Fisher terms to sum
    terms = np.zeros(11)
    terms[0] = np.dot(dy_didj.T, np.dot(Kinv, y_arr))
    terms[1] = -np.dot(dy_di.T, np.dot(Kinv, np.dot(dK_dj, np.dot(Kinv, y_arr))))
    terms[2] = np.dot(dy_di.T, np.dot(Kinv, dy_dj))
    terms[3] = -np.dot(dy_dj.T,np.dot(Kinv, np.dot(dK_di, np.dot(Kinv, y_arr))))
    terms[4] = np.dot(y_arr.T, np.dot(Kinv,
                                     np.dot(dK_dj,
                                            np.dot(Kinv,
                                                   np.dot(dK_di,
                                                          np.dot(Kinv, y_arr))))))
    terms[5] = -np.dot(y_arr.T, np.dot(Kinv, np.dot(dK_didj, np.dot(Kinv, y_arr))))
    terms[6] = np.dot(y_arr.T, np.dot(Kinv,
                                     np.dot(dK_di,
                                            np.dot(Kinv,
                                                   np.dot(dK_dj,
                                                          np.dot(Kinv, y_arr))))))
    terms[7] = -np.dot(y_arr.T, np.dot(Kinv, np.dot(dK_di, np.dot(Kinv, dy_dj))))
    terms[8] = np.dot(dy_dj.T, np.dot(Kinv, dy_di))
    terms[9] = -np.dot(y_arr.T, np.dot(Kinv, np.dot(dK_dj, np.dot(Kinv, dy_di))))
    terms[10] = np.dot(y_arr.T, np.dot(Kinv, dy_didj))

    return .5 * np.sum(terms)


def _keplerian(Krv, phase):
    return -Krv * np.sin(2*np.pi*phase)


def _intovector(x, size):
    '''
    Input from partial derivatives may be a 1D vector or a scalar. Change dimensions
    do future dot product calculations.
    '''
    if type(x) == int:  # x == 0
        return np.zeros((size, 1))
    else:
        return x.reshape(size, 1)
    
def _intomatrix(x, size):
    '''
    Input from partial derivatives may be a 1D vector or a scalar. Change dimensions
    do future dot product calculations.
    '''
    if type(x) == int:  # x == 0
        return np.zeros((size, size))
    else:
        return x
