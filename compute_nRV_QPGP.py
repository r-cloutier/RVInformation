from imports import *
from sympy.matrices import *
from sympy import *
from rvmodel import get_rv1

init_printing()

#t = np.linspace(0,100,101)
#y = np.sin(t)
#ey = np.repeat(.1, t.size)
#theta = 1, 100, .5, 10


def covariance_matrix(theta, t, sig=np.zeros(0)):
    a, l, G, P = theta
    assert t.shape == (t.size, 1)
    K = np.zeros ((t.size, t.size))
    K += -.5/l**2 * (np.tile(t[:,0], t.shape) - \
                     np.tile(t[:,0], t.shape).T)**2  # :)
    K += -G**2 * np.sin(np.pi/P * abs(np.tile(t[:,0], t.shape) - \
                                      np.tile(t[:,0], t.shape).T))**2
    K = a*a*np.exp(K)
    if sig.size == t.size:
        return K + sig**2 * np.eye(t.size)
    else:
        return K


def lnlikeold(N=100):
    t = MatrixSymbol('t', N, 1)
    y = MatrixSymbol('y', N, 1)
    a, l, G, P = symbols('a l G P')
    K = -.5*(t-t)**2 / l**2
    #lnl = -.5 * (Transpose(y) * K.inv('LU') * y) 
                

# try sympy
#t, a, l, G, P = symbols('t a l G P')
#diff(lnl(t


def partial_lnlike(theta, t, y, ey):
    '''See Eq 5.9 in Rasmussen & Williams to compute the partial 
    derivative of the lnlikelihood wrt the GP hyperparameters.'''
    t, y, ey = t.reshape(t.size,1), y.reshape(t.size,1), ey.reshape(t.size,1)
    K = covariance_matrix(theta, t, ey)
    Kinv = np.linalg.inv(K)
    alpha = np.dot(Kinv, y)
    #return .5 * np.trace((np.dot(alpha, alpha.T) - Kinv))





def get_timeseries():
    tarr = np.random.rand(100) * 100
    P, T0, K, sig = 10., 0., 3., 1.
    rvarr = get_rv1((P,T0,0,K,0,0), tarr) + np.random.randn(tarr.size)*sig
    ervarr = np.repeat(sig, tarr.size)
    return tarr, rvarr, ervarr


def compute_Fisher_information_GP(theta, t_arr, rv_arr, erv_arr):
    '''Compute the Fisher information matrix with very complex terms for a 
    circular keplerian RV model and a QP GP red noise model.
    theta = list of parameter values (P, T0, Krv, a, l, G, Pgp)
    '''
    # get orbital phase array
    P, T0, Krv, a, l, G, Pgp = theta
    phase_arr = foldAt(t_arr, P, T0)
    thetavals = theta[2:]
    
    # define variables
    Krv, a, l, G, P = symbols('Krv a lambda Gamma P')
    symbol_vals = Krv, a, l, G, P

    # define time-series symbols
    t, phi, rv, sig = symbols('t phi RV sigma')
    y = rv - (-Krv*sympy.sin(phi))
    
    # compute covariance function
    k = a*a*sympy.exp(-.5*t**2 / l**2 - G*G*sympy.sin(np.pi*abs(t)/P)**2)
    symbol_arrs = t, phi, rv, y, k

    # get arrays
    Kinv = np.linalg.inv(covariance_matrix(theta[3:], t_arr, erv_arr))
    thetaarrs = t_arr, phase_arr, rv_arr, Kinv

    # compute Fisher matrix
    Nparams = len(theta)
    B = np.zeros((Nparams, Nparams))
    for i in range(Nparams):
        for j in range(Nparams):
            if i > j:
                B[i,j] = B[j,i]
            else:
                B[i,j] = _compute_Fisher_entry(theta[i], theta[j], symbol_vals,
                                               symbol_arrs, thetavals, thetaarrs)
    return B
                

def _compute_Fisher_entry(symbol_i, symbol_j, symbols_values, symbol_arrays,
                          thetavals, thetaarrs):
    '''
    Compute the Fisher information entry for one pair of model 
    parameters. I.e. the partial of the lnlikelihood wrt to each 
    symbol.
    thetavals
    tarr = array of observation epochs (i.e. window function)
    y = RV residual symbol
    '''
    # compute partial expressions
    Krv_sym, a_sym, l_sym, G_sym, P_sym = symbol_values
    t_sym, phi_sym, rv_sym, y_sym, K_sym = symbol_arrays
    dy_didj = lambdify([Krv_sym, phi_sym, rv_sym],
                       diff(y_sym, symbol_i, symbol_j), 'numpy')
    dy_di = lambdify([Krv_sym, phi_sym, rv_sym],
                     diff(y_sym, symbol_i), 'numpy')
    dy_dj = lambdify([Krv_sym, phi_sym, rv_sym],
                     diff(y_sym, symbol_j), 'numpy')
    dK_didj = lambdify([a_sym, l_sym, G_sym, P_sym, t_sym],
                       diff(K_sym, symbol_i, symbol_j), 'numpy')
    dK_di = lambdify([a_sym, l_sym, G_sym, P_sym, t_sym],
                     diff(K_sym, symbol_i), 'numpy')
    dK_dj = lambdify([a_sym, l_sym, G_sym, P_sym, t_sym],
                     diff(K_sym, symbol_j), 'numpy')

    # evaluate partials at input values
    K_val, a_val, l_val, G_val, P_val = thetavals
    t_arr, phase_arr, rv_arr, Kinv =  thetaarrs
    y_arr = rv_arr - keplerian(K_val, phase_arr)
    dy_didj = dy_didj(K_val, phase_arr, rv_arr)
    dy_di = dy_di(K_val, phase_arr, rv_arr)
    dy_dj = dy_dj(K_val, phase_arr, rv_arr)
    dK_didj = dK_didj(a_val, l_val, G_val, P_val, t_arr)
    dK_di = dK_di(a_val, l_val, G_val, P_val, t_arr)
    dK_dj = dK_dj(a_val, l_val, G_val, P_val, t_arr)
    
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


def keplerian(Krv, phase):
    return -Krv * np.sin(phase)

    
def _residual_double_partials(theta, y):
    #Krv, a, l, G, P = theta
    Nparams = len(theta)
    matrix = np.zeros((Nparams, Nparams))
    for i in range(Nparams):
        for j in range(Nparams):
            if i > j:
                matrix[i,j] = matrix[j,i]
            else:
                matrix[i,j] = diff(y, theta[i], theta[j])
    # this matrix is all zeros so no need to run the loop
    return matrix


def _covariance_double_partials(theta, k):
    Nparams = len(theta)
    matrix = np.zeros((Nparams, Nparams))
    for i in range(Nparams):
        for j in range(Nparams):
            if i > j:
                #matrix[i,j] = matrix[j,i]
                print i,j,diff(k, theta[i], theta[j])
            else:
                #matrix[i,j] = diff(k, theta[i], theta[j])
                print i,j,diff(k, theta[i], theta[j])
    # this matrix is all zeros so no need to run the loop
    return matrix


def _covariance_single_partials(theta, k):
    Nparams = len(theta)
    vector = np.zeros((Nparams, 1))
    for i in range(Nparams):
        vector[i,0] = diff(k, theta[i])
    return vector


def _residual_single_partials(theta, y):
    Nparams = len(theta)
    vector = np.zeros((Nparams, 1))
    for i in range(Nparams):
        vector[i,0] = diff(y, theta[i])
    return vector
