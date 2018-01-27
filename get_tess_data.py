import numpy as np
import pylab as plt
from matplotlib.colors import LogNorm

# Read in the TESS simulation data from Sullivan et al. 2015
def get_TESS_data():
    '''--------------------------------------------------------------------------------
    Bytes Format Units  Label    Explanations                              Column Index
    --------------------------------------------------------------------------------
    1-  8 F8.4   deg    RAdeg    Right Ascension in decimal degrees (J2000) 0
    10- 18 F9.5   deg    DEdeg    Declination in decimal degrees (J2000) 1
    20- 26 F7.4   Rgeo   Rp       Planetary radius in Earth units 2
    28- 35 F8.4   d      Per      Period 3
    37- 45 E9.3   ---    S        Planetary insolation in Earth units 4
    47- 53 F7.4   m/s    K        Radial-velocity semiamplitude 5
    55- 60 F6.4   solRad R*       Stellar radius 6
    62- 70 E9.3   K      Teff     Stellar effective temperature 7
    72- 77 F6.3   mag    Vmag     Apparent V band magnitude 8
    79- 84 F6.3   mag    Imag     Apparent I_C_ band magnitude 9
    86- 91 F6.3   mag    Jmag     Apparent J band magnitude 10
    93- 98 F6.3   mag    Ksmag    Apparent K_s_ band magnitude 11
    100-104 F5.2   mag    DM       Distance Modulus 12 
    106-111 F6.3   ---    Dil      Dilution parameter for all detected eclipses 13
    113-118 F6.3   [-]    logsigV  Log standard deviation in relative flux (1) 14
    120-127 F8.3   ---    SNR      Signal-to-Noise Ratio 15
    129 I1     ---    Mult     Number of other planets orbiting the star 16
    --------------------------------------------------------------------------------'''
    d = np.loadtxt('input_data/apj516470t6_mrt.txt',skiprows=32)
    ra, dec, rp, per, F, K, Rs, Teff, Vmag, Imag, Jmag, Kmag, DM, dil, logsigV, snr, mult = d.T
    mult[mult == 0] = 1
    dist = 10**(1+DM/5.)
    return ra, dec, rp, per, F, K, Rs, Teff, Vmag, Imag, Jmag, Kmag, dist, dil, logsigV, snr, mult


def plot_Fdist(label=False, pltt=True):
    _,_,rp,P,S,K,_,Teff,_,_,_,_,dist,_,_,_,_ = get_TESS_data()
    good = Teff <= 38e2
    rp,P,S,K,dist=rp[good],P[good],S[good],K[good],dist[good]
    mp = np.zeros(rp.size)
    for i in range(rp.size):
    	mp[i] = rvs.rad2mass2(rp[i], mean_only=True)
    G, sigma, A = 6.67e-11, 5.67e-8, .3
    #Ms = np.sqrt(2*np.pi*G/rvs.days2sec(P)*(rvs.Mearth2kg(mp)/K)**3)
    #sma = rvs.semimajoraxis(P, Ms, mp)
    Teq = ((.25*S*1367.*(1-A)) / sigma)**(.25)
    
    fig = plt.figure(figsize=(6.5,5.8))
    ax = fig.add_subplot(111)
    minsize, maxsize = 1, 100
    #p = np.poly1d(np.polyfit([np.log10(K).min(),np.log10(K).max()],[minsize,maxsize],1))
    #cax = ax.scatter(dist, S, s=p(np.log10(K)), c='k', alpha=.2, edgecolors='none')
    #cbar = plt.colorbar(cax, orientation='horizontal', pad=.1)
    ax.plot([0,1e4], np.repeat(2,2), 'b--', lw=2)
    ax.plot([0,1e4], np.repeat(.2,2), 'b--', lw=2)

    bins = 15
    xbins, ybins = np.logspace(np.log10(dist.min()),np.log10(dist.max()),bins), np.logspace(np.log10(S.min()),np.log10(S.max()),bins)
    H, xedges, yedges = np.histogram2d(dist, S, bins=[xbins,ybins])
    cax = ax.pcolormesh(xedges, yedges, H.T, cmap=plt.get_cmap('Reds'), norm=LogNorm()) 
    ax.fill_between([1,1e3], .2,2, color='b', alpha=.1)
    cbar = plt.colorbar(cax, orientation='horizontal', pad=.11, ticks=np.append(range(1,11),range(20,50,10)))
    cbar.ax.set_xticklabels(['1','','','4','','','','','','10','','','40'])
    cbar.set_label('Number')
    
    ax.set_xscale('log'), ax.set_yscale('log')
    ax.set_xlabel('Distance (pc)'), ax.set_ylabel('Insolation ($S_{\oplus}$)')
    ax.set_xticks(np.logspace(0,3,4))
    ax.set_xticklabels(['%i'%i for i in np.logspace(0,3,4)])
    ax.set_yticks(np.logspace(-2,3,6))
    ax.set_yticklabels(['0.01','0.01','1','10','100','1000'])
    ax.set_xlim((1,1e3)), ax.set_ylim((1e-2,1e3))
    plt.subplots_adjust(left=.15, bottom=0, right=.94, top=.95)
    if label:
        plt.savefig('plots/Fdist_2dhist.png')
    if pltt:
        plt.show()
    plt.show()
    


if __name__ == '__main__':
    ra, dec, rp, per, S, K, Rs, Teff, Vmag, Imag, Jmag, Kmag, dist, dil, logsigV, snr, mult = get_TESS_data()
    #good = Teff <= 38e2
    #plot_Fdist()
