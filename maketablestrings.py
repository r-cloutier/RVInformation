from imports import *

def get_best50_table(self, N=50, sigma=5):
    '''create a latex table string with the best 50 small planets with either H or N'''
    # best planet parameters between either spectrograph
    rps = np.append(self.rps_med, self.rps_med)
    g = rps < 4
    rps = rps[g]
    assert g.size == 3968
    TOIs = np.append(self.starnums_med, self.starnums_med)[g]
    Ps = np.append(self.Ps_med, self.Ps_med)[g]
    mps = np.append(self.mps_med, self.mps_med)[g]
    Ks = np.append(self.Ks_med, self.Ks_med)[g]
    #Mss = np.append(self.Mss_med, self.Mss_med)[g]
    Teffs = np.append(self.Teffs_med, self.Teffs_med)[g]
    Vs = np.append(self.Vmags_med, self.Vmags_med)[g]
    Js = np.append(self.Jmags_med, self.Jmags_med)[g]
    #Nrvs = np.append(self.NrvGPs_med_H, self.NrvGPs_med_N)[g]
    corr = (.327/.189)**2 if sigma == 5 else 1
    tobs = np.append(self.tobsGPs_med_H, self.tobsGPs_med_N)[g] * corr
    spectrographs = np.append(np.repeat('opt',self.nstars), np.repeat('nIR',self.nstars))[g]
    
    # get best planets and their parameters
    s = np.argsort(tobs)[:int(N)]
    TOI50, P50, rp50, K50 = TOIs[s], Ps[s], rps[s], Ks[s]
    Teff50, V50, J50, tobs50 = Teffs[s], Vs[s], Js[s], tobs[s]
    spectrograph50 = spectrographs[s]
    print 'opt:', (spectrograph50=='opt').sum(), 'nIR:', (spectrograph50=='nIR').sum()

    # create string
    table = ''
    for i in range(int(N)):
	table += ' %.4d & %.3f & %.2f & %.2f & %i & %.2f & %.2f & %.2f & %s \\\\ \n'%(TOI50[i],P50[i],rp50[i], K50[i], 
										      Teff50[i],V50[i],J50[i],tobs50[i],
	 									      spectrograph50[i])

    # add to tex file
    f = open('paper/50random_template.tex', 'r')
    g = f.read()
    f.close()
    g = g.replace('<<table>>',table)
    h = open('paper/50random.tex', 'w')
    h.write(g)
    h.close()

    return table


def get_bestradvalley_table(self, N=10, sigma=5):
    '''create a latex table string with the best 50 small planets with either H or N'''
    # best planet parameters between either spectrograph
    rps = np.append(self.rps_med, self.rps_med)
    Ps = np.append(self.Ps_med, self.Ps_med)
    rp_upper = 10**(-0.09*np.log10(Ps)+.44)
    rp_lower = 10**(-0.09*np.log10(Ps)+.29)
    g1 = (rps <= 2) & (rps >= rp_lower) & (rps <= rp_upper)
    g2 = (rps > 2) & (rps >= rp_lower) & (rps <= rp_upper)

    # get planets on one side of the gap
    Ps, rps = Ps[g1], rps[g1]
    TOIs = np.append(self.starnums_med, self.starnums_med)[g1]
    Ks = np.append(self.Ks_med, self.Ks_med)[g1]
    Teffs = np.append(self.Teffs_med, self.Teffs_med)[g1]
    Vs = np.append(self.Vmags_med, self.Vmags_med)[g1]
    Js = np.append(self.Jmags_med, self.Jmags_med)[g1]
    corr = (.327/.189)**2 if sigma == 5 else 1
    tobs = np.append(self.tobsGPs_med_H, self.tobsGPs_med_N)[g1] * corr
    spectrographs = np.append(np.repeat('opt',self.nstars), np.repeat('nIR',self.nstars))[g1]
    
    # get best planets and their parameters for one side of the gap
    s = np.argsort(tobs)[:int(N)]
    TOI10, P10, rp10, K10 = TOIs[s], Ps[s], rps[s], Ks[s]
    Teff10, V10, J10, tobs10 = Teffs[s], Vs[s], Js[s], tobs[s]
    spectrograph10 = spectrographs[s]

    # get planets on one side of the gap
    rps = np.append(self.rps_med, self.rps_med)[g2]
    Ps = np.append(self.Ps_med, self.Ps_med)[g2]
    TOIs = np.append(self.starnums_med, self.starnums_med)[g2]
    Ks = np.append(self.Ks_med, self.Ks_med)[g2]
    Teffs = np.append(self.Teffs_med, self.Teffs_med)[g2]
    Vs = np.append(self.Vmags_med, self.Vmags_med)[g2]
    Js = np.append(self.Jmags_med, self.Jmags_med)[g2]
    tobs = np.append(self.tobsGPs_med_H, self.tobsGPs_med_N)[g2] * corr
    spectrographs = np.append(np.repeat('opt',self.nstars), np.repeat('nIR',self.nstars))[g2]
    
    # get best planets and their parameters for one side of the gap
    s = np.argsort(tobs)[:int(N)]
    TOI20, P20, rp20, K20 = np.append(TOI10,TOIs[s]), np.append(P10, Ps[s]), np.append(rp10,rps[s]), np.append(K10,Ks[s])
    Teff20, V20, J20, tobs20 = np.append(Teff10,Teffs[s]), np.append(V10,Vs[s]), np.append(J10,Js[s]), np.append(tobs10,tobs[s])
    spectrograph20 = np.append(spectrograph10,spectrographs[s])

    # sort by observing time
    s = np.argsort(tobs20)
    TOI20, P20, rp20, K20 = TOI20[s], P20[s], rp20[s], K20[s]
    Teff20, V20, J20, tobs20, spectrograph20 = Teff20[s], V20[s], J20[s], tobs20[s], spectrograph20[s] 
    
    # create string
    table = ''
    for i in range(int(N*2)):
	table += ' %.4d & %.3f & %.2f & %.2f & %i & %.2f & %.2f & %.2f & %s \\\\ \n'%(TOI20[i],P20[i],rp20[i], K20[i], #
										      Teff20[i],V20[i],J20[i],tobs20[i],
	 									      spectrograph20[i])

    # add to tex file
    f = open('paper/radgap_template.tex', 'r')
    g = f.read()
    f.close()
    g = g.replace('<<table>>',table)
    h = open('paper/radgap.tex', 'w')
    h.write(g)
    h.close()
    
    return table
