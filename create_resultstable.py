from RVInformationResults import *

header = '\clearpage\n\\begin{turnpage}\n\\begin{deluxetable}{ccccccccccc}\n\\tabletypesize{\scriptsize}\n\\tablecaption{Median radial velocity follow-up calculations for the \\citetalias{sullivan15} synthetic catalog including a GP red noise model \label{table:results}}\n\\tablewidth{0pt}\n\\tablehead{TOI & $ \sigma_{\\text{RV},\\text{opt}} $ & $ \sigma_{\\text{RV},\\text{nIR}} $ & $ \sigma_{\\text{act}} $ & $ \sigma_{\\text{planets}} $ & $ \sigma_{\\text{eff,opt}} $ & $ \sigma_{\\text{eff,nIR}} $ & $ n_{\\text{RV,opt}} $ & $ n_{\\text{RV,nIR}} $ & $ t_{\\text{obs,opt}} $ & $ t_{\\text{obs,nIR}} $ \\\ & $[$m s$^{-1}]$ & $[$m s$^{-1}]$ & $[$m s$^{-1}]$ & $[$m s$^{-1}]$ & $[$m s$^{-1}]$ & $[$m s$^{-1}]$ & & & [nights] & [nights]}\n'

footer = '\enddata\n\end{deluxetable}\n\end{turnpage}\n\clearpage\n\global\pdfpageattr\expandafter{\\the\pdfpageattr/Rotate 90}'

# get data
#self = loadpickle('pickles/RVInformationGP_Nrvge10')

def create_resultstable(self):
    body = ''
    for i in range(self.nstars):
        prefix = '' if i < 10 else '%'
        suffix = '' if i == self.nstars-1 else '\n'
        theta = prefix, i, self.sigmaRV_phot_med_H[i], \
                self.sigmaRV_phot_med_N[i], self.sigmaRV_acts_med_N[i], \
                self.sigmaRV_planets_med_N[i], self.sigmaRV_eff_med_H[i], \
                self.sigmaRV_eff_med_N[i], self.NrvGPs_med_H[i], \
                self.NrvGPs_med_N[i], self.tobsGPs_med_H[i]/7., \
                self.tobsGPs_med_N[i]/7., suffix
        body += '%s%.4d & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.1f & %.1f & %.1f & %.1f \\\\ %s'%theta
        
    # write table
    g = header + body + footer
    h = open('paper/resultstable.tex', 'w')
    h.write(g)
    h.close()
