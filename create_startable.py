from RVInformationResults import *

header = '\clearpage\n\\begin{turnpage}\n\\begin{deluxetable}{ccccccccccccccccc}\n\\tabletypesize{\scriptsize}\n\\tablecaption{Stellar parameters from the \\citetalias{sullivan15} synthetic catalog \label{table:stars}}\n\\tablewidth{0pt}\n\\tablehead{TOI & $\\alpha$ & $\\delta$ & $P$ & $m_p$ & $K$ & $S$ & $M_s$ & $T_{\\text{eff}}$ & Distance & $B$ & $V$ & $Y$ & $J$ & $H$ & $\\text{med}(v\\sin{i_s})$ & \\\ & $[$deg$]$ & $[$deg$]$ & $[$days$]$ & $[\\text{M}_{\oplus}]$ & $[$m s$^{-1}]$ & $[\\text{S}_{\oplus}]$ & $[\\text{M}_{\odot}]$ & $[$K$]$ & $[$pc$]$ & & & & & & $[$km s$^{-1}]$}\n'

footer = '\enddata\n\end{deluxetable}\n\end{turnpage}\n\clearpage\n\global\pdfpageattr\expandafter{\\the\pdfpageattr/Rotate 90}'

# get data
#self = loadpickle('pickles/RVInformationGP_Nrvge10')

def create_startable(self):
    body = ''
    for i in range(self.nstars):
        prefix = '' if i < 10 else '%'
        suffix = '' if i == self.nstars-1 else '\n'
        theta = prefix, i, self.ras_med[i], self.decs_med[i], self.Ps_med[i], \
                self.mps_med[i], self.Ks_med[i], self.Fs_med[i], \
                self.Mss_med[i], self.Teffs_med[i], self.dists_med[i], \
                self.Bmags_med[i], self.Vmags_med[i], self.Ymags_med[i], \
                self.Jmags_med[i], self.Hmags_med[i], self.vsinis_med[i], suffix
        body += '%s%.4d & %.2f & %.2f & %.3f & %.2f & %.2f & %.1f & %.2f & %i & %.1f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\ %s'%theta

        
    # write table
    g = header + body + footer
    h = open('paper/startable.tex', 'w')
    h.write(g)
    h.close()
