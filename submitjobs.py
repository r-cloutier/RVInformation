import os

nplanets = 1984
planetindices = range(nplanets)

for planetindex in planetindices:
    f = open('jobscript_template', 'r')
    g = f.read()
    f.close()

    g = g.replace('<<planetindex>>', '%i'%planetindex)

    h = open('jobscript', 'w')
    h.write(g)
    h.close()

    #os.system('qsub jobscript')
    #os.system('rm jobscript')
    os.system('cat jobscript')
