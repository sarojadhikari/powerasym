import sys
import os
import numpy as np
import healpy as hp
from ngSWgNL import gNLSWMap

NSIMS=10
LMAX=100
NSIDE=128
NEFOLDS=40  #

mapsdir="maps"+str(NEFOLDS)+"/"
datadir="gnldata"+str(NEFOLDS)+"/"

if not os.path.exists(datadir):
    os.makedirs(datadir)
if not os.path.exists(mapsdir):
    os.makedirs(mapsdir)
    
gNLlist=[1000, 10000, 50000, 100000]
NfNL=len(gNLlist)

usesavedmaps=True

A0G=[]; A0gNL=[[] for i in range(NfNL)]
AG=[]; AgNL=[[] for i in range(NfNL)]
AiG=[]; AigNL=[[] for i in range(NfNL)]

dG=[]; dgNL=[[] for i in range(NfNL)]
ClG=[]; ClgNL=[[] for i in range(NfNL)]

phisq=[]

# first either generate all the Gaussian maps or get the maps and compute <phi^2>
if not(usesavedmaps):
    print "use saved maps for gNL"
    #for sim in range(NSIMS):
        #swmap=SachsWolfeMap(fnls=fNLlist, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, mapsdir=mapsdir, N=NEFOLDS)
        #swmap.save_gaus_map(sim)
        #phisq0.append(np.mean(swmap.gausmap0*swmap.gausmap0))
        #phisq1.append(np.mean(swmap.gausmap1*swmap.gausmap1))
    #np.save(datadir+"phisq0.npy", phisq0)
    #np.save(datadir+"phisq1.npy", phisq1)
else:
    phisq=np.load(datadir+"phisq.npy")

print np.mean(phisq)
phisqsub=np.mean(phisq)

#sys.exit(0)

# now generate non-Gaussian maps
for sim in range(NSIMS):
    print sim
    swmap=gNLSWMap(gnls=gNLlist, LMAX=LMAX, NSIDE=NSIDE, readGmap=sim, mapsdir=mapsdir, N=NEFOLDS)    
    swmap.generate_gnl_maps(phisqsub)
    swmap.calculate()

    A0G.append(swmap.gausA0)
    AG.append(swmap.gausA)
    AiG.append(swmap.gausAi)
    dG.append(swmap.gausdipole)
    ClG.append(swmap.gausCls)
    
    for i in range(NfNL):
        A0gNL[i].append(swmap.gnlA0[i])
        AgNL[i].append(swmap.gnlA[i])
        AigNL[i].append(swmap.gnlAi[i])
        dgNL[i].append(swmap.gnldipole[i])
        ClgNL[i].append(swmap.gnlCls[i])
        
# SAVE 
np.save(datadir+"A0distG.npy", A0G)
np.save(datadir+"AdistG.npy", AG)
np.save(datadir+"AidistG.npy", AiG)
np.save(datadir+"dipolesG.npy", dG)
np.save(datadir+"ClsG.npy", ClG)

for i in range(NfNL):
    np.save(datadir+"A0distgNL"+str(gNLlist[i])+".npy", A0gNL[i])
    np.save(datadir+"AdistgNL"+str(gNLlist[i])+".npy", AgNL[i])
    np.save(datadir+"AidistgNL"+str(gNLlist[i])+".npy", AigNL[i])
    np.save(datadir+"dipolesgNL"+str(gNLlist[i])+".npy", dgNL[i])
    np.save(datadir+"ClsgNL"+str(gNLlist[i])+".npy", ClgNL[i])
