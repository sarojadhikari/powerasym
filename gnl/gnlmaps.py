import sys
import os
import numpy as np
import healpy as hp
from ngSWgNL import gNLSWMap

SIMStart=0
NSIMS=5000
LMAX=100
NSIDE=128
NEFOLDS=150  #

mapsdir="maps50/"
datadir="gnldata"+str(NEFOLDS)+"/"

if not os.path.exists(datadir):
    os.makedirs(datadir)
if not os.path.exists(mapsdir):
    os.makedirs(mapsdir)
    
gNLlist=[1000, 10000, 100000, 1000000]
NgNL=len(gNLlist)

usesavedmaps=True

A0G=[]; A0gNL=[[] for i in range(NgNL)]
AG=[]; AgNL=[[] for i in range(NgNL)]
AiG=[]; AigNL=[[] for i in range(NgNL)]

dG=[]; dgNL=[[] for i in range(NgNL)]
ClG=[]; ClgNL=[[] for i in range(NgNL)]

phisq=[]

# first either generate all the Gaussian maps or get the maps and compute <phi^2>
if not(usesavedmaps):
    #print "use saved maps for gNL"
    for sim in range(SIMStart, NSIMS):
        swmap=gNLSWMap(gnls=gNLlist, LMAX=LMAX, NSIDE=NSIDE, mapsdir=mapsdir, N=NEFOLDS)
        swmap.save_gaus_map(sim)
        phisq.append(np.mean(swmap.gausmap*swmap.gausmap))
    np.save(mapsdir+"phisq.npy", phisq)
else:
    for sim in range(0, 10000):
        swmap=gNLSWMap(gnls=gNLlist, LMAX=LMAX, NSIDE=NSIDE, readGmap=sim, mapsdir=mapsdir, N=NEFOLDS)
        phisq.append(np.mean(swmap.gausmap*swmap.gausmap))

np.save(datadir+"phisq.npy")
print np.mean(phisq)
phisqsub=np.mean(phisq)

sys.exit(0)

# now generate non-Gaussian maps
if (SIMStart>0):
    try:
        A0G=np.load(datadir+"A0distG.npy")[:SIMStart]
        AG=np.load(datadir+"AdistG.npy")[:SIMStart]
        AiG=np.load(datadir+"AidistG.npy")[:SIMStart]
        dG=np.load(datadir+"dipolesG.npy")[:SIMStart]
        ClG=np.load(datadir+"ClsG.npy")[:SIMStart]
    
        for i in range(NgNL):
            A0gNL[i]=np.load(datadir+"A0distgNL"+str(gNLlist[i])+".npy")[:SIMStart]
            AgNL[i]=np.load(datadir+"AdistgNL"+str(gNLlist[i])+".npy")[:SIMStart]
            AigNL[i]=np.load(datadir+"AidistgNL"+str(gNLlist[i])+".npy")[:SIMStart]
            dgNL[i]=np.load(datadir+"dipolesgNL"+str(gNLlist[i])+".npy")[:SIMStart]
            ClgNL[i]=np.load(datadir+"ClsgNL"+str(gNLlist[i])+".npy")[:SIMStart]
    except:
        print "could not read all the saved files"

for sim in range(SIMStart, NSIMS):
    print sim
    swmap=gNLSWMap(gnls=gNLlist, LMAX=LMAX, NSIDE=NSIDE, readGmap=sim, mapsdir=mapsdir, N=NEFOLDS)    
    swmap.generate_gnl_maps(phisqsub)
    swmap.calculate()

    A0G=np.append(A0G, swmap.gausA0)
    AG=np.append(AG, swmap.gausA)
    AiG=np.append(AiG, swmap.gausAi)
    dG=np.append(dG, swmap.gausdipole)
    ClG=np.append(ClG, swmap.gausCls)
    
    for i in range(NgNL):
        A0gNL[i]=np.append(A0gNL[i], swmap.gnlA0[i])
        AgNL[i]=np.append(AgNL[i], swmap.gnlA[i])
        AigNL[i]=np.append(AigNL[i], swmap.gnlAi[i])
        dgNL[i]=np.append(dgNL[i], swmap.gnldipole[i])
        ClgNL[i]=np.append(ClgNL[i], swmap.gnlCls[i])
        
# SAVE 
np.save(datadir+"A0distG.npy", A0G)
np.save(datadir+"AdistG.npy", AG)
np.save(datadir+"AidistG.npy", AiG)
np.save(datadir+"dipolesG.npy", dG)
np.save(datadir+"ClsG.npy", ClG)

for i in range(NgNL):
    np.save(datadir+"A0distgNL"+str(gNLlist[i])+".npy", A0gNL[i])
    np.save(datadir+"AdistgNL"+str(gNLlist[i])+".npy", AgNL[i])
    np.save(datadir+"AidistgNL"+str(gNLlist[i])+".npy", AigNL[i])
    np.save(datadir+"dipolesgNL"+str(gNLlist[i])+".npy", dgNL[i])
    np.save(datadir+"ClsgNL"+str(gNLlist[i])+".npy", ClgNL[i])
