import sys
import os
import numpy as np
import healpy as hp
from ngsachswolfe import SachsWolfeMap

NSIMS=1000
LMAX=100
NSIDE=128
NEFOLDS=40  #

mapsdir="maps"+str(NEFOLDS)+"/"
datadir="data"+str(NEFOLDS)+"/"

if not os.path.exists(datadir):
    os.makedirs(datadir)
if not os.path.exists(mapsdir):
    os.makedirs(mapsdir)
    
NODIPOLE=False

fNLlist=[100000, 500000, 1000000, 10000000]
NfNL=len(fNLlist)
#gNL=700000000

usesavedmaps=True
#usesavedmaps=False

A0G=[]; A0fNL=[[] for i in range(NfNL)]
AG=[]; AfNL=[[] for i in range(NfNL)]
AiG=[]; AifNL=[[] for i in range(NfNL)]

A2G=[]; A2fNL=[[] for i in range(NfNL)]
Ai2G=[]; Ai2fNL=[[] for i in range(NfNL)]

dG=[]; dfNL=[[] for i in range(NfNL)]
Cl0G=[]; Cl0fNL=[[] for i in range(NfNL)]
Cl1G=[]; Cl1fNL=[[] for i in range(NfNL)]

phisq0=[]; phisq1=[]

# first either generate all the Gaussian maps or get the maps and compute <phi^2>
if not(usesavedmaps):
    for sim in range(NSIMS):
        swmap=SachsWolfeMap(fnls=fNLlist, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, mapsdir=mapsdir, N=NEFOLDS)
        swmap.save_gaus_map(sim)
        phisq0.append(np.mean(swmap.gausmap0*swmap.gausmap0))
        phisq1.append(np.mean(swmap.gausmap1*swmap.gausmap1))
    np.save(datadir+"phisq0.npy", phisq0)
    np.save(datadir+"phisq1.npy", phisq1)
else:
    phisq0=np.load(datadir+"phisq0.npy")
    phisq1=np.load(datadir+"phisq1.npy")

print np.mean(phisq0), np.mean(phisq1)
phisqsub0=np.mean(phisq0); phisqsub1=np.mean(phisq1) # global

#sys.exit(0)

# now generate non-Gaussian maps
for sim in range(NSIMS):
    print sim
    swmap=SachsWolfeMap(fnls=fNLlist, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, readGmap=sim, mapsdir=mapsdir, N=NEFOLDS)
    
    swmap.generate_fnl_maps(phisqsub0, phisqsub1)
    #swmap.generate_gnl_map(phisqsub0, phisqsub1)
    swmap.calculate()

    A0G.append(swmap.gausA0)
    AG.append(swmap.gausA)
    AiG.append(swmap.gausAi)
    dG.append(swmap.gausdipole)
    A2G.append(swmap.gausA2)
    Ai2G.append(swmap.gausAi2)
    Cl0G.append(swmap.gausCls0)
    Cl1G.append(swmap.gausCls1)
    
    for i in range(NfNL):
        A0fNL[i].append(swmap.fnlA0[i])
        AfNL[i].append(swmap.fnlA[i])
        AifNL[i].append(swmap.fnlAi[i])
        dfNL[i].append(swmap.fnldipole[i])
        A2fNL[i].append(swmap.fnlA2[i])
        Ai2fNL[i].append(swmap.fnlAi2[i])
        Cl0fNL[i].append(swmap.fnlCls0[i])
        Cl1fNL[i].append(swmap.fnlCls1[i])
        
# SAVE 
np.save(datadir+"A0distG.npy", A0G)
np.save(datadir+"AdistG.npy", AG)
np.save(datadir+"AidistG.npy", AiG)
np.save(datadir+"Ai2distG.npy", Ai2G) #
np.save(datadir+"A2distG.npy", A2G) #
np.save(datadir+"dipolesG.npy", dG)

np.save(datadir+"Cls0G.npy", Cl0G)
np.save(datadir+"Cls1G.npy", Cl1G)

for i in range(NfNL):
    np.save(datadir+"A0distfNL"+str(fNLlist[i])+".npy", A0fNL[i])
    np.save(datadir+"AdistfNL"+str(fNLlist[i])+".npy", AfNL[i])
    np.save(datadir+"AidistfNL"+str(fNLlist[i])+".npy", AifNL[i])
    np.save(datadir+"dipolesfNL"+str(fNLlist[i])+".npy", dfNL[i])
    np.save(datadir+"Ai2distfNL"+str(fNLlist[i])+".npy", Ai2fNL[i])
    np.save(datadir+"A2distfNL"+str(fNLlist[i])+".npy", A2fNL[i])
    np.save(datadir+"Cls0fNL"+str(fNLlist[i])+".npy", Cl0fNL[i])
    np.save(datadir+"Cls1fNL"+str(fNLlist[i])+".npy", Cl1fNL[i])
