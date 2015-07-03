import sys
import os
import numpy as np
import healpy as hp
from ngsachswolfe import SachsWolfeMap

NSIMS=10000
LMAX=100
NSIDE=64
NEFOLDS=50  #

mapsdir="maps"+str(NEFOLDS)+"/"
datadir="data"+str(NEFOLDS)+"/"

if not os.path.exists(datadir):
    os.makedirs(datadir)
if not os.path.exists(mapsdir):
    os.makedirs(mapsdir)
    
NODIPOLE=False

fNL=250
gNL=10000

usesavedmaps=False
#usesavedmaps=False

A0G=[]; A0fNL=[]; A0gNL=[]
AG=[]; AfNL=[]; AgNL=[]
AiG=[]; AifNL=[]; AigNL=[]

A2G=[]; A2fNL=[]; A2gNL=[]
Ai2G=[]; Ai2fNL=[]; Ai2gNL=[]

dG=[]; dfNL=[]; dgNL=[]
ClG=[]; ClfNL=[]; ClgNL=[]

phisq0=[]
phisq1=[]
# first either generate all the Gaussian maps or get the maps and compute <phi^2>
for sim in range(NSIMS):
    if not(usesavedmaps):
        swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, mapsdir=mapsdir, N=NEFOLDS)
        swmap.save_gaus_map(sim)
    else:
        swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, readGmap=sim, mapsdir=mapsdir, N=NEFOLDS)
    
    phisq0.append(np.mean(swmap.gausmap0*swmap.gausmap0))
    phisq1.append(np.mean(swmap.gausmap1*swmap.gausmap1))

print np.mean(phisq0), np.mean(phisq1)
np.save(datadir+"phisq0.npy", phisq0)
np.save(datadir+"phisq1.npy", phisq1)
phisqsub0=np.mean(phisq0); phisqsub1=np.mean(phisq1) # global

#sys.exit(0)

# now generate non-Gaussian maps
for sim in range(NSIMS):
    print sim
    swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, readGmap=sim, mapsdir=mapsdir, N=NEFOLDS)
    
    swmap.generate_fnl_map(phisqsub0, phisqsub1)
    swmap.generate_gnl_map(phisqsub0, phisqsub1)
    swmap.calculate()
    
    A0G.append(swmap.gausA0); A0fNL.append(swmap.fnlA0); A0gNL.append(swmap.gnlA0)
    AG.append(swmap.gausA); AfNL.append(swmap.fnlA); AgNL.append(swmap.gnlA)
    AiG.append(swmap.gausAi); AifNL.append(swmap.fnlAi); AigNL.append(swmap.gnlAi)
    dG.append(swmap.gausdipole); dfNL.append(swmap.fnldipole); dgNL.append(swmap.gnldipole)

    A2G.append(swmap.gausA2); A2fNL.append(swmap.fnlA2); A2gNL.append(swmap.gnlA2)
    Ai2G.append(swmap.gausAi2); Ai2fNL.append(swmap.fnlAi2); Ai2gNL.append(swmap.gnlAi2)

    #ClG.append(swmap.gausCls0); ClfNL.append(swmap.fnlCls); ClgNL.append(swmap.gnlCls)
    
# SAVE 
if not(usesavedmaps):
    np.save(datadir+"A0distG.npy", A0G)
    np.save(datadir+"AdistG.npy", AG)
    np.save(datadir+"AidistG.npy", AiG)
    np.save(datadir+"Ai2distG.npy", Ai2G) #
    np.save(datadir+"A2distG.npy", A2G) #
    np.save(datadir+"dipolesG.npy", dG)
    np.save(datadir+"ClsG.npy", ClG)

np.save(datadir+"A0distfNL"+str(fNL)+".npy", A0fNL)
np.save(datadir+"A0distgNL"+str(gNL)+".npy", A0gNL)

np.save(datadir+"AdistfNL"+str(fNL)+".npy", AfNL)
np.save(datadir+"AdistgNL"+str(gNL)+".npy", AgNL)

np.save(datadir+"AidistfNL"+str(fNL)+".npy", AifNL)
np.save(datadir+"AidistgNL"+str(gNL)+".npy", AigNL)

np.save(datadir+"dipolesfNL"+str(fNL)+".npy", dfNL)
np.save(datadir+"dipolesgNL"+str(gNL)+".npy", dgNL)

np.save(datadir+"Ai2distfNL"+str(fNL)+".npy", Ai2fNL)
np.save(datadir+"Ai2distgNL"+str(gNL)+".npy", Ai2gNL)

np.save(datadir+"A2distfNL"+str(fNL)+".npy", A2fNL)
np.save(datadir+"A2distgNL"+str(gNL)+".npy", A2gNL)
    
#np.save(datadir+"ClsfNL"+str(fNL)+".npy", ClfNL)
#np.save(datadir+"ClsgNL"+str(gNL)+".npy", ClgNL)
