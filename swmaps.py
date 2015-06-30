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

fNL=100
gNL=10000

# test: 30000 has -3 while 20000 has -6 subtraction in gnl

usesavedmaps=False
#usesavedmaps=False

A0G=[]; A0fNL=[]; A0gNL=[]
AG=[]; AfNL=[]; AgNL=[]
AiG=[]; AifNL=[]; AigNL=[]
dG=[]; dfNL=[]; dgNL=[]
ClG=[]; ClfNL=[]; ClgNL=[]

phisq=[]
# first either generate all the Gaussian maps or get the maps and compute <phi^2>
for sim in range(NSIMS):
    if not(usesavedmaps):
        swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, mapsdir=mapsdir, N=NEFOLDS)
        swmap.save_gaus_map(sim)
    else:
        swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, readGmap=sim, mapsdir=mapsdir, N=NEFOLDS)
    
    phisq.append(np.mean(swmap.gausmap*swmap.gausmap))

print np.mean(phisq)
np.save(datadir+"phisq.npy", phisq)
phisqsub=np.mean(phisq) # global

#sys.exit(0)

# now generate non-Gaussian maps
for sim in range(NSIMS):
    swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, readGmap=sim, mapsdir=mapsdir, N=NEFOLDS)
    
    swmap.generate_fnl_map(phisqsub)
    swmap.generate_gnl_map(phisqsub)
    swmap.calculate()
    
    A0G.append(swmap.gausA0); A0fNL.append(swmap.fnlA0); A0gNL.append(swmap.gnlA0)
    AG.append(swmap.gausA); AfNL.append(swmap.fnlA); AgNL.append(swmap.gnlA)
    AiG.append(swmap.gausAi); AifNL.append(swmap.fnlAi); AigNL.append(swmap.gnlAi)
    dG.append(swmap.gausdipole); dfNL.append(swmap.fnldipole); dgNL.append(swmap.gnldipole)
    ClG.append(swmap.gausCls); ClfNL.append(swmap.fnlCls); ClgNL.append(swmap.gnlCls)
    
# SAVE 
if not(usesavedmaps):
    np.save(datadir+"A0distG.npy", A0G)
    np.save(datadir+"AdistG.npy", AG)
    np.save(datadir+"AidistG.npy", AiG)
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
    
np.save(datadir+"ClsfNL"+str(fNL)+".npy", ClfNL)
np.save(datadir+"ClsgNL"+str(gNL)+".npy", ClgNL)
