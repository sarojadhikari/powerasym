import os
import numpy as np
import healpy as hp
from ngsachswolfe import SachsWolfeMap

NSIMS=5000
LMAX=100
NSIDE=64
NEFOLDS=50

mapsdir="maps"+str(NEFOLDS)+"/"
datadir="data"+str(NEFOLDS)+"/"

if not os.path.exists(datadir):
    os.makedirs(datadir)
if not os.path.exists(mapsdir):
    os.makedirs(mapsdir)
    
NODIPOLE=True

fNL=500
gNL=5000000

usesavedmaps=False
#usesavedmaps=True

A0Gprev=np.array([]); A0fNLprev=np.array([]); A0gNLprev=np.array([])
AGprev=np.array([]); AfNLprev=np.array([]); AgNLprev=np.array([])
AiGprev=np.array([]); AifNLprev=np.array([]); AigNLprev=np.array([])  
dGprev=np.array([]); dfNLprev=np.array([]); dgNLprev=np.array([])
ClGprev=np.array([]); ClfNLprev=np.array([]); ClgNLprev=np.array([])

if not(usesavedmaps):
    try:
        A0Gprev=np.load(datadir+"A0distG.npy")
        A0fNLprev=np.load(datadir+"A0distfNL"+str(fNL)+".npy")
        A0gNLprev=np.load(datadir+"A0distgNL"+str(gNL)+".npy")
        
        AGprev=np.load(datadir+"AdistG.npy")
        AfNLprev=np.load(datadir+"AdistfNL"+str(fNL)+".npy")
        AgNLprev=np.load(datadir+"AdistgNL"+str(gNL)+".npy")
    
        AiGprev=np.load(datadir+"AidistG.npy")
        AifNLprev=np.load(datadir+"AidistfNL"+str(fNL)+".npy")
        AigNLprev=np.load(datadir+"AidistgNL"+str(gNL)+".npy")
    
        dGprev=np.load(datadir+"dipolesG.npy")
        dfNLprev=np.load(datadir+"dipolesfNL"+str(fNL)+".npy")
        dgNLprev=np.load(datadir+"dipolesgNL"+str(gNL)+".npy")
    
        ClGprev=np.load(datadir+"ClsG.npy")
        ClfNLprev=np.load(datadir+"ClsfNL"+str(fNL)+".npy")
        ClgNLprev=np.load(datadir+"ClsgNL"+str(gNL)+".npy")
        
    except:
        print "file read error"
        
startN=len(ratiosGprev); print startN

A0G=[]; A0fNL=[]; A0gNL=[]
AG=[]; AfNL=[]; AgNL=[]
AiG=[]; AifNL=[]; AigNL=[]
dG=[]; dfNL=[]; dgNL=[]
ClG=[]; ClfNL=[]; ClgNL=[]

for sim in range(NSIMS):
    if not(usesavedmaps):
        swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, mapsdir=mapsdir)
    else:
        print "reading map"
        swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, readGmap=sim, mapsdir=mapsdir)
        
    swmap.generate_fnl_map()
    swmap.generate_gnl_map()
    swmap.calculate()
    
    A0G.append(swmap.gausA0); A0fNL.append(swmap.fnlA0); A0gNL.append(swmap.gnlA0)
    AG.append(swmap.gausA); AfNL.append(swmap.fnlA); AgNL.append(swmap.gnlA)
    AiG.append(swmap.gausAi); AifNL.append(swmap.fnlAi); AigNL.append(swmap.gnlAi)
    dG.append(swmap.gausdipole); dfNL.append(swmap.fnldipole); dgNL.append(swmap.gnldipole)
    ClG.append(swmap.gausCls); ClfNL.append(swmap.fnlCls); ClgNL.append(swmap.gnlCls)
    
    if not(usesavedmaps):
        swmap.save_gaus_map(startN+sim)

# SAVE 
if not(usesavedmaps):
    np.save(datadir+"A0distG.npy", np.append(A0Gprev, A0G))
    np.save(datadir+"AdistG.npy", np.append(AGprev, AG))
    np.save(datadir+"AidistG.npy", np.append(AiGprev, AiG))
    np.save(datadir+"dipolesG.npy", np.append(dGprev, dG))
    np.save(datadir+"ClsG.npy", np.append(ClsGprev, ClsG))

np.save(datadir+"A0distfNL"+str(fNL)+".npy", np.append(A0fNLprev, A0fNL))
np.save(datadir+"A0distgNL"+str(gNL)+".npy", np.append(A0gNLprev, A0gNL))

np.save(datadir+"AdistfNL"+str(fNL)+".npy", np.append(AfNLprev, AfNL))
np.save(datadir+"AdistgNL"+str(gNL)+".npy", np.append(AgNLprev, AgNL))

np.save(datadir+"AidistfNL"+str(fNL)+".npy", np.append(AifNLprev, AifNL))
np.save(datadir+"AidistgNL"+str(gNL)+".npy", np.append(AigNLprev, AigNL))

np.save(datadir+"dipolesfNL"+str(fNL)+".npy", np.append(dfNLprev, dfNL))
np.save(datadir+"dipolesgNL"+str(gNL)+".npy", np.append(dgNLprev, dgNL))
    
np.save(datadir+"ClsfNL"+str(fNL)+".npy", np.append(ClfNLprev, ClfNL))
np.save(datadir+"ClsgNL"+str(gNL)+".npy", np.append(ClgNLprev, ClgNL))
