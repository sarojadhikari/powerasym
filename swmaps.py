import numpy as np
import healpy as hp
from ngsachswolfe import SachsWolfeMap

NSIMS=10000
LMAX=100
NSIDE=64
datadir="datanodipole/"
NODIPOLE=True

fNL=500
gNL=5000000

#usesavedmaps=False
usesavedmaps=True

ratiosGprev=np.array([]); ratiosfNLprev=np.array([]); ratiosgNLprev=np.array([])
rGprev=np.array([]); rfNLprev=np.array([]); rgNLprev=np.array([])  
dGprev=np.array([]); dfNLprev=np.array([]); dgNLprev=np.array([])
ClGprev=np.array([]); ClfNLprev=np.array([]); ClgNLprev=np.array([])

if not(usesavedmaps):
    try:
        ratiosGprev=np.load(datadir+"AdistG.npy")
        ratiosfNLprev=np.load(datadir+"AdistfNL"+str(fNL)+".npy")
        ratiosgNLprev=np.load(datadir+"AdistgNL"+str(gNL)+".npy")
    
        rGprev=np.load(datadir+"AidistG.npy")
        rfNLprev=np.load(datadir+"AidistfNL"+str(fNL)+".npy")
        rgNLprev=np.load(datadir+"AidistgNL"+str(gNL)+".npy")
    
        dGprev=np.load(datadir+"dipolesG.npy")
        dfNLprev=np.load(datadir+"dipolesfNL"+str(fNL)+".npy")
        dgNLprev=np.load(datadir+"dipolesgNL"+str(gNL)+".npy")
    
        ClGprev=np.load(datadir+"ClsG.npy")
        ClfNLprev=np.load(datadir+"ClsfNL"+str(fNL)+".npy")
        ClgNLprev=np.load(datadir+"ClsgNL"+str(gNL)+".npy")
        
    except:
        print "file read error"
        
startN=len(ratiosGprev); print startN

ratiosG=[]; ratiosfNL=[]; ratiosgNL=[]
rG=[]; rfNL=[]; rgNL=[]
dG=[]; dfNL=[]; dgNL=[]
ClG=[]; ClfNL=[]; ClgNL=[]

for sim in range(NSIMS):
    if not(usesavedmaps):
        swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE)
    else:
        print "reading map"
        swmap=SachsWolfeMap(fnl=fNL, gnl=gNL, LMAX=LMAX, NSIDE=NSIDE, nodipole=NODIPOLE, readGmap=sim)
        
    swmap.generate_fnl_map()
    swmap.generate_gnl_map()
    swmap.calculate()
    
    ratiosG.append(swmap.gausA); ratiosfNL.append(swmap.fnlA); ratiosgNL.append(swmap.gnlA)
    rG.append(swmap.gausAi); rfNL.append(swmap.fnlAi); rgNL.append(swmap.gnlAi)
    dG.append(swmap.gausdipole); dfNL.append(swmap.fnldipole); dgNL.append(swmap.gnldipole)
    ClG.append(swmap.gausCls); ClfNL.append(swmap.fnlCls); ClgNL.append(swmap.gnlCls)
    
    #if not(usesavedmaps):
        #swmap.save_gaus_map(startN+sim)

# SAVE 

if not(usesavedmaps):
    np.save(datadir+"AdistG.npy", np.append(ratiosGprev, ratiosG))
np.save(datadir+"AdistfNL"+str(fNL)+".npy", np.append(ratiosfNLprev, ratiosfNL))
np.save(datadir+"AdistgNL"+str(gNL)+".npy", np.append(ratiosgNLprev, ratiosgNL))

if not(usesavedmaps):
    np.save(datadir+"AidistG.npy", np.append(rGprev, rG))
np.save(datadir+"AidistfNL"+str(fNL)+".npy", np.append(rfNLprev, rfNL))
np.save(datadir+"AidistgNL"+str(gNL)+".npy", np.append(rgNLprev, rgNL))

if not(usesavedmaps):
    np.save(datadir+"dipolesG.npy", np.append(dGprev, dG))
np.save(datadir+"dipolesfNL"+str(fNL)+".npy", np.append(dfNLprev, dfNL))
np.save(datadir+"dipolesgNL"+str(gNL)+".npy", np.append(dgNLprev, dgNL))
    
if not(usesavedmaps):
    np.save(datadir+"ClsG.npy", np.append(ClGprev, ClG))
np.save(datadir+"ClsfNL"+str(fNL)+".npy", np.append(ClfNLprev, ClfNL))
np.save(datadir+"ClsgNL"+str(gNL)+".npy", np.append(ClgNLprev, ClgNL))
