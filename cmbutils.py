import healpy as hp
import numpy as np

def get_hem_Cls(skymap, direction, LMAX=256):
    """
    from the given healpix skymap, return Cls for two hemispheres defined by the
    direction given, useful to study the possible scale dependence of power modulation
    
    direction should be a unit vector
    """
    # generate hemispherical mask
    NPIX=len(skymap)
    NSIDE=hp.npix2nside(NPIX)
    maskp=np.array([0.]*NPIX)
    disc=hp.query_disc(nside=NSIDE, vec=direction, radius=0.0174532925*90.)
    maskp[disc]=1.
    map1=hp.ma(skymap)
    map1.mask=maskp
    Clsp=hp.anafast(map1, lmax=LMAX)
    map1.mask=np.logical_not(maskp)
    Clsm=hp.anafast(map1, lmax=LMAX)

    return [Clsp, Clsm]  

def get_A0(Clrealization, Clinput):
    NCls1=len(Clrealization)
    NCls2=len(Clinput)
    if (NCls1!=NCls2):
        print "different number of Cl values in the two spectra. Please check!"
    
    # also do a weighted average 
    total=np.sum([2*i+1 for i in range(1, NCls1)])
    A0sum=np.sum([(2*i+1)*(Clrealization[i]/Clinput[i]-1) for i in range(1, NCls1)])
    return A0sum/total

def weightedA(Clp, Cln):
    ls=range(len(Clp))
    total=np.sum([2*i+1 for i in ls])
    num=np.sum([(2*i+1)*((Clp[i]-Cln[i])/(Clp[i]+Cln[i])) for i in ls])
    return num/total

def Ais(map1, LMAX):
    """
    get the Cl power asymmetry in a particular direction 
    """
    dir1=np.array([1.,0.,0.])
    dir2=np.array([0.,1.,0.])
    dir3=np.array([0.,0.,1.])
    dirs=[dir1, dir2, dir3]
    A123=[]
    for direction in dirs:
        Clsp, Clsm = get_hem_Cls(map1, direction, LMAX)
        A123.append(weightedA(Clsp[2:],Clsm[2:]))
        
    return A123

def AistoA(Ais):
    return np.sqrt(Ais[0]**2.0+Ais[1]**2.0+Ais[2]**2.0)

def get_dipole(map1):
    uw=np.array([1.]*len(map1))   
    try:
        ddir=hp.remove_dipole(map1, fitval=True)[2] #direction vector array)
    except:
        ddir=hp.remove_dipole(map1, uw, fitval=True)[2]
    
    return ddir