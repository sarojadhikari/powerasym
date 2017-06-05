import healpy as hp
import numpy as np

def get_hem_Cls(skymap, direction, LMAX=256, deg=90.):
    """
    from the given healpix skymap, return Cls for two hemispheres defined by the
    direction given, useful to study the possible scale dependence of power modulation

    direction should be a unit vector
    """
    # generate hemispherical mask
    NPIX=len(skymap)
    NSIDE=hp.npix2nside(NPIX)
    maskp=np.array([0.]*NPIX)
    disc=hp.query_disc(nside=NSIDE, vec=direction, radius=0.0174532925*deg)
    maskp[disc]=1.
    #skymap=hp.remove_monopole(skymap)
    map1=hp.ma(skymap);
    map1.mask=maskp;
    Clsp=hp.anafast(map1, lmax=LMAX);
    if (deg<90.):
        maskm=np.array([0.]*NPIX)
        disc=hp.query_disc(nside=NSIDE, vec=-direction, radius=0.0174532925*deg);
        maskm[disc]=1.
        map1.mask=maskm;
    else:
        map1.mask=np.logical_not(maskp);

    Clsm=hp.anafast(map1, lmax=LMAX);

    return [Clsp, Clsm]

def get_A0(Clrealization, Clinput):
    NCls1=len(Clrealization)
    NCls2=len(Clinput)
    if (NCls1!=NCls2):
        print ("different number of Cl values in the two spectra. Please check!")

    # also do a weighted average
    total=np.sum([2*i+1 for i in range(2, NCls1)])
    A0sum=np.sum([(2*i+1)*(Clrealization[i]/Clinput[i]-1) for i in range(2, NCls1)])
    return A0sum/total

def weightedA(Clp, Cln, al=0):
    ls = range(len(Clp))
    if al==0:
        al=range(len(Clp))
    total=np.sum([2*al[i]+1 for i in ls])
    num=np.sum([(2*al[i]+1)*(Clp[i]-Cln[i])/(Clp[i]+Cln[i]) for i in ls])
    return num/total

def A_wrt_dir(map1, direction, LMAX, deg=90.):
    """
    get the Cl power asymmetry wrt to the direction specified
    """
    if len(direction)==2:
        direction = hp.dir2vec(direction[0], direction[1], lonlat=True)
        
    Clsp, Clsm = get_hem_Cls(map1, direction, LMAX, deg)
    return weightedA(Clsp[2:], Clsm[2:])

def Ais(map1, LMAX, deg=90.):
    """
    get the Cl power asymmetry in a particular direction
    """
    dir1=np.array([1.,0.,0.])
    dir2=np.array([0.,-1.,0.])
    dir3=np.array([0.,0.,1.])
    dirs=[dir1, dir2, dir3]
    A123=[]
    for direction in dirs:
        Clsp, Clsm = get_hem_Cls(map1, direction, LMAX, deg);
        A123.append(weightedA(Clsp[2:], Clsm[2:]))

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

def cmb_parity(Cls, lmax=28):
    n_plus = np.sum(np.array([(1+(-1)**l) for l in range(2, lmax+1)]))
    n_minus = np.sum(np.array([(1-(-1)**l) for l in range(2, lmax+1)]))
    p_plus = np.sum(np.array([(1+(-1)**l)*l*(l+1)/(4.*np.pi)*Cls[l] for l in range(2, lmax+1)]))/n_plus
    p_minus= np.sum(np.array([(1-(-1)**l)*l*(l+1)/(4.*np.pi)*Cls[l] for l in range(2, lmax+1)]))/n_minus
    #print (p_plus, p_minus)
    return p_minus/p_plus
