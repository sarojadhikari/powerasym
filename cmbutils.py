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
    Clsp=hp.anafast(map1, lmax=LMAX*2.0);
    if (deg<90.):
        maskm=np.array([0.]*NPIX)
        disc=hp.query_disc(nside=NSIDE, vec=-direction, radius=0.0174532925*deg);
        maskm[disc]=1.
        map1.mask=maskm;
    else:
        map1.mask=np.logical_not(maskp);

    Clsm=hp.anafast(map1, lmax=LMAX*2.0);

    return [Clsp[0:LMAX+1], Clsm[0:LMAX+1]]

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

def A_wrt_dir(map1, direction, LMAX, LMIN=2, deg=90.):
    """
    get the Cl power asymmetry wrt to the direction specified
    """
    if len(direction)==2:
        direction = hp.dir2vec(direction[0], direction[1], lonlat=True)
        
    Clsp, Clsm = get_hem_Cls(map1, direction, LMAX, deg)
    return weightedA(Clsp[LMIN:], Clsm[LMIN:])

def Ais(map1, LMAX, LMIN=2, deg=90.):
    """
    get the Cl power asymmetry in a particular direction
    """
    dir1=np.array([1.,0.,0.])
    dir2=np.array([0.,-1.,0.])
    dir3=np.array([0.,0.,1.])
    dirs=[dir1, dir2, dir3]
    A123=[]
    for direction in dirs:
        A123.append(A_wrt_dir(map1, direction, LMAX, LMIN, deg))

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

"""
Dipole modulation power asymmetry estimators and other utilities below
"""

def sigma_X(Cllist, lmin=2, lmax=100):
    """
    the cosmic variance of the modulation amplitude estimator
    from 1704.03143
    """
    sxarray = np.array([(l+1)*4*(Cllist[l]+Cllist[l+1])**2.0/Cllist[l]/Cllist[l+1]
                        for l in range(lmin, lmax)])
    sigmaXsq = 12./(np.sum(sxarray))
    return np.sqrt(sigmaXsq)

def _Alm(l, m):
    num = (l+1)**2.0 - m * m
    den = (2*l+1)*(2*l+3)
    return np.sqrt(num/den)

def _Blm(l, m):
    num = (l+m+1)*(l+m+2)
    den = 2*(2*l+1)*(2*l+3)
    return np.sqrt(num/den)

def DeltaXM(alms, Cllist, lmin=2, lmax=100):
    """
    given the alms, and theory Cls, estimate the 
    dipole modulation parameters DeltaX0, DeltaX(+-1) and return them
    
    assumptions: full-sky, noiseless CMB, small modulation (linear order)
    See: Planck 2015 isotropy paper (1506.07135, Appendix C) or 1512.02618
    """
    num0 = 0.; num1 = 0.;
    den = 0.
    alm_lmax = hp.Alm.getlmax(len(alms))
    
    for l in range(lmin, lmax):
        Cll1inv = 1./(Cllist[l]*Cllist[l+1])
        dCl = (Cllist[l] + Cllist[l+1])*2
        den = den + dCl**2.0*(l+1)*Cll1inv
        for m in range(0, l+1):
            alm = alms[hp.Alm.getidx(lmax=alm_lmax, l=l, m=m)]
            alplus1m = alms[hp.Alm.getidx(lmax=alm_lmax, l=l+1, m=m)]
            al1m1 = alms[hp.Alm.getidx(lmax=alm_lmax, l=l+1, m=m+1)]
            common_factor = 6. * dCl * alm.conjugate() * Cll1inv
            num0 = num0 + common_factor * _Alm(l,m) * alplus1m
            num1 = num1 + common_factor * _Blm(l,m) * al1m1
            
        for m in range(-l, 0):
            alm = (-1)**(-m)*alms[hp.Alm.getidx(lmax=alm_lmax, l=l, m=-m)].conjugate()
            alplus1m = (-1)**(-m)*alms[hp.Alm.getidx(lmax=alm_lmax, l=l+1, m=-m)].conjugate()
            al1m1 = (-1)**(-m+1)*alms[hp.Alm.getidx(lmax=alm_lmax, l=l+1, m=-m+1)]
            common_factor = 6. * dCl * alm.conjugate() * Cll1inv
            num0 = num0 + common_factor * _Alm(l,m) * alplus1m
            num1 = num1 + common_factor * _Blm(l,m) * al1m1
            
    return num0/den, num1/den
            
def Athetaphi(alms, Cllist, lmin=2, lmax=100):
    """
    use this function to obtain the
        * Amplitude A of dipole modulation
        * (theta, phi) direction of the modulation
        
    the inputs are
        * alms - the alms from the input CMB (T/E etc) map
        * Cllist - the theory Cls -- Cl(TT) Cl(EE) etc
        * (lmin, lmax) the multipole range to use for calculating
          (A, theta, phi)
    """
    
    DX0, DX1 = DeltaXM(alms, Cllist, lmin, lmax)
    A = np.sqrt(DX0.conjugate()*DX0+2.*DX1.conjugate()*DX1).real
    theta = np.arccos(DX0.real/A)
    phi = np.arctan2(DX1.imag, DX1.real)
    
    return A, theta, phi
