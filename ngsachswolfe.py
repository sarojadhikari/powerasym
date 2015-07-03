"""
Created on Fri Apr 10 14:54:45 2015

local non-gaussianity and modulations. 
for scale invariant power spectrum with amplitude A (Bardeen potential):
l(l+1)/2pi C_l = A/9 = constant
A(\phi) = (3./5)^2 * 2.2 x 10^{-9}   (Planck 2015)

The expression for non-scale invariant case is slightly more complicated, but
mathematica can do the integration; this code uses the output from the 
mathematica script SachsWolfeCls.nb

Since we are subsampling a larger universe to CMB skies of the size of our
observable universe, we need to keep the C_0 and C_1 terms, as we want to
study possible monopole and dipole modulations.

@author: sarojadhikari
"""
import numpy as np
import healpy as hp
from cmbutils import get_hem_Cls, Ais, AistoA, get_dipole, get_A0
from scipy.special import gamma, gammaln

class SachsWolfeMap(object):
    """
    defines a healpix cmb map, given the temperature anisotropy values
    for each pixel, assumed Gaussian. The pixel values are calculated using 
    another class/functions beforehand (basically healpix synfast function). 
    This class can then also generate simple local type non-Gaussian maps 
    (fNL, gNL etc) for the Sachs Wolfe case.
    """
    def __init__(self, fnl=100, gnl=1.0E5, LMAX=100, NSIDE=64, N=50, readGmap=-1, nodipole=False, mapsdir="maps/"):
        """
        * LMAX   is the maximum l to be used to compute power asymmetry and other quantities whereas
                 the maps are generated using LMAX*3 multipoles
        * N      is the number of extra efolds assumed for C_0 calculation, if N=0, C0=0.0
        """  
        self.gausmap=None
        self.fnl=fnl
        self.gnl=gnl
        self.fnlmap=None
        self.gnlmap=None
        self.lmax=LMAX  # this is the lmax to compute Ais and As; for mapmaking use 3*lmax
        self.NSIDE=NSIDE
        self.efolds=N
        self.nodipole=nodipole
        self.mapsdir=mapsdir
        self.get_SWCls()
        
        if (readGmap>-1):
            self.gausmap=hp.read_map(self.mapsdir+"gmap_"+str(readGmap)+".fits")
            if (nodipole):
                uw=np.array([1.]*len(self.gausmap))
                self.gausmap=hp.remove_dipole(self.gausmap, uw)
        else:
            self.generate_gaus_map()
    
    def generate_gaus_map(self):
        self.gausmap=hp.synfast(self.inputCls, nside=self.NSIDE)
    
    def save_gaus_map(self, mapN):
        print "writing "+str(mapN)
        hp.write_map(self.mapsdir+"gmap_"+str(mapN)+".fits", self.gausmap)
    
    def generate_fnl_map(self, phisq=0.):
        self.fnlmap=self.gausmap+3.*self.fnl*(np.power(self.gausmap, 2.0)-phisq)
        
    def generate_gnl_map(self, phisq=0.):
        self.gnlmap=self.gausmap+9.*self.gnl*(np.power(self.gausmap, 3.0)-3.0*self.gausmap*phisq)

    def get_SWCls(self, Aphi=7.94E-10, ns=0.965):
        """
        return the Cl value at a multipole assuming scale independent power spectrum
        with amplitude Asq and Saches Wolfe effect only
        Asq=amplitude of fluctuations in curvature perturbation R
        """
        if (ns==1):
            self.inputCls=np.array([Aphi*2.0*np.pi/(9.0*l*(l+1)) for l in range(1, self.lmax*3)])
            C0=4*np.pi*Aphi*self.efolds/9.0
        else:
            print ns
            k0rcmbfactor=1.2135 #((pi./(rcmb*k0))**(ns-1)
            nfactor=k0rcmbfactor*4.*np.pi*np.pi*np.power(2.0, ns-4.0)/9.0
            self.inputCls=np.array([Aphi*nfactor*np.exp(gammaln(l+ns/2.0-0.5)-gammaln(l+5./2-ns/2.0))*gamma(3.0-ns)/np.power(gamma(2.0-ns/2.0), 2.0) for l in range(1, self.lmax*3)])
            C0factor=k0rcmbfactor*4.*np.pi*Aphi/9.0
            C0=C0factor*(1.0-np.exp(-(ns-1)*self.efolds))/(ns-1)
            print C0
        self.inputCls=np.append(C0, self.inputCls)
        
        if (self.nodipole):
            self.inputCls[1]=0.0
    
    def calculate(self):
        """
        calculate various properties of the gaussian, fnl, gnl maps i.e. compute
        * Ais
        * As
        * Cls
        * dipoles using remove_dipole
        """
        if (self.gausmap!=None):
            self.gausCls=hp.anafast(self.gausmap, lmax=self.lmax)
            self.gausA0=get_A0(self.gausCls, self.inputCls[:self.lmax+1])
            self.gausAi=Ais(self.gausmap, self.lmax)
            self.gausA=AistoA(self.gausAi)
            self.gausmp=hp.remove_monopole(self.gausmap, fitval=True)[1]
            self.gausdipole=get_dipole(self.gausmap)
            
        if (self.fnlmap!=None):
            self.fnlCls=hp.anafast(self.fnlmap, lmax=self.lmax)
            self.fnlA0=get_A0(self.fnlCls, self.inputCls[:self.lmax+1])
            
            self.fnlAi=Ais(self.fnlmap, self.lmax)
            self.fnlA=AistoA(self.fnlAi)
            self.fnlmp=hp.remove_monopole(self.fnlmap, fitval=True)[1]
            self.fnldipole=get_dipole(self.fnlmap)
        
        if (self.gnlmap!=None):
            self.gnlCls=hp.anafast(self.gnlmap, lmax=self.lmax)
            self.gnlA0=get_A0(self.gnlCls, self.inputCls[:self.lmax+1])
            
            self.gnlAi=Ais(self.gnlmap, self.lmax)
            self.gnlA=AistoA(self.gnlAi)
            self.gnlmp=hp.remove_monopole(self.gnlmap, fitval=True)[1]
            self.gnldipole=get_dipole(self.gnlmap)