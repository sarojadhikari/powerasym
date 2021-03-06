"""
Created on Fri Apr 10 14:54:45 2015

This file defines the SachsWolfeMap class to generate Gaussian and local-type
non-Gaussian CMB maps, to compute the power spectrum modulation (both
monopole and dipole).

Since we are subsampling a larger Universe to CMB skies of the size of our
observable universe, we need to keep the C_0 and C_1 terms, as we want to
study possible monopole and dipole modulations.

In the local non-Gaussian model, since the L=0 and L=1 modulations are independent,
we can include only the C_1 term to study the dipole modulations.
@author: sarojadhikari
"""
import numpy as np
import healpy as hp
from .cmbutils import get_hem_Cls, Ais, AistoA, get_dipole, get_A0
from scipy.special import gamma, gammaln    #gammaln useful for large l, gamma(l) is too large

class SachsWolfeMap(object):
    """
    Defines SachsWolfe CMB maps. First computes the Gaussian Sachs Wolfe map
    corresponding to the input primordial power spectrum (given A and ns). Then,
    compute the corresponding non-Gaussian maps for the fNL values specified.

    In the Sachs-Wolfe approximation, one can simply square the temperature field
    and add the appropriately weighted temperature squared field to the temperature
    field to obtain the CMB temperature map for the primordial local-type
    non-Gaussianity.
    """
    def __init__(self, fnls=[100], LMAX=100, NSIDE=128, N=50, readGmap=-1, nodipole=False, mapsdir="maps/"):
        """
        * LMAX   is the maximum l to be used to compute power asymmetry and other quantities whereas
                 the maps are generated using LMAX*3 multipoles
        * N      is the number of extra efolds assumed for C_0 calculation, if N=0, C0=0.0
        """
        self.gausmap=None
        self.fnls=fnls
        self.Nfnls=len(fnls)
        self.fnlmaps0=None
        self.lmax=LMAX  # this is the lmax to compute Ais and As; for mapmaking use 3*lmax
        self.NSIDE=NSIDE
        self.efolds=N
        self.nodipole=nodipole
        self.mapsdir=mapsdir
        self.get_SWCls()
        self.generate_gaus_map(readGmap)

    def generate_gaus_map(self, readGmap=-1):
        if (readGmap<0):
            self.gausalm0=hp.synalm(self.inputCls)
        else:
            self.gausalm0=hp.read_alm(self.mapsdir+"gmap_"+str(readGmap)+".fits")

        self.gausalm1=np.copy(self.gausalm0); self.gausalm1[0]=0.0
        if (self.nodipole):
            ndxfl=np.ones(len(self.inputCls))
            ndxfl[1]=0.0
            hp.almxfl(self.gausalm1, ndxfl, inplace=True)
            hp.almxfl(self.gausalm0, ndxfl, inplace=True)

        self.gausmap0=hp.alm2map(self.gausalm0, nside=self.NSIDE) # includes the monopole bit
        self.gausmap1=hp.alm2map(self.gausalm1, nside=self.NSIDE) # does not include the monopole

    def save_gaus_map(self, mapN):
        print ("writing "+str(mapN))
        hp.write_alm(self.mapsdir+"gmap_"+str(mapN)+".fits", self.gausalm0)
        #hp.write_map(self.mapsdir+"gmap_"+str(mapN)+".fits", self.gausmap)

    def generate_fnl_maps(self, phisq0=0., phisq1=0.):
        #self.fnlmaps0=[]
        self.fnlmaps1=[]
        
        #sqmap0=np.power(self.gausmap0, 2.0)
        sqmap1=np.power(self.gausmap1, 2.0)

        for i in range(self.Nfnls):
            #self.fnlmaps0.append(self.gausmap0+3.*self.fnls[i]*(sqmap0-phisq0))
            self.fnlmaps1.append(self.gausmap1+3.*self.fnls[i]*(sqmap1-phisq1))
            
        self.fnlalm1 = hp.map2alm(self.fnlmaps1[0])

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
            k0rcmbfactor=1.2135 #((pi./(rcmb*k0))**(ns-1)
            nfactor=k0rcmbfactor*4.*np.pi*np.pi*np.power(2.0, ns-4.0)/9.0
            self.inputCls=np.array([Aphi*nfactor*np.exp(gammaln(l+ns/2.0-0.5)-gammaln(l+5./2-ns/2.0))*gamma(3.0-ns)/np.power(gamma(2.0-ns/2.0), 2.0) for l in range(1, self.lmax*3)])
            C0factor=k0rcmbfactor*4.*np.pi*Aphi/9.0
            C0=C0factor*(1.0-np.exp(-(ns-1)*self.efolds))/(ns-1)
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
        inpCls=self.inputCls[:self.lmax+1]

        if (self.gausmap0!=None):
            self.gausCls0=hp.alm2cl(self.gausalm0)[0:self.lmax+1]
            self.gausCls1=hp.alm2cl(self.gausalm1)[0:self.lmax+1]
            self.gausA0=get_A0(self.gausCls0[1:], inpCls[1:])
            self.gausAi=Ais(self.gausmap1, self.lmax); self.gausA=AistoA(self.gausAi)
            self.gausAi2=Ais(self.gausmap0, self.lmax); self.gausA2=AistoA(self.gausAi2)
            self.gausmp=hp.remove_monopole(self.gausmap0, fitval=True)[1]
            self.gausdipole=get_dipole(self.gausmap1)

        if (self.fnlmaps1!=None):
            self.fnlCls0=[]; self.fnlCls1=[]; self.fnlA0=[]; self.fnlAi=[]; self.fnlAi2=[]
            self.fnlmp=[]; self.fnldipole=[]; self.fnlA=[]; self.fnlA2=[]

            for i in range(self.Nfnls):
                #self.fnlCls0.append(hp.anafast(self.fnlmaps0[i], nspec=self.lmax))
                self.fnlCls1.append(hp.anafast(self.fnlmaps1[i], nspec=self.lmax))
                #self.fnlA0.append(get_A0(self.fnlCls0[i][1:], inpCls[1:]))
                self.fnlAi.append(Ais(self.fnlmaps1[i], self.lmax)); self.fnlA.append(AistoA(self.fnlAi[i]))
                #self.fnlAi2.append(Ais(self.fnlmaps0[i], self.lmax)); self.fnlA2.append(AistoA(self.fnlAi2[i]))
                #self.fnlmp.append(hp.remove_monopole(self.fnlmaps0[i], fitval=True)[1])
                self.fnldipole.append(get_dipole(self.fnlmaps1[i]))
