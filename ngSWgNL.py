"""
"""
import numpy as np
import healpy as hp
from cmbutils import get_hem_Cls, Ais, AistoA, get_dipole, get_A0
from scipy.special import gamma, gammaln    #gammaln useful for large l, gamma(l) is too large

class gNLSWMap(object):
    """
    """
    def __init__(self, gnls=[10000], LMAX=100, NSIDE=128, N=50, readGmap=-1, mapsdir="maps/"):
        """
        * LMAX   is the maximum l to be used to compute power asymmetry and other quantities whereas
                 the maps are generated using LMAX*3 multipoles
        * N      is the number of extra efolds assumed for C_0 calculation, if N=0, C0=0.0
        """  
        self.gausmap=None
        self.gnls=gnls
        self.Ngnls=len(gnls)
        self.gnlmaps=None
        self.lmax=LMAX  # this is the lmax to compute Ais and As; for mapmaking use 3*lmax
        self.NSIDE=NSIDE
        self.efolds=N
        self.mapsdir=mapsdir
        self.get_SWCls()
        self.generate_gaus_map(readGmap)
    
    def generate_gaus_map(self, readGmap=-1):
        if (readGmap<0):
            self.gausalm=hp.synalm(self.inputCls)
        else:
            self.gausalm=hp.read_alm(self.mapsdir+"gmap_"+str(readGmap)+".fits")
                        
        self.gausmap=hp.alm2map(self.gausalm, nside=self.NSIDE) # includes the monopole bit
    
    def save_gaus_map(self, mapN):
        print "for the gNL we will use the maps generated earlier, so no need to save"
        #print "writing "+str(mapN)
        #hp.write_alm(self.mapsdir+"gmap_"+str(mapN)+".fits", self.gausalm0)
    
    def generate_gnl_maps(self, phisq=0.):
        self.gnlmaps=[]
        
        cubemap=np.power(self.gausmap, 3.0)
        
        for i in range(self.Ngnls):
            self.gnlmaps.append(self.gausmap+9.*self.gnls[i]*(cubemap-3.0*phisq*self.gausmap))
                    
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
            
    def calculate(self):
        """
        calculate various properties of the gaussian, fnl, gnl maps i.e. compute
        * Ais
        * As
        * Cls
        * dipoles using remove_dipole
        """
        inpCls=self.inputCls[:self.lmax+1]
        
        if (self.gausmap!=None):
            self.gausCls=hp.alm2cl(self.gausalm)[0:self.lmax+1]
            self.gausA0=get_A0(self.gausCls[1:], inpCls[1:])
            self.gausAi=Ais(self.gausmap, self.lmax); self.gausA=AistoA(self.gausAi)
            self.gausmp=hp.remove_monopole(self.gausmap, fitval=True)[1]
            self.gausdipole=get_dipole(self.gausmap)
            
        if (self.gnlmaps!=None):
            self.gnlCls=[]; self.gnlA0=[]; self.gnlAi=[]
            self.gnlmp=[]; self.gnldipole=[]; self.gnlA=[]
            
            for i in range(self.Ngnls):
                self.gnlCls.append(hp.anafast(self.gnlmaps[i], nspec=self.lmax))
                self.gnlA0.append(get_A0(self.gnlCls[i][1:], inpCls[1:]))
                self.gnlAi.append(Ais(self.gnlmaps[i], self.lmax)); self.gnlA.append(AistoA(self.gnlAi[i]))
                self.gnlmp.append(hp.remove_monopole(self.gnlmaps[i], fitval=True)[1])
                self.gnldipole.append(get_dipole(self.gnlmaps[i]))
        