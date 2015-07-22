"""
compute the power modulation amplitues in local and scale-dependent local models
"""
import numpy as np

from scipy import integrate
from scipy.special import sph_jn

class PowerModulation(object):
    def __init__(self, ns=0.965, A=7.94E-10, Lbox=14800., nfNL=0.0, k0=0.05, klmin=1E-40):
        self.ns=ns
        self.A=A
        self.Lbox=Lbox
        self.klmin=klmin
        self.nfNL=nfNL
        self.k0=k0
    
    def ASq(self, L):
        integrand = lambda kl: np.power(sph_jn(L, self.Lbox*kl)[0][-1], 2.0)*np.power(kl, self.ns-2.0)
        const = 64.*np.pi*self.A/(np.power(self.k0, self.ns-1.0))
        klmax = np.pi/self.Lbox
        intresult = integrate.quad(integrand, self.klmin, klmax)
        self.asq0=const*intresult[0]
        self.L=L
        return self.asq0
    
    def ASqk(self, k, L):
        # first get asq0
        self.ASq(L)
        
        
        