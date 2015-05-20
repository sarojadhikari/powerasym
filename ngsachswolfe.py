# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:54:45 2015

local non-gaussianity and modulations. simply assume that:
l(l+1)/2pi C_l = A^2/25 = constant
A^2(\zeta) =2.2 x 10^{-9}

we will always keep the cosmic variance primordial dipole, as we want to study
statistics related to it.

keeping or ignoring the monopole is a choice at this stage as local non-g
only shifts the monopole; we will make it 0 for the local case

modifications: May 17
=====================
save 10000 base Gaussian maps for the time being, so that we can reproduce the exact
maps and statistics later if necessary; and generate statistics for the following:
* Gaussian case
* local non-Gaussian with fNL=100, 250, 500
* local non-Gaussian without background dipole i.e. C_1 = 0
* local non-Gaussian with gNL=1E6, 2.5E6, 5E6
* in each case, save the C_ls of each map upto C_l=100


@author: sarojadhikari
"""
import numpy as np
import healpy as hp
from cmbutils import get_hem_Cls, Ais, AistoA, get_dipole

class SachsWolfeMap(object):
    """
    defines a healpix cmb map, given the temperature anisotropy values
    for each pixel, assumed Gaussian. The pixel values are calculated using 
    another class/functions beforehand (basically healpix synfast function). 
    This class can then also generate simple local type non-Gaussian maps 
    (fNL, gNL etc) for the Sachs Wolfe case.
    """
    def __init__(self, fnl=100, gnl=1.0E6, LMAX=100, NSIDE=64, readGmap=-1, nodipole=False):
        self.gausmap=None
        self.fnl=fnl
        self.gnl=gnl
        self.fnlmap=None
        self.gnlmap=None
        self.lmax=LMAX  # this is the lmax to compute Ais and As; for mapmaking use 3*lmax
        self.NSIDE=NSIDE
        self.nodipole=nodipole
        self.get_SWCls()
        #self.inputCls=np.array([self.Cl(i) for i in range(LMAX)])
        
        if (readGmap>-1):
            self.gausmap=hp.read_map("maps/gmap_"+str(readGmap)+".fits")
        else:
            self.generate_gaus_map()
    
    def generate_gaus_map(self):
        self.gausmap=hp.synfast(self.inputCls, nside=self.NSIDE)
    
    def save_gaus_map(self, mapN):
        print "writing "+str(mapN)
        hp.write_map("maps/gmap_"+str(mapN)+".fits", self.gausmap)
    
    def generate_fnl_map(self):
        self.fnlmap=hp.remove_monopole(self.gausmap+3.*self.fnl*self.gausmap*self.gausmap)
        
    def generate_gnl_map(self):
        self.gnlmap=hp.remove_monopole(self.gausmap+9.*self.gnl*np.power(self.gausmap, 3.0))

    def get_SWCls(self, Asq=2.2e-9):
        """
        return the Cl value at a multipole assuming scale independent power spectrum
        with amplitude Asq and Saches Wolfe effect only
        Asq=amplitude of fluctuations in curvature perturbation R
        """
        self.inputCls=np.array([Asq*2.0*np.pi/(25.0*l*(l+1)) for l in range(1, self.lmax*3)])
        self.inputCls=np.append(0.0, self.inputCls)
        
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
            self.gausAi=Ais(self.gausmap, self.lmax)
            self.gausA=AistoA(self.gausAi)
            self.gausCls=hp.anafast(self.gausmap, lmax=self.lmax)
            self.gausdipole=get_dipole(self.gausmap)
        
        if (self.fnlmap!=None):
            self.fnlAi=Ais(self.fnlmap, self.lmax)
            self.fnlA=AistoA(self.fnlAi)
            self.fnlCls=hp.anafast(self.fnlmap, lmax=self.lmax)
            self.fnldipole=get_dipole(self.fnlmap)
            
        if (self.gnlmap!=None):
            self.gnlAi=Ais(self.gnlmap, self.lmax)
            self.gnlA=AistoA(self.gnlAi)
            self.gnlCls=hp.anafast(self.gnlmap, lmax=self.lmax)
            self.gnldipole=get_dipole(self.gnlmap)
        