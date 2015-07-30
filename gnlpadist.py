"""
padists.py [PowerAsymmetryDistribution] for gNL
"""

from padists import PowerAsymmetryDistribution
import numpy as np
from scipy.special import kn

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.rcParams.update({'ytick.major.pad': 9})
matplotlib.rcParams.update({'xtick.major.pad': 7})


class gNLPowerAsymmetryDist(PowerAsymmetryDistribution):
    def __init__(self, typ='gNL', datafolder='gnldata', NMAPS=10000, efolds=50, ns=0.965, lmax=100, theory=True):
        self.TYPE=typ
        self.efolds=efolds
        self.ns=ns
        self.basedir=datafolder+str(efolds)+"/"
        self.fgnls=None
        self.nmaps=NMAPS
        self.lmax=lmax
        self.read_data()
        self.clrs=['b', 'g', 'r', 'black']
        self.ls=[':', '--', '-', 'o']
        self.theoryplot=theory
        self.A0const=7.94E-10
        if (self.ns==1.0):
            self.Nconst=self.efolds
        else:
            self.Nconst=1.2135*(1.0-np.exp(-(self.ns-1)*self.efolds))/(self.ns-1)
            
    def read_data(self):
        self.read_gaus_data()
        self.read_ngaus_data()
        if self.TYPE=='gNL':
            self.TYPELABEL=r'$g_{\rm NL}$'
        