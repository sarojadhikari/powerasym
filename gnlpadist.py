"""
padists.py [PowerAsymmetryDistribution] for gNL
"""

from padists import PowerAsymmetryDistribution, NtoSTR, plot_hist
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
            
        self.phisqs=np.load(self.basedir+"phisq.npy")
        self.phisq=np.mean(self.phisqs)
        
    def plot_A0(self):
        if (self.TYPE!='gNL'):
            plot_hist(plt, self.gA0, clr=self.clrs[0], alp=0.5, ht='step')
        
        sG=np.sqrt(np.var(self.gA0))
        mG=np.mean(self.gA0)
        amin, amax=plt.xlim()
        alist=np.arange(3*amin, amax*3, amax/250.)
        theoryGdist=norm.pdf(alist, loc=mG, scale=sG)
        
        for i in range(len(self.fgnls)):
            if (self.theoryplot==False):
                lbl=self.TYPELABEL+"$=$"+NtoSTR(self.fgnls[i])
            else:
                lbl=None
            
            if (self.TYPE=='gNL'):
                #plot_hist(plt, self.fgNLA0[i]-self.gA0+6*self.fgnls[i]*self.efolds*self.A0const-6*self.fgnls[i]*self.phisq, clr=self.clrs[i+1], alp=ALPHA, labl=lbl, ht='stepfilled')
                neg=np.min(self.fgNLA0[i]-self.gA0)
                const=9*6.*self.fgnls[i]*self.phisq
                #scl=6.0*self.A0const*self.Nconst*self.fgnls[i]
                plot_hist(plt, (self.fgNLA0[i]-self.gA0+const)/(1-9.*3.*self.fgnls[i]*1.5E-8), clr=self.clrs[i+1], alp=ALPHA, labl=lbl, ht='stepfilled')
                
            amin, amax=plt.xlim()
            alist=np.arange(3*amin, amax*3, 0.001)
            if (self.theoryplot):
                if (self.TYPE=='fNL'):
                    theorynGdist=norm.pdf(alist, loc=mG, scale=np.sqrt(16.*self.A0const*self.Nconst*self.fgnls[i]**2.0+sG**2.0))
                if (self.TYPE=='gNL'):
                    scl=6.0*self.A0const*self.Nconst*self.fgnls[i]
                    #scl=2.0*self.A0const*self.efolds*self.fgnls[i]
                    theorynGdist=np.sign(self.fgnls[i])*chi2.pdf(alist,1, scale=scl)
                    # new
                    """
                    sigmaphi00=np.sqrt(self.A0const*self.efolds)
                    a0=24.*np.pi*self.fgnls[i]
                    theorynGdist=1./a0/np.sqrt(2.*np.pi)/sigmaphi00 * np.sqrt(a0/(alist+a0*sigmaphi00^2))*np.exp(-(alist+a0*sigmaphi00**2.0)/(2.*a0*sigmaphi00**2.0))
                    self.theorynGdist=theorynGdist
                    """
                    
                plt.plot(alist, theorynGdist, self.clrs[i+1], linestyle=self.ls[i+1], linewidth=LW, label=self.TYPELABEL+"="+NtoSTR(self.fgnls[i]))
            
        plt.xlabel(r'$A_0$')
        plt.ylabel(r'$p(A_0)$')
        plt.yscale('log')
        if (self.TYPE=='fNL'):
            plt.xlim(-1.0, 1.0)
            plt.ylim(0.1, 100)
        if (self.TYPE=='gNL'):
            #plt.xscale('log')
            plt.ylim(0.05, 1000)
            plt.xlim(-0.02, 0.25)
        
        plt.legend()
        
        
    def plot_Ai(self, histtype="stepfilled"):
        if histtype=="step":
            ALPHAG=0.5
        else:
            ALPHAG=ALPHA
        
        sG=np.sqrt(np.var(self.gAi.flatten()))
        print sG
        mG=np.mean(self.gAi.flatten())
        amin, amax = np.min(self.gAi), np.max(self.gAi)
        alist=np.arange(2*amin, amax*2, amax/250.)
        theoryGdist=norm.pdf(alist, loc=mG, scale=sG)
        if (self.TYPE!='gNL'):
            plt.plot(alist, theoryGdist,self.clrs[0], linestyle=self.ls[0], linewidth=LW, label=self.TYPELABEL+"=0")
            plot_hist(plt, self.gAi.flatten(), clr=self.clrs[0], alp=ALPHAG, ht=histtype)

        for i in range(len(self.fgnls)):
            if (self.theoryplot==False):
                lbl=self.TYPELABEL+"="+NtoSTR(self.fgnls[i])
            else:
                lbl=None
                
            plot_hist(plt, self.fgNLAi[i].flatten(), clr=self.clrs[i+1], alp=ALPHA, labl=lbl)
            
            if (self.theoryplot):
                if (self.TYPE=='fNL'):
                    theorynGdist=norm.pdf(alist, loc=mG, scale=np.sqrt((self.A1const*self.fgnls[i])**2.0+sG**2.0))
                elif (self.TYPE=='gNL'):
                    sigmax0=np.sqrt(self.A0const*self.Nconst)
                    sigmax1=8.*self.fgnls[i]*np.sqrt(3.*np.pi*self.A1const)
                    theorynGdist=kn(0, np.abs(alist)/sigmax0/sigmax1)/np.pi/sigmax0/sigmax1
                    
                plt.plot(alist, theorynGdist, self.clrs[i+1], linestyle=self.ls[i+1], linewidth=LW, label=self.TYPELABEL+"="+NtoSTR(self.fgnls[i]))
        
        plt.xlim(-0.15, 0.15)
        plt.xlabel(r'$A_i$')
        plt.ylabel(r'$p(A_i)$')
        plt.legend(fontsize=20)
        
    def plot_A(self, histtype="stepfilled"):
        if histtype=="step":
            ALPHAG=0.5
        else:
            ALPHAG=ALPHA
            
        plot_hist(plt, self.gA, clr=self.clrs[0], alp=ALPHAG, ht=histtype)

        sG=np.sqrt(np.var(self.gAi))
        mG=np.mean(self.gAi)
        amin, amax = plt.xlim()
        alist=np.arange(0.0, np.max(self.fgNLA), 0.001)
        theoryGdist=chi.pdf(alist, 3,scale=sG)
        plt.plot(alist, theoryGdist,self.clrs[0], linestyle=self.ls[0], linewidth=LW, label=self.TYPELABEL+r"$=0$")

        for i in range(len(self.fgnls)):
            if (self.theoryplot==False):
                lbl=self.TYPELABEL+"="+NtoSTR(self.fgnls[i])
            else:
                lbl=None
                
            plot_hist(plt, self.fgNLA[i], clr=self.clrs[i+1], alp=ALPHA, labl=lbl)
            
            if (self.theoryplot):
                theorynGdist=chi.pdf(alist, 3, scale=np.sqrt((self.A1const*self.fgnls[i])**2.0+sG**2.0))
                plt.plot(alist, theorynGdist, self.clrs[i+1], linestyle=self.ls[i+1], linewidth=LW, label=self.TYPELABEL+"="+NtoSTR(self.fgnls[i]))
        
        # theory plots
        
        #plt.xlim(0.0, 0.2)
        #plt.xticks([0.02, 0.04, 0.06, 0.08, 0.10])
        plt.xlabel(r'$A$')
        plt.ylabel(r'$p(A)$')
        plt.legend(fontsize=20)