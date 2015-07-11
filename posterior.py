"""
posterior class to plot fNL posterior and other relevant quantities
using A, A0 as the data (with additional cosmological inputs like Aphi, ns, N)
"""
import numpy as np
from scipy import integrate
from scipy.stats import chi, norm

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.rcParams.update({'ytick.major.pad': 9})
matplotlib.rcParams.update({'xtick.major.pad': 7})

class PosteriorfNL(object):
    """
    """
    def __init__(self, ns=0.965, N=50):
        self.ns=ns
        self.efolds=N
        self.A0const=7.94E-10
        self.A1const=0.0258/500. # for ns=0.965
        self.sigmaG=0.0137876
        self.pfNL0=self.pvalue_norm(0.0, mu=40, sigma=20)
        
    def pdf(self, A, fNL):
        """
        for power asymmetry amplitude A as the variable, and fNL as the parameter
        the pdf is a Maxwell distribution
        """
        return chi.pdf(A, 3, scale=np.sqrt((self.A1const*fNL)**2.0+self.sigmaG**2.0))
    
    def pvalue(self, A, fNL):
        """
        p-value for A-observed=A given fNL=fNL is the non-Gaussian universe
        """
        fnc = lambda A: self.pdf(A, fNL)
        pvalue = integrate.quad(fnc, A, np.infty)[0]
        return pvalue
    
    def pvalue_norm(self, x, mu=0, sigma=1):
        """
        return the p-value of a normal distribution N(mu, sigma)
        """
        fnc = lambda x: norm.pdf(x, loc=mu, scale=sigma)
        if (x>=mu):
            return integrate.quad(fnc, x, np.infty)[0]
        else:
            return integrate.quad(fnc, -np.infty, x)[0]
    
    def sddr(self, A):
        return self.pvalue(A, 0)/self.pfNL0
    
    def posterior_norm(self, A):
        """
        get the normalization for the posteriro pdf
        """
        fnc = lambda fNL: self.pdf(A, fNL)
        norm = integrate.quad(fnc, 0.0, np.infty)[0]
        return norm
    
    def posterior(self, fNLlist, A):
        """
        given an array of fNL values return the posterior likelihood pdf after normalizing
        """
        # avoid too small values of A
        if (A<0.000001):
            A=0.000001
            
        norm = self.posterior_norm(A)
        print norm
        return np.array([self.pdf(A, fNL)/norm for fNL in fNLlist])
    
    def plot_posterior(self, A, ls="solid", lw=2, clr="b", fNLmin=0, fNLmax=2000, ymin=1E-6, ALP=0.02):
        """
        plot a posterior pdf curve
        """
        fNLlist = np.arange(fNLmin, fNLmax, 1)
        plt.plot(fNLlist, self.posterior(fNLlist, A), linestyle=ls, linewidth=lw, label=r"$A="+str(A)+"$", color=clr)
        
    def plot_posteriors(self, Alist=[0.0, 0.04, 0.055, 0.06], 
                              llist=["dashed", "dashdot", "dotted", "solid"],
                              clist=["black", "red", "blue", "green"],
                              LW=2, fNLmax=2000, ymin=5E-6, ymax=5E-6, ALP=0.02):
        """
        plot multiple posterior curves
        """
        for i in range(0, len(Alist)):
            self.plot_posterior(Alist[i], ls=llist[i], lw=LW, clr=clist[i], fNLmax=fNLmax, ymin=ymin, ALP=ALP)
        
        if (ymin!=ymax):
            plt.ylim(ymin, ymax)
            
        plt.xlabel(r"$f_{\rm NL}$")
        plt.ylabel(r"$p(f_{\rm NL})$")
        plt.legend()
        plt.yscale('log')
    
    def plot_pvalue(self, A, ls="solid", lw=2, clr="b", fNLmin=0, fNLmax=500):
        fNLlist=np.arange(fNLmin, fNLmax, 10)
        pvalues=np.array([self.pvalue(A, fNL) for fNL in fNLlist])
        plt.plot(fNLlist, pvalues, linestyle=ls, linewidth=lw, label=r"$A="+str(A)+"$", color=clr)
        
    def plot_pvalues(self, Alist=[0.0, 0.040, 0.055, 0.060],
                           llist=["dashed", "dashdot", "dotted", "solid"],
                           clist=["black", "red", "blue", "green"],
                           LW=2, fNLmax=500):
        for i in range(0, len(Alist)):
            self.plot_pvalue(Alist[i], ls=llist[i], lw=LW, clr=clist[i], fNLmax=fNLmax)
        
        plt.axhline(0.005, ls="dashed", lw=1.5, color="black", alpha=0.35)
        plt.text(220, 0.0055, r"$2.8\sigma$")
        plt.axhline(0.001, ls="dashed", lw=1.5, color="black", alpha=0.35)
        plt.text(220, 0.0011, r"$3.3\sigma$")
        plt.xlabel(r"$f_{\rm NL}$")
        plt.ylabel(r"$p{\rm-value}$")
        plt.legend(loc=0)
        plt.ylim(1E-4, 1.1)
        plt.yscale('log')
