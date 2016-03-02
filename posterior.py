"""
posterior class to plot fNL posterior and other relevant quantities
using A, A0 as the data (with additional cosmological inputs like Aphi, ns, N)
"""
import numpy as np
from scipy import integrate
from scipy.stats import chi, norm, foldnorm

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams.update({'figure.autolayout': False})
matplotlib.rcParams.update({'ytick.major.pad': 9})
matplotlib.rcParams.update({'xtick.major.pad': 7})
lf=18

class PosteriorfNL(object):
    """
    """
    def __init__(self, ns=0.965, N=50):
        self.ns=ns
        self.efolds=N
        self.A0const=7.94E-10
        self.A1const=0.0258/500. # for ns=0.965
        self.sigmaG=0.0137876
        self.sigma0G=0.0142527
        self.Nconst=1.2135*(1.0-np.exp(-(self.ns-1)*self.efolds))/(self.ns-1)

    def pdf(self, A, fNL):
        """
        for power asymmetry amplitude A as the variable, and fNL as the parameter
        the pdf is a Maxwell distribution
        """
        if (fNL>=0.0):
            return chi.pdf(A, 3, scale=np.sqrt((self.A1const*fNL)**2.0+self.sigmaG**2.0))
        else:
            return 0.0

    def pdfA0(self, A0, fNL):
        self.Nconst=1.2135*(1.0-np.exp(-(self.ns-1)*self.efolds))/(self.ns-1)
        return self.pdf_fold(A0, 0.0, np.sqrt(16.*self.Nconst*self.A0const*fNL**2.0+self.sigma0G**2.0))
        #return norm.pdf(A0, loc=0.0, scale=np.sqrt(a0var**2.0+self.sigma0G**2.0))

    def post_A01(self, fNLlist, A, A0):
        post = lambda fNL: self.pdfA0(A0, fNL)*self.pdf(A, fNL)
        pnorm = integrate.quad(post, 0.0, np.infty)[0]
        return np.array([post(fNL)/pnorm for fNL in fNLlist])

    def pvalue(self, A, fNL=0.0):
        """
        p-value for A-observed=A given fNL=fNL is the non-Gaussian universe
        """
        fnc = lambda A: self.pdf(A, fNL)
        pvalue = integrate.quad(fnc, A, np.infty)[0]
        return pvalue

    def pvalue_fold(self, x, mu=-100, sigma=100):
        pdf = lambda x: self.pdf_fold(x, mu, sigma)
        return integrate.quad(pdf, x, np.infty)[0]

    def pvalue_norm(self, x, mu=0, sigma=1):
        """
        return the p-value of a normal distribution N(mu, sigma)
        """
        fnc = lambda x: norm.pdf(x, loc=mu, scale=sigma)
        if (x>=mu):
            return integrate.quad(fnc, x, np.infty)[0]
        else:
            return integrate.quad(fnc, -np.infty, x)[0]

    def sddr(self, A, fNL=-100, sigmafNL=100):
        pfold=self.fNL_bispectrum_pdf([0.0], mean=fNL, sigma=sigmafNL)[0]
        return self.combined_posterior([0.0], A, mean=fNL, sigma=sigmafNL)[0]/pfold

    def sddr_withA0(self, A, A0, N, fNL=-100, sigmafNL=100):
        pfold=self.fNL_bispectrum_pdf([0.0], mean=fNL, sigma=sigmafNL)[0]
        self.efolds=N
        return self.combined_posteriors_withA0([0.0], A, A0, mean=fNL, sigma=sigmafNL)[0]/pfold

    def posterior_norm(self, A):
        """
        get the normalization for the posteriro pdf
        """
        fnc = lambda fNL: self.pdf(A, fNL)
        norm = integrate.quad(fnc, 0.0, np.infty)
        return norm[0]

    def posterior(self, fNLlist, A):
        """
        given an array of fNL values return the posterior likelihood pdf after normalizing
        """
        # avoid too small values of A
        if (A<0.000001):
            A=0.000001

        norm = self.posterior_norm(A)
        #print norm
        return np.array([self.pdf(A, fNL)/norm for fNL in fNLlist])

    def pdf_fold(self, x, mean=-100, sigma=100):
        if x<0.0:
            return 0.0
        else:
            scl=1.0/np.sqrt(2.*np.pi)/sigma
            return scl*(np.exp(-(x-mean)**2.0/2.0/sigma**2.0) + np.exp(-(x+mean)**2.0/2.0/sigma**2.0))

    def fNL_bispectrum_pdf(self, fNLlist, mean=-100, sigma=100):
        """
        return a half gaussian pdf with sigma as the pdf for the magnitude of fNL from
        bispectrum measurements
        """
        return np.array([self.pdf_fold(fNL, mean=mean, sigma=sigma) for fNL in fNLlist])

    def combined_posteriors_withA0(self, fNLlist, A, A0=0.3, mean=-100, sigma=100):
        if (A<0.000001):
            A=0.000001
        #scl = 1.0/np.sqrt(2.*np.pi)/sigma
        if (A0!=0.0):
            post = lambda fNL: self.pdf_fold(fNL, mean, sigma)*self.pdf(A, fNL)*self.pdfA0(A0, fNL)
        else:
            post = lambda fNL: self.pdf_fold(fNL, mean, sigma)*self.pdf(A, fNL)

        norm = integrate.quad(post, 0.0, 5000.)[0]
        return np.array([post(fNL)/norm for fNL in fNLlist])

    def combined_posterior(self, fNLlist, A, mean=-100, sigma=100):
        if (A<0.000001):
            A=0.000001
        #scl = 1.0/np.sqrt(2.*np.pi)/sigma
        post = lambda fNL: self.pdf_fold(fNL, mean, sigma)*self.pdf(A, fNL)
        norm = integrate.quad(post, 0.0, np.infty)[0]
        return np.array([post(fNL)/norm for fNL in fNLlist])

    def plot_posterior(self, A, ls="solid", lw=2, clr="b", fNLmin=0, fNLmax=2000, ymin=1E-6, ALP=0.02):
        """
        plot a posterior pdf curve
        """
        fNLlist = np.arange(fNLmin, fNLmax, 1)
        plt.plot(fNLlist, self.posterior(fNLlist, A), linestyle=ls, linewidth=lw, label=r"$A="+str(A)+"$", color=clr)

    def plot_posteriors(self, Alist=[0.0, 0.04, 0.055, 0.06],
                              llist=["dashed", "dashdot", "dotted", "solid"],
                              clist=["black", "red", "blue", "green"],
                              LW=2.5, fNLmax=2000, ymin=2E-6, ymax=5E-6, ALP=0.02):
        """
        plot multiple posterior curves
        """
        for i in range(0, len(Alist)):
            self.plot_posterior(Alist[i], ls=llist[i], lw=LW, clr=clist[i], fNLmax=fNLmax, ymin=ymin, ALP=ALP)

        if (ymin!=ymax):
            plt.ylim(ymin, ymax)

        plt.xlabel(r"$|f_{\rm NL}|$")
        plt.ylabel(r"$p(|f_{\rm NL}|)$")
        plt.legend(fontsize=lf)
        plt.yscale('log')

    def plot_combined_posteriors(self, Alist=[0.0, 0.04, 0.055, 0.06],
                              llist=["dashed", "dashdot", "dotted", "solid"],
                              clist=["black", "red", "blue", "green"],
                              LW=2.5, fNLmax=800, ALP=0.02, mean=-100, sigma=100):
        """
        plot multiple posterior curves
        """
        fNLlist=np.arange(0.0, fNLmax, 1)
        for i in range(0, len(Alist)):
            plt.plot(fNLlist, self.combined_posterior(fNLlist, Alist[i], mean=mean, sigma=sigma), linestyle=llist[i], linewidth=LW, color=clist[i], label=r"$A="+str(Alist[i])+"$")

        plt.xlabel(r"$|f_{\rm NL}|$")
        plt.ylabel(r"$p(|f_{\rm NL}|)$")
        plt.legend(loc=3, fontsize=lf)

    def plot_combined_posteriors_withA0(self, A=0.055, A0list=[0.0, 0.02, 0.04, 0.04],
                                        Nlist=[10, 50, 50, 100],
                              llist=["dashed", "dashdot", "dotted", "solid"],
                              clist=["black", "red", "blue", "green"],
                              LW=2.5, fNLmax=800, ALP=0.02, mean=-100, sigma=100):
        """
        plot multiple posterior curves
        """
        fNLlist=np.arange(0.0, 50, 0.1)
        fNLlist=np.append(fNLlist, np.arange(50, fNLmax, 1))
        for i in range(0, len(A0list)):
            if (A0list[i]!=0.0):
                labl=r"$A_0="+str(A0list[i])+r",\;N_{\rm extra}="+str(Nlist[i])+"$"
            else:
                labl=r"$A_0\,{\rm not\,used}$"
            self.efolds=Nlist[i]
            plt.plot(fNLlist, self.combined_posteriors_withA0(fNLlist, A, A0list[i], mean=mean, sigma=sigma), linestyle=llist[i], linewidth=LW, color=clist[i], label=labl)
            if (A0list[i]>0.0):
                sd=self.sddr_withA0(A, A0list[i], Nlist[i])
            else:
                sd=self.sddr(A)
            print (sd, np.log(sd))

        plt.xlabel(r"$|f_{\rm NL}|$")
        plt.ylabel(r"$p(|f_{\rm NL}|)$")
        plt.title(r"$A="+str(A)+"$")
        plt.legend(loc=1, fontsize=lf)

    def plot_pvalue(self, A, ls="solid", lw=2, clr="b", fNLmin=0, fNLmax=500):
        fNLlist=np.arange(fNLmin, fNLmax, 1)
        pvalues=np.array([self.pvalue(A, fNL) for fNL in fNLlist])
        plt.plot(fNLlist, pvalues, linestyle=ls, linewidth=lw, label=r"$A="+str(A)+"$", color=clr)

    def plot_pvalues(self, Alist=[0.0, 0.040, 0.055, 0.060],
                           llist=["dashed", "dashdot", "dotted", "solid"],
                           clist=["black", "red", "blue", "green"],
                           LW=2.5, fNLmax=500):
        for i in range(0, len(Alist)):
            self.plot_pvalue(Alist[i], ls=llist[i], lw=LW, clr=clist[i], fNLmax=fNLmax)

        plt.axhline(0.005, ls="dashed", lw=1.5, color="black", alpha=0.35)
        plt.text(220, 0.0055, r"$2.8\sigma$")
        plt.axhline(0.001, ls="dashed", lw=1.5, color="black", alpha=0.35)
        plt.text(220, 0.0011, r"$3.3\sigma$")
        plt.xlabel(r"$|f_{\rm NL}|$")
        plt.ylabel(r"$p{\rm-value}$")
        plt.legend(loc=0, fontsize=lf)
        plt.ylim(1E-4, 1.)
        plt.yscale('log')
