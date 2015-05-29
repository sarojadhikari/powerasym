"""
class to read the power asymmetry data, manipulate if necessary and 
make meaningful plots
"""

import sys
import numpy as np
from scipy.stats import chi, norm, chi2
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.rcParams.update({'ytick.major.pad': 9})
matplotlib.rcParams.update({'xtick.major.pad': 7})

ALPHA=0.12
LW=2.0

def plot_hist(pls, data, BINS=40, ht='stepfilled', clr='b', labl="", alp=0.1):
    pls.hist(data, histtype=ht, bins=BINS, normed=True, color=clr, label=labl, linewidth=LW, alpha=alp)

def NtoSTR(num):
    if (num<100000):
        return "$"+str(num)+"$"
    else:
        index=int(np.log10(num))
        base=num/np.power(10, index)
        if base>1:
            return r"${0}\times 10^{1}$".format(base, index)
        else:
            return r"$10^{0}$".format(index)
    

class PowerAsymmetryDistribution(object):
    def __init__(self, typ='fNL', datafolder='data', NMAPS=10000, efolds=50, theory=True):
        """
        typ specifies where it is fNL or gNL; the gaussian distribution is read in both cases
        """
        self.TYPE=typ
        self.efolds=efolds
        self.basedir=datafolder+str(efolds)+"/"
        self.fgnls=None
        self.nmaps=NMAPS
        self.read_data()
        self.clrs=['b', 'g', 'r', 'black']
        self.ls=[':', '--', '-', 'o']
        self.theoryplot=theory
    
    def set_fgnl(self, fgnllist):
        self.fgnls=fgnllist
        
    def read_data(self):
        try:
            self.gA0=np.load(self.basedir+"A0distG.npy")
            self.gAi=np.load(self.basedir+"AidistG.npy")
            self.NSIMS=len(self.gAi)/3
            self.gA=np.load(self.basedir+"AdistG.npy")
            self.gCls=np.load(self.basedir+"ClsG.npy")[0:self.nmaps*101].reshape(self.nmaps, 101)
        except:
            print "cannot read the asymmetry distribution for the Gaussian case!"
            sys.exit(1)
        
        if self.fgnls != None:
            self.fgNLA0=[]
            self.fgNLAi=[]
            self.fgNLA=[]
            self.fgNLCls=[]
            
            for fgnl in self.fgnls:
                if self.TYPE=='fNL':
                    self.TYPELABEL=r'$f_{\rm NL}$'
                    self.A0const=0.00011257 # to get sigma(A0), multiply by fNL and sqrt(N)
                    self.A1const=0.04621/500.
                    self.fgNLA0.append(np.load(self.basedir+"A0distfNL"+str(fgnl)+".npy")[0:self.nmaps])
                    self.fgNLAi.append(np.load(self.basedir+"AidistfNL"+str(fgnl)+".npy")[0:3*self.nmaps])
                    self.fgNLA.append(np.load(self.basedir+"AdistfNL"+str(fgnl)+".npy")[0:self.nmaps])
                    self.fgNLCls.append(np.load(self.basedir+"ClsfNL"+str(fgnl)+".npy")[0:self.nmaps*101].reshape(self.nmaps, 101)) # 101 because the Cls are saved upto LMAX=100
                elif self.TYPE=='gNL':
                    #self.A0const=3.95E-08
                    self.A0const=7.47562E-9 # to get scale multiply by gNL^0.5 and N
                    self.A1const=0.023/5000000.
                    self.TYPELABEL=r'$g_{\rm NL}$'
                    self.fgNLA0.append(np.load(self.basedir+"A0distgNL"+str(fgnl)+".npy")[0:self.nmaps])
                    self.fgNLAi.append(np.load(self.basedir+"AidistgNL"+str(fgnl)+".npy")[0:3*self.nmaps])
                    self.fgNLA.append(np.load(self.basedir+"AdistgNL"+str(fgnl)+".npy")[0:self.nmaps])
                else:
                    print "type must be fNL or gNL"
                    sys.exit(1)

    def plot_A0(self):
        plot_hist(plt, self.gA0, clr=self.clrs[0], alp=1.0, ht='step')
        
        sG=np.sqrt(np.var(self.gA0))
        mG=np.mean(self.gA0)
        amin, amax=plt.xlim()
        alist=np.arange(3*amin, amax*3, amax/250.)
        theoryGdist=norm.pdf(alist, loc=mG, scale=sG)
        plt.plot(alist, theoryGdist, self.clrs[0], linestyle=self.ls[0], linewidth=LW, label=self.TYPELABEL+"$=0$")
        
        for i in range(len(self.fgnls)):
            if (self.theoryplot==False):
                lbl=self.TYPELABEL+"$=$"+NtoSTR(self.fgnls[i])
            else:
                lbl=None
            
            plot_hist(plt, self.fgNLA0[i], clr=self.clrs[i+1], alp=1.0, labl=lbl, ht='step')
            amin, amax=plt.xlim()
            alist=np.arange(3*amin, amax*3, amax/250.)
            if (self.theoryplot):
                if (self.TYPE=='fNL'):
                    theorynGdist=norm.pdf(alist, loc=mG, scale=np.sqrt((self.A0const*self.fgnls[i])**2.0*self.efolds+sG**2.0))
                if (self.TYPE=='gNL'):
                    scl=self.A0const*self.efolds*self.fgnls[i]
                    print scl
                    ngpdf=chi2.pdf(alist, 1, loc=mG, scale=scl)
                    theorynGdist=ngpdf
                    
                plt.plot(alist, theorynGdist, self.clrs[i+1], linestyle=self.ls[i+1], linewidth=LW, label=self.TYPELABEL+"="+NtoSTR(self.fgnls[i]))
            
        plt.xlabel(r'$A_0$')
        plt.ylabel(r'$f_{\rm norm}$')
        plt.yscale('log')
        plt.ylim(0.01, 100)
        plt.legend()
        
        
    def plot_Ai(self):
        plot_hist(plt, self.gAi, clr=self.clrs[0], alp=ALPHA)

        sG=np.sqrt(np.var(self.gAi))
        mG=np.mean(self.gAi)
        amin, amax = plt.xlim()
        alist=np.arange(2*amin, amax*2, amax/250.)
        theoryGdist=norm.pdf(alist, loc=mG, scale=sG)
        plt.plot(alist, theoryGdist,self.clrs[0], linestyle=self.ls[0], linewidth=LW, label=self.TYPELABEL+"=0")

        for i in range(len(self.fgnls)):
            if (self.theoryplot==False):
                lbl=self.TYPELABEL+"="+NtoSTR(self.fgnls[i])
            else:
                lbl=None
                
            plot_hist(plt, self.fgNLAi[i], clr=self.clrs[i+1], alp=ALPHA, labl=lbl)
            
            if (self.theoryplot):
                theorynGdist=norm.pdf(alist, loc=mG, scale=np.sqrt((self.A1const*self.fgnls[i])**2.0+sG**2.0))
                plt.plot(alist, theorynGdist, self.clrs[i+1], linestyle=self.ls[i+1], linewidth=LW, label=self.TYPELABEL+"="+NtoSTR(self.fgnls[i]))
        
        plt.xlim(-0.15, 0.15)
        plt.xlabel(r'$A_i$')
        plt.ylabel(r'$f_{\rm norm}$')
        plt.legend()
        
    def plot_A(self):
        plot_hist(plt, self.gA, clr=self.clrs[0], alp=ALPHA)

        sG=np.sqrt(np.var(self.gAi))
        mG=np.mean(self.gAi)
        amin, amax = plt.xlim()
        alist=np.arange(amin, amax*3, amax/250.)
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
        
        plt.xlim(0.0, 0.2)
        plt.xticks([0.05, 0.1, 0.15, 0.2])
        plt.xlabel(r'$A$')
        plt.ylabel(r'$f_{\rm norm}$')
        plt.legend()
    
        
        
    def plot_hist_scatter(self, dist1, dist2, BINS=40, label1=None, label2=None):
        """
        given two list of numbers dist1 and dist2 generate a scatter plot
        and also draw histograms on the two axes
        =======================================
        INCOMPLETE: DOES NOT WORK AT THE MOMENT
        =======================================
        """
        from matplotlib.ticker import NullFormatter
        
        nullfmt=NullFormatter()
        
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left+width+0.02
        
        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]
        
        plt.figure(1)
        
        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        
        axScatter.scatter(dist1, dist2, 'bo', alpha=0.12)
        
        axHistx.hist(dist1, bins=BINS, histtype='stepfilled', alpha=0.15)
        axHisty.hist(dist2, bins=BINS, orientation='horizontal', histtype='stepfilled', alpha=0.15)
        
        axHistx.set_xlim(axScatter.get_xlim())
        axHisty.set_ylim(axScatter.get_ylim())
        
        plt.show()
        
        
