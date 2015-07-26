"""
script to generate useful plots
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.rcParams.update({'ytick.major.pad': 9})
matplotlib.rcParams.update({'xtick.major.pad': 7})
#matplotlib.rc('tick_params', color="gray")
matplotlib.rc('axes', edgecolor="gray")

from padists import PowerAsymmetryDistribution

MAPS=10000

# generate fNLlocal Ai and A plots #
#----------------------------------#
xsize=7.0
ysize=5.6

def figurebox():
    plt.figure(num=None, figsize=(xsize, ysize))
    plt.tick_params(which='both', color="gray")

pad = PowerAsymmetryDistribution(NMAPS=MAPS)
pad.set_fgnl([250, 500]); pad.read_data()
figurebox()
pad.plot_A(histtype='step')
plt.xlim(0.0, 0.12)
plt.yticks([0, 10, 20, 30, 40])
plt.savefig("plots/Adist.pdf")
plt.close()

figurebox()
pad.plot_Ai(histtype='step')
plt.xlim(-0.1, 0.1)
plt.ylim(0.0, 32)
plt.yticks([0, 10, 20, 30])
plt.savefig("plots/Aidist.pdf")
plt.close()

pad.set_fgnl([50, 100]); pad.read_data()
figurebox()
pad.plot_A0()
plt.xlim(-0.4, 0.4); plt.ylim(0.3, 50);
plt.xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
plt.title(r"$N="+str(pad.efolds)+"$")
plt.savefig("plots/A0fNL.pdf")
plt.close()

# generate fNLlocal nodipole plots #
#----------------------------------#
pad = PowerAsymmetryDistribution(NMAPS=MAPS)
pad.set_fgnl([500]); pad.read_data()
figurebox()
pad.plot_A(histtype='step')
plt.xlim(0.0, 0.12)
plt.yticks([0,10, 20, 30, 40])
Anodipole=np.load("dataND50/AdistfNL500.npy")
plt.hist(Anodipole, histtype='step', bins=50, lw=2.0, color="r", normed=True, label=r"$f_{\rm NL}=500,\; C_1=0$")
plt.legend(fontsize=19)
plt.xlabel(r"$A$")
plt.ylabel(r"$p(A)$")
plt.savefig("plots/Anodipole.pdf")
plt.close()

# generate gNLlocal plots #
#-------------------------#
#pad = PowerAsymmetryDistribution(typ='gNL', theory=True, NMAPS=MAPS)
#pad.set_fgnl([10000, 50000]); pad.read_data()
#plt.figure(num=None, figsize=(8,5))
#pad.plot_A0()
#plt.xlim(-0.01, 0.25)
#plt.ylim(0.1, 100)
#plt.title(r"$N="+str(pad.efolds)+"$")
#plt.savefig("plots/A0gNL.pdf")
#plt.close()

# generate C_1-A correlation #
#----------------------------#
#pad = PowerAsymmetryDistribution(typ='fNL', theory=False)
#pad.set_fgnl([500]); pad.read_data()

NPTS=2000
pad.set_fgnl([500]); pad.read_data()
fNLC1s=np.array([pad.fgNLCls[0][i][1] for i in range(NPTS)])
figurebox()
plt.plot(pad.gA[:NPTS], fNLC1s, 'bo', alpha=0.1, label=r"$f_{\rm NL}=0$")
plt.plot(pad.fgNLA[-1][:NPTS], fNLC1s, 'g*', alpha=0.1, label=r"$f_{\rm NL}=500$")
plt.yscale('log')
plt.xlabel(r'$A$')
plt.ylabel(r'$C_1$')
plt.ylim(1E-13, 1E-8)
plt.xticks([0.0, 0.04, 0.08, 0.12])
plt.legend(loc=0)

plt.savefig("plots/AC1corr.pdf")
plt.close()

# generate correlation plot between G and fNLlocal case
#======================================================
AiG=np.load("data50/AidistG.npy")[0:NPTS]
AinG=np.load("data50/AidistfNL500.npy")[0:NPTS]
figurebox()
plt.plot(AiG, AinG, 'go', alpha=0.1)
plt.xlabel(r"$A_i(f_{\rm NL}=0)$")
plt.ylabel(r"$A_i(f_{\rm NL}=500)$")
plt.xticks([-0.06,-0.03, 0.0, 0.03, 0.06])
plt.yticks([-0.12, -0.06, 0.0, 0.06, 0.12])
plt.savefig("plots/AiGnG.pdf")
plt.close()

AG=np.load("data50/AdistG.npy")[0:NPTS]
AnG=np.load("data50/AdistfNL500.npy")[0:NPTS]
figurebox()
plt.plot(AG, AnG, 'bo', alpha=0.1)
plt.xlabel(r"$A (f_{\rm NL}=0)$")
plt.ylabel(r"$A (f_{\rm NL}=500)$")
plt.savefig("plots/AGnG.pdf")
#plt.close()


# generate posterior plots
# ========================
from posterior import PosteriorfNL

pf = PosteriorfNL()
figurebox()
pf.plot_posteriors(ymax=0.01)
plt.savefig("plots/posteriorA1.pdf")
plt.close

figurebox()
pf.plot_pvalues()
plt.xlim(0, 500)
plt.savefig("plots/pvalues.pdf")
plt.close()

figurebox()
pf.plot_combined_posteriors()
plt.xlim(0, 500)
plt.yscale('log')
plt.ylim(1E-5, 0.01)
plt.savefig("plots/comb_posteriors.pdf")
plt.close()

figurebox()
pf.plot_combined_posteriors_withA0(A0list=[0.0, 0.02, 0.04, 0.04], Nlist=[10, 50, 50, 100])
plt.legend(loc=0, fontsize=18)
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.00001, 0.01)
plt.savefig("plots/comb_A0_posteriors.pdf")