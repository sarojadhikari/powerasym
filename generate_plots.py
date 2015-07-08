"""
script to generate useful plots
"""
import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'figure.autolayout': True})
matplotlib.rcParams.update({'ytick.major.pad': 9})
matplotlib.rcParams.update({'xtick.major.pad': 7})

from padists import PowerAsymmetryDistribution

MAPS=10000

# generate fNLlocal Ai and A plots #
#----------------------------------#
pad = PowerAsymmetryDistribution(NMAPS=MAPS)
pad.set_fgnl([250, 500]); pad.read_data()

pad.plot_A(histtype='step')
plt.xlim(0.0, 0.12)
plt.savefig("plots/Adist.pdf")
plt.close()

pad.plot_Ai(histtype='step')
plt.xlim(-0.08, 0.08)
plt.savefig("plots/Aidist.pdf")
plt.close()

pad.set_fgnl([50, 100]); pad.read_data()
plt.figure(num=None, figsize=(8, 5))
pad.plot_A0()
plt.xlim(-0.4, 0.4); plt.ylim(0.1, 50);
plt.title(r"$N="+str(pad.efolds)+"$")
plt.savefig("plots/A0fNL.pdf")
plt.close()

# generate fNLlocal nodipole plots #
#----------------------------------#
#pad = PowerAsymmetryDistribution(NMAPS=MAPS)
#pad.set_fgnl([500]); pad.read_data()
#pad.clrs=['w', 'g', 'r', 'black']
#pad.plot_A()
#AG=np.load("data50/AdistG.npy")
#plt.hist(AG, bins=50, histtype='step', color='b', normed=True)

#plt.savefig("plots/Anodipole.pdf")
#plt.close()

#pad.plot_Ai()
#plt.savefig("plots/Ainodipole.pdf")
#plt.close()

# generate gNLlocal plots #
#-------------------------#
pad = PowerAsymmetryDistribution(typ='gNL', theory=True, NMAPS=MAPS)
pad.set_fgnl([10000, 50000]); pad.read_data()
plt.figure(num=None, figsize=(8,5))
pad.plot_A0()
plt.xlim(-0.01, 0.25)
plt.ylim(0.1, 100)
plt.title(r"$N="+str(pad.efolds)+"$")
plt.savefig("plots/A0gNL.pdf")
plt.close()

# generate C_1-A correlation #
#----------------------------#
pad = PowerAsymmetryDistribution(typ='fNL', theory=False)
pad.set_fgnl([500]); pad.read_data()

NPTS=2000

fNLC1s=np.array([pad.fgNLCls[0][i][1] for i in range(NPTS)])
plt.plot(pad.gA[:NPTS], fNLC1s, 'bo', alpha=0.2, label=r"$f_{\rm NL}=0$")
plt.plot(pad.fgNLA[0][:NPTS], fNLC1s, 'g*', alpha=0.2, label=r"$f_{\rm NL}=500$")
plt.yscale('log')
plt.xlabel(r'$A$')
plt.ylabel(r'$C_1$')
plt.ylim(1E-13, 1E-8)
plt.legend(loc=0)

plt.savefig("plots/AC1corr.pdf")
plt.close()

# generate correlation plot between G and fNLlocal case
#======================================================
AiG=np.load("data50/AidistG.npy")[0:NPTS]
AinG=np.load("data50/AidistfNL500.npy")[0:NPTS]
plt.plot(AiG, AinG, 'go', alpha=0.2)
plt.xlabel(r"$A_i(f_{\rm NL}=0)$")
plt.ylabel(r"$A_i(f_{\rm NL}=500)$")
plt.savefig("plots/AiGnG.pdf")
plt.close()

AG=np.load("data50/AdistG.npy")[0:NPTS]
AnG=np.load("data50/AdistfNL500.npy")[0:NPTS]
plt.plot(AG, AnG, 'bo', alpha=0.1)
plt.xlabel(r"$A (f_{\rm NL}=0)$")
plt.ylabel(r"$A (f_{\rm NL}=500)$")
plt.savefig("plots/AGnG.pdf")
plt.close()


# generate posterior plots
# ========================

from posterior import PosteriorfNL

pf = PosteriorfNL()
plt.figure(num=None, figsize=(8,5))
pf.plot_posteriors(ymax=0.01)
plt.savefig("plots/posteriorA1.pdf")

