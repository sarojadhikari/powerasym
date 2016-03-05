"""
script to generate useful plots and save them
"""
import sys
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rcParams.update({'figure.autolayout': False})
matplotlib.rcParams.update({'ytick.major.pad': 9})
matplotlib.rcParams.update({'xtick.major.pad': 7})
matplotlib.rc('axes', edgecolor="gray")

from padists import PowerAsymmetryDistribution

MAPS=10000

# generate fNLlocal Ai and A plots #
#----------------------------------#
xsize=7.0
ysize=5.6

def figurebox(xs=xsize, ys=ysize):
    plt.figure(num=None, figsize=(xs, ys))
    plt.tick_params(which='both', color="gray")

pad = PowerAsymmetryDistribution(NMAPS=MAPS)
pad.set_fgnl([250, 500]); pad.read_data()
figurebox()
pad.plot_A(histtype='stepfilled')
plt.xlim(0.0, 0.12)
plt.yticks([0, 10, 20, 30, 40])
plt.savefig("plots/Adist.pdf")
plt.close()

figurebox()
pad.plot_Ai(histtype='stepfilled')
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
plt.title(r"$N_{\rm extra}="+str(pad.efolds)+"$", y=1.02)
plt.savefig("plots/A0fNL.pdf")
plt.close()

# generate posterior plots
# ========================
from posterior import PosteriorfNL

pf = PosteriorfNL()
figurebox()
pf.plot_posteriors(ymax=0.02)
plt.xticks([0, 400, 800, 1200, 1600, 2000])
plt.tight_layout()
plt.savefig("plots/posteriorA1.pdf")
plt.close()

figurebox()
pf.plot_pvalues()
plt.xlim(0, 500)
plt.ylim(2.E-4, 2.0)
plt.tight_layout()
plt.savefig("plots/pvalues.pdf")
plt.close()

figurebox()
pf.plot_combined_posteriors()
plt.xlim(0, 500)
plt.yscale('log')
plt.ylim(2E-5, 0.02)
plt.tight_layout()
plt.savefig("plots/comb_posteriors.pdf")
plt.close()

sys.exit(0)
figurebox()
pf.plot_combined_posteriors_withA0(A0list=[0.0, 0.02, 0.04, 0.04], Nlist=[10, 50, 50, 100])
plt.legend(loc=0, fontsize=18)
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.00001, 0.01)
plt.savefig("plots/comb_A0_posteriors.pdf")
