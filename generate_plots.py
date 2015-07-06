"""
script to generate useful plots
"""
import numpy as np
import matplotlib.pyplot as plt
from padists import PowerAsymmetryDistribution

MAPS=10000

# generate fNLlocal Ai and A plots #
#----------------------------------#
pad = PowerAsymmetryDistribution(NMAPS=MAPS)
pad.set_fgnl([250, 500]); pad.read_data()

pad.plot_A()
plt.xlim(0.0, 0.12)
plt.savefig("plots/Adist.pdf")
plt.close()

pad.plot_Ai()
plt.xlim(-0.08, 0.08)
plt.savefig("plots/Aidist.pdf")
plt.close()

pad.set_fgnl([100, 250]); pad.read_data()
plt.figure(num=None, figsize=(8, 5))
pad.plot_A0()
plt.xlim(-0.75, 0.75); plt.ylim(0.1, 50);
plt.title(r"$N="+str(pad.efolds)+"$")
plt.savefig("plots/A0fNL.pdf")
plt.close()

# generate fNLlocal nodipole plots #
#----------------------------------#
#pad = PowerAsymmetryDistribution(datafolder="datanomonopoles/datanodipole", NMAPS=MAPS)
#pad.set_fgnl([500]); pad.read_data()

#pad.plot_A()
#plt.savefig("plots/Anodipole.pdf")
#plt.close()

#pad.plot_Ai()
#plt.savefig("plots/Ainodipole.pdf")
#plt.close()

# generate gNLlocal plots #
#-------------------------#
pad = PowerAsymmetryDistribution(typ='gNL', theory=True, NMAPS=MAPS)
pad.set_fgnl([10000, 100000]); pad.read_data()

pad.plot_A0()
plt.xlim(-0.01, 0.75)
plt.savefig("plots/A0gNLdist.pdf")
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
plt.ylim(2E-12, 3E-9)
plt.legend(loc=0)

plt.savefig("plots/AC1corr.pdf")
plt.close()