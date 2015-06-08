"""
script to generate useful plots
"""
import numpy as np
import matplotlib.pyplot as plt
from padists import PowerAsymmetryDistribution

# generate fNLlocal Ai and A plots #
#----------------------------------#
pad = PowerAsymmetryDistribution()
pad.set_fgnl([250, 500]); pad.read_data()

pad.plot_A()
plt.savefig("plots/Adist.pdf")
plt.close()

pad.plot_Ai()
plt.savefig("plots/Aidist.pdf")
plt.close()

# generate fNLlocal nodipole plots #
#----------------------------------#
pad = PowerAsymmetryDistribution(datafolder="datanodipole")
pad.set_fgnl([500]); pad.read_data()

pad.plot_A()
plt.savefig("plots/Anodipole.pdf")
plt.close()

pad.plot_Ai()
plt.savefig("plots/Ainodipole.pdf")
plt.close()

# generate gNLlocal plots #
#-------------------------#
pad = PowerAsymmetryDistribution(typ='gNL', theory=False)
pad.set_fgnl([1000000, 5000000, 10000000]); pad.read_data()

pad.plot_A()
plt.savefig("plots/AgNLdist.pdf")
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