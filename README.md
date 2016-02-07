# Power Asymmetry from non-Gaussian Modulation #
![map.jpg](https://github.com/sarojadhikari/powerasym/blob/master/figures/map.jpg)

This project contains the code necessary to produce some of the numerical results in the paper [arXiv:1508.06489](http://arxiv.org/abs/1508.06489), in which we extrapolate the power spectrum to superhorizon scales and show that the presence of the background density dipole, in the presence of large enough local type non-Gaussianity, generates a power asymmetry in the CMB power spectrum.

The basic idea is very simple: generate a large number of CMB maps with the Sachs-Wolfe spectrum (including the dipole term C_1, which is usually dropped, as the temperature dipole measured is dominated by the Doppler dipole). Then, assume local type non-Gaussianity and generate non-Gaussian temperature fluctuation maps. This is simple when only working in a Sachs-Wolfe universe. This couples the small scale power spectrum with the background dipole and one therefore gets a power asymmetry whose strength depends on (i) the strength of the non-linearity (fNL), and (ii) the magnitude of the background dipole (which depends on the particular realization of the modes for a CMB sky)

The following are some useful files:

* **swmaps.py** generates the necessary Gaussian maps and calculates the relevant data: asymmetry amplitudes, background dipole in each map etc

* **ngsachswolfe.py** contains a useful class SachsWolfeMap that is used for the calculations done in swmaps

* **cmbutils.py** some utilities for use in the above two files but that may be useful at other places later.

* **padists.py** contains a class called PowerAsymmetryDistribution that can read the data generated by swmaps.py to generate useful plots

* **generate_plots.py** script to generate some useful plots
