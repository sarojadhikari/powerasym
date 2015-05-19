# Power Asymmetry from non-Gaussian Modulation #

This project contains the code necessary to produce some of the numerical results in the paper arXiv:150x.xxxxx

The basic idea is very simple: generate a large number of CMB maps with the Sachs Wolfe spectrum (including the dipole term C_1, which is usually as the temperature dipole measured is dominated by the Doppler dipole). Then, assume local type non-Gaussianity and generate non-Gaussian maps. This couples the small scale power spectrum with the background dipole and one therefore gets a power asymmetry whose strength depends on (i) the strength of the non-linearity (fNL, gNL), and (ii) the magnitude of the background dipole (which depends on the particular realization of C_1)

The following are some useful files:

* **swmaps.py** generates the necessary Gaussian maps and calculates the relevant data: asymmetry amplitudes, background dipole in each map etc

* **ngsachswolfe.py** contains a useful class SachsWolfeMap that is used for the calculations done in swmaps

* **cmbutils.py** some utilities for use in the above two files but that may be useful at other places later.