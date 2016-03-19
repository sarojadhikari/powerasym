""" LOCAL VARIANCE TOOLS for Power Asymmetry Analysis in CMB data
    Some of the code was used for the work arXiv:1408.5396
"""

import healpy as hp
import numpy as np

def nside(deg):
    """ return the NSIDE for placing the center of the disk of radius [deg] that covers the whole sky
    """
    if (deg>16.0): return 4
    elif (deg>=4.0): return 16
    elif (deg>=2.0): return 32
    elif (deg>=1.0): return 64
    elif (deg>=0.4): return 128
    elif (deg>=0.17): return 256
    else: return 512

def deg2rad(deg):
    """ convert degree to radians
    """
    return 0.0174532925*deg

def local_power_spectrum(map1, mask1, deg, LMAX=256):
    """ return the local power spectrum for each disk given by the radius [deg]
    """
    mp1=hp.ma(map1)
    mp1.mask=np.logical_not(mask1)
    NSIDE1=hp.npix2nside(len(map1))
    NSIDE2=nside(deg)
    NPIX2=hp.nside2npix(NSIDE2)
    NPIX1=hp.nside2npix(NSIDE1)
    mp2=[0]*NPIX2
    for pix2 in range(0, NPIX2):
        newmask=np.array([0.]*NPIX1)
        disc1=hp.query_disc(nside=NSIDE1, vec=hp.pix2vec(NSIDE2, pix2), radius=deg2rad(deg))
        newmask[disc1]=1.
        newmask=newmask*mask1
        mp1.mask=np.logical_not(mask1)
        mp2[pix2]=hp.anafast(mp1, lmax=LMAX)
    return mp2


def local_mean_map(map1, mask1, deg):
    """ return the local mean map for map1 with mask mask1
    """
    mp1=hp.ma(map1)
    mp1.mask=np.logical_not(mask1)
    NSIDE1=hp.npix2nside(len(map1))
    NSIDE2=nside(deg)
    NPIX2=hp.nside2npix(NSIDE2)
    mp2=np.array([0.]*NPIX2)
    for pix2 in range(0, NPIX2):
        disc1=hp.query_disc(nside=NSIDE1, vec=hp.pix2vec(NSIDE2, pix2), radius=deg2rad(deg))
        mp2[pix2]=np.mean(mp1[disc1])
    return mp2

def local_variance_map(map1, mask1, deg):
    """ return the local variance map for map1 with mask mask1, given the error tolerance tol for mask2
    """
    mp1=hp.ma(map1)
    mp1.mask=np.logical_not(mask1)

    NSIDE1=hp.npix2nside(len(map1))
    NSIDE2=nside(deg)
    NPIX2=hp.nside2npix(NSIDE2)
    #mask2=np.round(hp.ud_grade(mask1, nside_out=NSIDE2)+tol)
    mp2=np.array([0.]*NPIX2)
    for pix2 in range(0, NPIX2):
        disc1=hp.query_disc(nside=NSIDE1, vec=hp.pix2vec(NSIDE2, pix2), radius=deg2rad(deg))
        mp2[pix2]=np.var(mp1[disc1])
    #varmp2.mask=np.logical_not(mask2)
    return mp2

def remove_low_l(map1, mask1, LMAX):
    """ return low_l removed map(masked)
    """
    mp1=hp.ma(map1)
    mp1.mask=np.logical_not(mask1)
    NSIDE1=hp.npix2nside(len(map1))
    alms=hp.map2alm(mp1, lmax=3072)
    alms[0:(LMAX+1)*(LMAX+2)/2]=0.
    newmap=hp.ma(hp.alm2map(alms, nside=NSIDE1))
    newmap.mask=np.logical_not(mask1)
    return newmap

def cos_modulation_map(nside, dvec):
    NSIDE=nside
    NPIX=hp.nside2npix(NSIDE)
    mmap=[np.dot(dvec, hp.pix2vec(NSIDE,i)) for i in range(NPIX)]
    return mmap

def add_doppler(map1, A, mmap):
    """ return a map that has the doppler dipole added to the temperature map given by the amplitude A and direction (angle); the amplitdue A should inlclude the frequency dependent boost factor b_nu.
    """
    NPIX=len(map1)
    NSIDE=hp.npix2nside(NPIX)
    newmap=np.array([0.]*NPIX)
    for i in range(NPIX):
        pixvec=hp.pix2vec(NSIDE, i)
        newmap[i]=map1[i]*(1+A*np.dot(pixvec, dvec))
    return newmap

def add_doppler_smoothed(map1, A, mmap, deg):
    """ return a map that has the doppler dipole added but only for the smoothed temperature maps with fwhm
    """
    map2=hp.smoothing(map1, fwhm=deg2rad(deg))
    newmap=map1-map2+add_doppler(map2, A, mmap)
    return newmap
