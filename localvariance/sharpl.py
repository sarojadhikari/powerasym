# remove low l features from a map, use of mask is permitted
import healpy as hp
import numpy as np
UNSEEN=-1.6375e30
LMAX=4000
lREM=0
def remove_low_l(map1, mask1):
    map1m=hp.ma(map1)
    map1m.mask=np.logical_not(mask1)
    map2=np.ma.masked_values(map1, UNSEEN)
    alm2=hp.map2alm(map2, lmax=LMAX)
    fac2=np.array([1.]*(LMAX+1))
    fac2[0:lREM]=0.
    alm2n=hp.sphtfunc.almxfl(alm2, fac2, inplace=False)
    map2new=hp.ma(hp.alm2map(alm2n, 2048))
    map2new.mask=np.logical_not(mask1)
    return map2new
