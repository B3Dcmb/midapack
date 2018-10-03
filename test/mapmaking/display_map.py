#!/usr/bin/python

# Import basic modules
from __future__ import division, absolute_import, print_function
import numpy as np
import healpy as hp
import pylab as pl

# Import array map from binary file
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll1.dat",'rb') as file:
    map_rand = np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll12.dat",'rb') as file:
    map = np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll.dat",'rb') as file:
    map24 = np.fromfile(file)

# print(map)
mask = map_rand == 0
map_rand[mask] = np.nan
# mask = map == 0
map[mask] = np.nan
map24[mask] = np.nan


# Load input map
from s4cmb.input_sky import HealpixFitsMap
#Sky in
nside_out = 512
path_to_cls = '/global/homes/e/elbouha/s4cmb/s4cmb/data/test_data_set_lensedCls.dat'
sky_in = HealpixFitsMap(path_to_cls, do_pol=False, fwhm_in=3.5,
                        nside_in=nside_out, map_seed=5843787,
                        verbose=False, no_ileak=False, no_quleak=False)


# Define Planck color map
from matplotlib.colors import ListedColormap
planck_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
planck_cmap.set_bad("gray")
planck_cmap.set_under("white")
cmap = planck_cmap

# Plot the maps
pl.figure(1,figsize=(8, 8))
xsize = 700
sky_in.I[mask] = np.nan
hp.gnomview(sky_in.I, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=331,
            title='Input', notext=True, min=-250, max=250, cmap=cmap)
hp.gnomview(map_rand, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 332, title="Observed map 1 nces with identity Nt", notext=True,max = 250, min = -250, cmap=cmap)
hp.gnomview(sky_in.I - map_rand, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=333,
            title='Difference', notext=True, max=25, min=-25, cmap = cmap)

hp.gnomview(sky_in.I, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=334,
            title='Input', notext=True, min=-250, max=250, cmap=cmap)
hp.gnomview(map, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 335, title="Observed map 12 nces with identity Nt", notext=True, min=-250, max=250, cmap=cmap)
hp.gnomview(sky_in.I - map, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=336,
            title='Difference', notext=True, max=25, min=-25, cmap = cmap)

hp.gnomview(sky_in.I, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=337,
            title='Input', notext=True, min=-250, max=250, cmap=cmap)
hp.gnomview(map24, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 338, title="Observed map 12 nces-fs = 100 Hz with identity Nt", notext=True, min=-250, max=250, cmap=cmap)
hp.gnomview(sky_in.I - map24, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=339,
            title='Difference', notext=True, cmap = cmap)

# pl.show()

a = sky_in.I - map
a = a[~np.isnan(a)]

a_rand = sky_in.I - map_rand
a_rand = a_rand[~np.isnan(a_rand)]

a24 = sky_in.I - map24
a24 = a24[~np.isnan(a24)]

pl.figure(2,figsize=(9, 9))

pl.hist(a_rand, bins='auto' , label = "Single CES, sigma = {}".format(np.std(a_rand)), histtype = 'step')
pl.hist(a, bins='auto', label = "12 CES, sigma = {}".format(np.std(a)), histtype = 'step')
pl.hist(a24, bins='auto' , label = "12 CES - fs = 100 Hz, sigma = {}".format(np.std(a24)), histtype = 'step')

pl.title("Histogram of temperature maps residuals")
pl.grid(True)
pl.legend()
pl.show()
