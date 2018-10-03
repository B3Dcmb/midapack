#!/usr/bin/python

# Import basic modules
from __future__ import division, absolute_import, print_function
import numpy as np
import healpy as hp
import pylab as pl

# Import array map from binary file
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_hits0.dat",'rb') as file:
    map_mask = np.fromfile(file, dtype = 'int32')
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_U0.dat",'rb') as file:
    map0 = np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_U1.dat",'rb') as file:
    map1 = np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_U2.dat",'rb') as file:
    map2 = np.fromfile(file)



# print(map)
mask = map_mask < 5
maskI = mask
maskP = mask
# maskI = hits<10
map0[maskI] = np.nan
# maskP = hits<10
# mask = map == 0
map1[maskP] = np.nan
map2[maskP] = np.nan

# Load input map
from s4cmb.input_sky import HealpixFitsMap
#Sky in
nside_out = 512
path_to_cls = '/global/homes/e/elbouha/s4cmb/s4cmb/data/test_data_set_lensedCls.dat'
sky_in = HealpixFitsMap(path_to_cls, do_pol=True, fwhm_in=3.5,
                        nside_in=nside_out, map_seed=5843787,
                        verbose=False, no_ileak=False, no_quleak=False)


# Define Planck color map
from matplotlib.colors import ListedColormap
planck_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
planck_cmap.set_bad("gray")
planck_cmap.set_under("white")
cmap = planck_cmap

# Plot the maps
# pl.figure(1,figsize=(8, 8))
xsize = 700

pl.figure(1,figsize=(8, 8))
sky_in.U[maskI] = np.nan
hp.gnomview(sky_in.U, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=331,
            title='Input U map', notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(map0, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 332, title="Observed U map - run0", notext=True,max = 15, min = -15, cmap=cmap)
hp.gnomview(sky_in.U - map0, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=333,
            title='Noise map', notext=True, min=-25, max=25, cmap = cmap)
# sky_in.U = np.zeros_like(sky_in.U)
sky_in.U[maskP] = np.nan
hp.gnomview(sky_in.U, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=334,
            title='Input U map', notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(map1, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 335, title="Observed U map - run1", notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(sky_in.U - map1, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=336,
            title='Noise map', notext=True, min=-25, max=25, cmap = cmap)
# sky_in.U = np.zeros_like(sky_in.U)
sky_in.U[maskP] = np.nan
hp.gnomview(sky_in.U, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=337,
            title='Input U map', notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(map2, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 338, title="Observed U map - run2", notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(sky_in.U - map2, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=339,
            title='Noise map', notext=True, min=-25, max=25, cmap = cmap)

# pl.show()

a1 = sky_in.U - map1
a1 = a1[~np.isnan(a1)]

a0 = sky_in.U - map0
a0 = a0[~np.isnan(a0)]

a2 = sky_in.U - map2
a2 = a2[~np.isnan(a2)]

pl.figure(2,figsize=(9, 9))

pl.hist(a0, bins='auto' , label = "run0, sigma = {}".format(np.std(a0)), histtype = 'step')
pl.hist(a1, bins='auto', label = "run1, sigma = {}".format(np.std(a1)), histtype = 'step')
pl.hist(a2, bins='auto' , label = "run2, sigma = {}".format(np.std(a2)), histtype = 'step')

pl.title("Histogram of Q noise maps")
pl.grid(True)
pl.legend()
pl.show()
