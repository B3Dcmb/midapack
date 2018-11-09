#!/usr/bin/python

# Import basic modules
from __future__ import division, absolute_import, print_function
import numpy as np
import healpy as hp
import pylab as pl

# Import array map from binary file
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_hits0.dat",'rb') as file:
    map_mask = np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_I.dat",'rb') as file:
    mapI = np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_Q.dat",'rb') as file:
    mapQ = np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_U.dat",'rb') as file:
    mapU = np.fromfile(file)
# with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_hits.dat",'rb') as file:
#     hits = np.fromfile(file, dtype = 'int32')
# with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_cond.dat",'rb') as file:
#     cond = np.fromfile(file)

# print(map)
# mask = map_mask == 0
maskI = mapI == 0
maskP = maskI
# maskI = hits<10
mapI[maskI] = np.nan
# maskP = hits<10
# mask = map == 0
mapQ[maskP] = np.nan
mapU[maskP] = np.nan
# cond[maskP] =np.nan

# M = np.sum(hits)
# print("total time samples: {}".format(M))

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
# hp.gnomview(hits, rot=[0, -57.5], xsize=xsize, reso=6.9,
#             title='nhits', notext=True, cmap=cmap)
# hp.gnomview(cond-hits, rot=[0, -57.5], xsize=xsize, reso=6.9,
#             title='nhits(Cos^2+sin^2-1) ', notext=True, cmap=cmap)
# pl.show()

pl.figure(1,figsize=(8, 8))
sky_in.I[maskI] = np.nan
hp.gnomview(sky_in.I, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=331,
            title='Input Temperature', notext=True, min=-250, max=250, cmap=cmap)
hp.gnomview(mapI, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 332, title="Reconstructed Temperature map", notext=True,max = 250, min = -250, cmap=cmap)
hp.gnomview(np.log10(np.abs(sky_in.I - mapI)), rot=[0, -57.5], xsize=xsize, reso=6.9, sub=333,
            title='Log difference', notext=True, min=-16 ,max=0, cmap = cmap)
# sky_in.Q = np.zeros_like(sky_in.I)
sky_in.Q[maskP] = np.nan
hp.gnomview(sky_in.Q, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=334,
            title='Input Q', notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(mapQ, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 335, title="Reconstructed Q map", notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(np.log10(np.abs(sky_in.Q - mapQ)), rot=[0, -57.5], xsize=xsize, reso=6.9, sub=336,
            title='Log difference', notext=True, min=-16 ,max=0, cmap = cmap)
# sky_in.U = np.zeros_like(sky_in.I)
sky_in.U[maskP] = np.nan
hp.gnomview(sky_in.U, rot=[0, -57.5], xsize=xsize, reso=6.9, sub=337,
            title='Input U', notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(mapU, rot=[0, -57.5], xsize = xsize,reso=6.9, sub = 338, title="Reconstructed U map", notext=True, min=-15, max=15, cmap=cmap)
hp.gnomview(np.log10(np.abs(sky_in.U - mapU)), rot=[0, -57.5], xsize=xsize, reso=6.9, sub=339,
            title='Log difference', notext=True, min=-16 ,max=0, cmap = cmap)

pl.show()

# aQ = sky_in.Q - mapQ
# aQ = aQ[~np.isnan(aQ)]
#
# aI = sky_in.I - mapI
# aI = aI[~np.isnan(aI)]
#
# aU = sky_in.U - mapU
# aU = aU[~np.isnan(aU)]
#
# pl.figure(2,figsize=(9, 9))
#
# pl.hist(aI, bins='auto' , label = "Temperature, sigma = {}".format(np.std(aI)), histtype = 'step')
# pl.hist(aQ, bins='auto', label = "Q, sigma = {}".format(np.std(aQ)), histtype = 'step')
# pl.hist(aU, bins='auto' , label = "U, sigma = {}".format(np.std(aU)), histtype = 'step')
#
# pl.title("Histogram of maps residuals")
# pl.grid(True)
# pl.legend()
# pl.show()
