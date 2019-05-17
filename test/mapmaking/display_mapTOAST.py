#!/usr/bin/python

# Import basic modules
from __future__ import division, absolute_import, print_function
import numpy as np
import healpy as hp
import pylab as pl

# Import array map from binary file
# with open("/global/cscratch1/sd/elbouha/data_clean_nside128/point_data_0_66_scan_0.dat",'rb') as file:
#     point = np.fromfile(file, dtype = 'int32')
# with open("/global/cscratch1/sd/elbouha/data_TOAST/test0/pixels_0.dat",'rb') as file:
#     point = np.fromfile(file, dtype = 'int32')

# correct_hits =hp.fitsfunc.read_map("/global/homes/e/elbouha/toast/examples/out_tiny_ground_simple_cori-haswell/out/madam_hmap_003.fits", nest=True)

with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_hits0.dat",'rb') as file:
    map_mask = np.fromfile(file, dtype = 'int32')
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_I.dat",'rb') as file:
    mapI = 1e6*np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_Q.dat",'rb') as file:
    mapQ = 1e6*np.fromfile(file)
with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_U.dat",'rb') as file:
    mapU = 1e6*np.fromfile(file)
# with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_hits.dat",'rb') as file:
#     hits = np.fromfile(file, dtype = 'int32')
# with open("/global/homes/e/elbouha/midapack/test/mapmaking/mapoutAll_cond.dat",'rb') as file:
#     cond = np.fromfile(file)

# print(mapI)
# mask = map_mask == 0
# maskI = mask
maskI = mapI == 0
# maskP = maskI
# maskI = hits<300
# maskI = cond < 1e-3
mapI[maskI] = np.nan
# maskP = hits<300
maskP = maskI
# mask = map == 0
mapQ[maskP] = np.nan
mapU[maskP] = np.nan
# cond[mask] =np.nan
# map = np.zeros_like(hits)
# for i in point:
#     map[int(i/3)] = int(i/3)
# index = np.nonzero(correct_hits)
# for i in index:
#     map[i] = i
# print(index)
# print(np.sum(hits))



# M = np.sum(hits)
# print("total time samples: {}".format(M))

# Load input map
from s4cmb.input_sky import HealpixFitsMap
#Sky in
nside_out = 512
path_to_in = "/global/homes/e/elbouha/toast/examples/data/ffp10_lensed_scl_100_nside0512.fits"
sky_in = 1e6*hp.fitsfunc.read_map(path_to_in, field = None, nest=True)



# Define Planck color map
from matplotlib.colors import ListedColormap
planck_cmap = ListedColormap(np.loadtxt("Planck_Parchment_RGB.txt")/255.)
planck_cmap.set_bad("gray")
planck_cmap.set_under("white")
cmap = planck_cmap

# Plot the maps
# pl.figure(2,figsize=(8, 8))
xsize = 700
# hp.gnomview(np.log10(hits+1), rot=[0, -30], xsize=xsize, reso=6.9,
#             title='nhits', notext=True, cmap=cmap)
# hp.gnomview(1e-1*hits, rot=[0, -30], xsize=xsize, reso = 1,
#             title='hits', notext=True, cmap=cmap, nest=True)
# pl.show()
# hp.gnomview(cond, rot=[0, -30], xsize=xsize, reso = 1.5,
#             title='rcond', notext=True, cmap=cmap, nest=True)
# pl.show()

pl.figure(1,figsize=(8, 8))
sky_in[0][maskI] = np.nan
hp.gnomview(sky_in[0], rot=[0, -30], xsize=xsize, sub=331, reso = 1.3,
            title='Input Temperature', notext=True, cmap=cmap, nest=True)
hp.gnomview(mapI, rot=[0, -30], xsize = xsize, sub = 332, reso = 1.3, title="Reconstructed Temperature map", notext=True, cmap=cmap, nest=True)
hp.gnomview(sky_in[0] - mapI, rot=[0, -30], xsize=xsize, sub=333, reso = 1.3,
            title='difference', notext=True, cmap = cmap, nest=True)
# sky_in.Q = np.zeros_like(sky_in.I)
sky_in[1][maskP] = np.nan
hp.gnomview(sky_in[1], rot=[0, -30], xsize=xsize, sub=334, reso = 1.3,
            title='Input Q', notext=True, cmap=cmap, nest=True)
hp.gnomview(mapQ, rot=[0, -30], xsize = xsize, sub = 335, reso = 1.3, title="Reconstructed Q map", notext=True, cmap=cmap, nest=True)
hp.gnomview(sky_in[1] - mapQ, rot=[0, -30], xsize=xsize, sub=336, reso = 1.3,
            title='difference', notext=True, cmap = cmap, nest=True)
# sky_in.U = np.zeros_like(sky_in.I)
sky_in[2][maskP] = np.nan
hp.gnomview(sky_in[2], rot=[0, -30], xsize=xsize,sub=337, reso = 1.3,
            title='Input U', notext=True, cmap=cmap, nest=True)
hp.gnomview(mapU, rot=[0, -30], xsize = xsize, sub = 338, reso = 1.3, title="Reconstructed U map", notext=True, cmap=cmap, nest=True)
hp.gnomview(sky_in[2] - mapU, rot=[0, -30], xsize=xsize, sub=339, reso = 1.3,
            title='difference', notext=True, cmap = cmap, nest=True)

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
