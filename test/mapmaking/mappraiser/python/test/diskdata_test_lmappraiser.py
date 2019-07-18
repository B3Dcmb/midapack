# This is a test code for the MAPPRAISER python wrapper using data from TOAST
# stored on the disk. This is a beta version of the test: dependent on the data format,
# a specific data distribution configuration (however general enough for our purposes),
# and requires manual post operations to verify the accuracy of the map reconstruction.
# More systematic unit test cases will replace this code in the future.

#@Author: Hamza El Bouhargani <elbouha@apc.in2p3.fr>
#@Date: May 2019

# Import basic modules
from __future__ import division, absolute_import, print_function

from mpi4py import MPI
import ctypes as ct
import numpy as np

import mappraiser_wrapper as mappraiser

# Define MPI parameters
rank = MPI.COMM_WORLD.rank
print("My rank is {}".format(rank))
size = MPI.COMM_WORLD.size
print("Size of MPI processes is {}".format(size))

# Define paths to main data files
In_dir = "/global/cscratch1/sd/elbouha/data_TOAST/test4_clean/" # Path to the data (user dependent)

pix_name = "pixels_" # Reference to the pixel indices array (user dependent)
pixweights_name = "weights_" # Reference to the pixel weights array (user dependent)
signal_name = "pure_signal_" # Reference to the signal array (user dependent)
invtt_name = "inv_tt_x3" # Reference to the first row of the inverse noise time correlation matrix (user dependent)

ndet = size # number of detectors (each MPI process reads one detector full timestream)
Lambda = 2**8 # half-bandwidth of inverse noise time correlation matrix
Nnz = 3 # Solve for I, Q, and U parameters

# Read disk data
with open(In_dir+pix_name+str(rank)+".dat","rb") as file:
    pixels = np.fromfile(file, dtype=mappraiser.PIXEL_TYPE)
with open(In_dir+pixweights_name+str(rank)+".dat","rb") as file:
    pixweights = np.fromfile(file, dtype=mappraiser.WEIGHT_TYPE)
with open(In_dir+signal_name+str(rank)+".dat","rb") as file:
    signal = np.fromfile(file, dtype=mappraiser.SIGNAL_TYPE)
with open(In_dir+invtt_name+".bin","rb") as file:
    invtt = np.fromfile(file, dtype=mappraiser.INVTT_TYPE)
    invtt = invtt[:Lambda]

#Simple binning case
if Lambda == 1:
    invtt[0] = 1

# Run the map-making
mappraiser.MLmap(MPI.COMM_WORLD, ndet, len(signal), Nnz, pixels, pixweights, signal, Lambda, invtt)
