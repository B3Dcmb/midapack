from __future__ import division
from __future__ import print_function

import ctypes as ct
import ctypes.util as ctu
import os
import sys

import numpy as np
import numpy.ctypeslib as npc

from mpi4py import MPI

SIGNAL_TYPE = np.float64
PIXEL_TYPE = np.int32
WEIGHT_TYPE = np.float64
INVTT_TYPE = np.float64
TIMESTAMP_TYPE = np.float64
PSD_TYPE = np.float64

try:
    _mappraiser = ct.CDLL("libmappraiser.so")
except OSError:
    path = ctu.find_library("mappraiser")
    if path is not None:
        _mappraiser = ct.CDLL(path)

available = _mappraiser is not None

try:
    if MPI._sizeof(MPI.Comm) == ct.sizeof(ct.c_int):
        MPI_Comm = ct.c_int
    else:
        MPI_Comm = ct.c_void_p
except Exception as e:
    raise Exception(
        'Failed to set the portable MPI communicator datatype: "{}". '
        "MPI4py is probably too old. ".format(e)
    )


def encode_comm(comm):
    comm_ptr = MPI._addressof(comm)
    return MPI_Comm.from_address(comm_ptr)

_mappraiser.MLmap.restype = None
_mappraiser.MLmap.argtypes =[
    MPI_Comm, #comm
    ct.c_char_p, #ref
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), #data_size_proc
    ct.c_int, #nb_blocks_loc
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), #local_blocks_sizes
    ct.c_int, #Nnz
    npc.ndpointer(dtype=PIXEL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=WEIGHT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_int, #lambda
    npc.ndpointer(dtype=INVTT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
]

def MLmap(comm, ref, data_size_proc, nb_blocks_loc, local_blocks_sizes, Nnz, pixels, pixweights, signal, noise, Lambda, invtt):
    """
    Compute the MLMV solution of the GLS estimator, assuming uniform detector weighting and a single PSD
    For all stationary intervals. (These assumptions will be removed in future updates)
    """
    if not available:
        raise RuntimeError("No libmappraiser available, cannot reconstruct the map")
    ref = ref.encode('ascii')
    _mappraiser.MLmap(encode_comm(comm), ref, data_size_proc, nb_blocks_loc, local_blocks_sizes, Nnz, pixels, pixweights, signal, noise, Lambda, invtt)
    return
