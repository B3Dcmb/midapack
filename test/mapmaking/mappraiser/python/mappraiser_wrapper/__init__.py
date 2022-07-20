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
_mappraiser.MLmap.argtypes = [
    MPI_Comm,  # comm
    ct.c_char_p,  # outpath
    ct.c_char_p,  # ref
    ct.c_int,  # solver
    ct.c_int,  # precond
    ct.c_int,  # Z_2lvl
    ct.c_int,  # pointing_commflag
    ct.c_double,  # tol
    ct.c_int,  # maxIter
    ct.c_int,  # enlFac
    ct.c_int,  # ortho_alg
    ct.c_int,  # bs_red
    ct.c_int,  # nside
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # data_size_proc
    ct.c_int,  # nb_blocks_loc
    # local_blocks_sizes
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_int,  # Nnz
    npc.ndpointer(dtype=PIXEL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=WEIGHT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_int,  # lambda
    npc.ndpointer(dtype=INVTT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
]


def MLmap(
    comm,
    params,
    data_size_proc,
    nb_blocks_loc,
    local_blocks_sizes,
    nnz,
    pixels,
    pixweights,
    signal,
    noise,
    Lambda,
    invtt,
):
    """
    Compute the MLMV solution of the GLS estimator, assuming uniform detector weighting and a single PSD for all stationary intervals.
    (These assumptions will be removed in future updates)

    Args:
        comm: communicator over which data is distributed
        params: parameter dictionary
        data_size_proc: data sizes in full communicator
        nb_blocks_loc: nbr of local observations
        local_blocks_sizes: local data sizes
        nnz: number of non-zero elements per row
        pixels: pixel indices of non-zero values
        pixweights: corresponding matrix values
        signal: signal buffer
        noise: noise buffer
        Lambda: toeplitz matrix half-bandwidth
        invtt: inverse noise weights

    """
    if not available:
        raise RuntimeError("No libmappraiser available, cannot reconstruct the map")

    outpath = params["path_output"].encode("ascii")
    ref = params["ref"].encode("ascii")

    _mappraiser.MLmap(
        encode_comm(comm),
        outpath,
        ref,
        params["solver"],
        params["precond"],
        params["Z_2lvl"],
        params["pointing_commflag"],
        params["tol"],
        params["maxiter"],
        params["enlFac"],
        params["ortho_alg"],
        params["bs_red"],
        params["nside"],
        data_size_proc,
        nb_blocks_loc,
        local_blocks_sizes,
        nnz,
        pixels,
        pixweights,
        signal,
        noise,
        Lambda,
        invtt,
    )

    return
