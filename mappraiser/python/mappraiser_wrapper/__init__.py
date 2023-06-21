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


############################################################
# MLmap routine
############################################################

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
    ct.c_int,  # gap_stgy
    ct.c_uint64,  # realization
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # data_size_proc
    ct.c_int,  # nb_blocks_loc
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # local_blocks_sizes
    ct.c_double,  # sample_rate
    npc.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS"),  # detindxs
    npc.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS"),  # obsindxs
    npc.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS"),  # telescopes
    ct.c_int,  # Nnz
    npc.ndpointer(dtype=PIXEL_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # pixels
    npc.ndpointer(dtype=WEIGHT_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # pixweights
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # signal
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # noise
    ct.c_int,  # lambda
    npc.ndpointer(dtype=INVTT_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # inv_tt
    npc.ndpointer(dtype=INVTT_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # tt
]


def MLmap(
        comm,
        params,
        data_size_proc,
        nb_blocks_loc,
        local_blocks_sizes,
        detindxs,
        obsindxs,
        telescopes,
        nnz,
        pixels,
        pixweights,
        signal,
        noise,
        inv_tt,
        tt,
):
    """
    Compute the MLMV solution of the GLS estimator, assuming uniform detector weighting and a single PSD
    for all stationary intervals.
    (These assumptions will be removed in future updates)

    Parameters
    ----------
    * `comm`: Communicator over which data is distributed
    * `params`: Parameter dictionary
    * `data_size_proc`: Data sizes in full communicator
    * `nb_blocks_loc`: Number of local observations
    * `local_blocks_sizes`: Local data sizes
    * `nnz`: Number of non-zero elements per row
    * `pixels`: Pixel indices of non-zero values
    * `pixweights`: Corresponding matrix values
    * `signal`: Signal buffer
    * `noise`: Noise buffer
    * `inv_tt`: Inverse noise correlation
    * `tt`: Noise autocorrelation

    """
    if not available:
        raise RuntimeError("No libmappraiser available, cannot reconstruct the map")

    outpath = params["path_output"].encode("ascii")
    ref = params["ref"].encode("ascii")

    comm.Barrier()

    _mappraiser.MLmap(
        encode_comm(comm),
        outpath,
        ref,
        params["solver"],
        params["precond"],
        params["Z_2lvl"],
        params["ptcomm_flag"],
        params["tol"],
        params["maxiter"],
        params["enlFac"],
        params["ortho_alg"],
        params["bs_red"],
        params["nside"],
        params["gap_stgy"],
        params["realization"],
        data_size_proc,
        nb_blocks_loc,
        local_blocks_sizes,
        params["fsample"],
        detindxs,
        obsindxs,
        telescopes,
        nnz,
        pixels,
        pixweights,
        signal,
        noise,
        params["Lambda"],
        inv_tt,
        tt,
    )

    return


############################################################
# gap filling routine
############################################################

_mappraiser.gap_filling.restype = None
_mappraiser.gap_filling.argtypes = [
    MPI_Comm,  # comm
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # data_size_proc
    ct.c_int,  # nb_blocks_loc
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),  # local_blocks_sizes
    ct.c_int,  # Nnz
    npc.ndpointer(dtype=INVTT_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # tt
    npc.ndpointer(dtype=INVTT_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # inv_tt
    ct.c_int,  # lambda
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # noise
    npc.ndpointer(dtype=PIXEL_TYPE, ndim=1, flags="C_CONTIGUOUS"),  # indices
    ct.c_uint64,  # realization
    npc.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS"),  # detindxs
    npc.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS"),  # obsindxs
    npc.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS"),  # telescopes
    ct.c_double,  # sample_rate
]


def gap_filling(
        comm,
        data_size_proc,
        nb_blocks_loc,
        local_blocks_sizes,
        params,
        detindxs,
        obsindxs,
        telescopes,
        nnz,
        pixels,
        signal,
        inv_tt,
        tt,
):
    if not available:
        raise RuntimeError("No libmappraiser available, cannot perform gap-filling")

    comm.Barrier()

    _mappraiser.gap_filling(
        encode_comm(comm),
        data_size_proc,
        nb_blocks_loc,
        local_blocks_sizes,
        nnz,
        tt,
        inv_tt,
        params["Lambda"],
        signal,
        pixels,
        params["realization"],
        detindxs,
        obsindxs,
        telescopes,
        params["fsample"],
    )

    return


############################################################
# PSD routine
############################################################

_mappraiser.psd_from_tt.restype = None
_mappraiser.psd_from_tt.argtypes = [
    ct.c_int,  # fftlen
    ct.c_int,  # lambda
    ct.c_int,  # psdlen
    npc.ndpointer(dtype=np.double, ndim=1, flags="C_CONTIGUOUS"),  # tt
    npc.ndpointer(dtype=np.double, ndim=1, flags="C_CONTIGUOUS"),  # psd
]


def psd_from_tt(
        fftlen,
        Lambda,
        psdlen,
        tt,
        psd,
):
    # if not available:
    #     raise RuntimeError("No libmappraiser available, cannot call psd_from_tt")

    _mappraiser.psd_from_tt(
        fftlen,
        Lambda,
        psdlen,
        tt,
        psd,
    )

    return


############################################################
# Timestream generation routine
############################################################

_mappraiser.sim_noise_tod.restype = None
_mappraiser.sim_noise_tod.argtypes = [
    ct.c_int,  # samples
    ct.c_int,  # lambda
    npc.ndpointer(dtype=np.double, ndim=1, flags="C_CONTIGUOUS"),  # tt
    npc.ndpointer(dtype=np.double, ndim=1, flags="C_CONTIGUOUS"),  # buf
    ct.c_uint64,  # realization
    ct.c_uint64,  # detindx
    ct.c_uint64,  # obsindx
    ct.c_uint64,  # telescope
    ct.c_double,  # sample_rate
]


def sim_noise_tod(
        samples,
        Lambda,
        tt,
        buf,
        realization,
        detindx,
        obsindx,
        telescope,
        sample_rate,
):
    # if not available:
    #     raise RuntimeError("No libmappraiser available, cannot call psd_from_tt")

    _mappraiser.sim_noise_tod(
        samples,
        Lambda,
        tt,
        buf,
        realization,
        detindx,
        obsindx,
        telescope,
        sample_rate,
    )

    return
