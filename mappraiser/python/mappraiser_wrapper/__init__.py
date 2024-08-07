from __future__ import division
from __future__ import print_function

import ctypes as ct
import ctypes.util as ctu
import os
import sys

import numpy as np
import numpy.ctypeslib as npc

from mpi4py import MPI

# A wrapper for ndpointer allowing to pass NULL pointers
def wrapped_ndptr(*args, **kwargs):
  base = npc.ndpointer(*args, **kwargs)
  def from_param(cls, obj):
    if obj is None:
      return obj
    return base.from_param(obj)
  return type(base.__name__, (base,), {'from_param': classmethod(from_param)})

SIGNAL_TYPE = np.float64
PIXEL_TYPE = np.int32
WEIGHT_TYPE = np.float64
INVTT_TYPE = np.float64
TIMESTAMP_TYPE = np.float64
PSD_TYPE = np.float64

#array_ptrs_type = npc.ndpointer(dtype=np.uintp, ndim=1, flags='C')
array_ptrs_type = wrapped_ndptr(dtype=np.uintp, ndim=1, flags='C')

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
    ct.c_char_p, #outpath
    ct.c_char_p, #ref
    ct.c_int, #solver
    ct.c_int,  # precond
    ct.c_int,  # Z_2lvl
    ct.c_int, #pointing_commflag
    ct.c_double, #tol
    ct.c_int, #maxIter
    ct.c_int, #enlFac
    ct.c_int, #ortho_alg
    ct.c_int, #bs_red
    ct.c_int, #nside
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), #data_size_proc
    ct.c_int, #nb_blocks_loc
    # local_blocks_sizes
    npc.ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_int, #Nnz
    npc.ndpointer(dtype=PIXEL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=WEIGHT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    npc.ndpointer(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_int, #lambda
    npc.ndpointer(dtype=INVTT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
]

_mappraiser.MTmap.restype = None
_mappraiser.MTmap.argtypes =[
    MPI_Comm, #comm
    ct.c_char_p, #outpath
    ct.c_char_p, #ref
    ct.c_int, #solver
    ct.c_int,  # precond
    ct.c_int,  # Z_2lvl
    ct.c_int, #pointing_commflag
    ct.c_int, #npoly
    ct.c_int, #nhwp
    ct.c_double, #hwpss-base
    ct.c_int, #sss
    ct.c_int, #sbins
    ct.c_double, #tol
    ct.c_int, #maxIter
    ct.c_int, #enlFac
    ct.c_int, #ortho_alg
    ct.c_int, #bs_red
    ct.c_int, #nside
    array_ptrs_type, #sweeptstamps
    wrapped_ndptr(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), #nsweeps
    array_ptrs_type, #az
    wrapped_ndptr(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"), #az_min
    wrapped_ndptr(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"), #az_max
    array_ptrs_type, # hwp_angle
    ct.c_int, #nces
    wrapped_ndptr(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), #data_size_proc
    ct.c_int, #nb_blocks_loc
    wrapped_ndptr(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"), #local_blocks_sizes
    ct.c_int, #Nnz
    wrapped_ndptr(dtype=PIXEL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    wrapped_ndptr(dtype=WEIGHT_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    wrapped_ndptr(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    wrapped_ndptr(dtype=SIGNAL_TYPE, ndim=1, flags="C_CONTIGUOUS"),
    ct.c_double, #samplerate
    npc.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
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
    invtt,
):
    """
    Compute the GLS estimator assuming a Toeplitz-structured noise covariance.

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
    * `invtt`: Inverse noise weights

    """
    if not available:
        raise RuntimeError("No libmappraiser available, cannot reconstruct the map")

    outpath = params["path_output"].encode("ascii")
    ref = params["ref"].encode("ascii")

    # Call C function
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
        data_size_proc,
        nb_blocks_loc,
        local_blocks_sizes,
        nnz,
        pixels,
        pixweights,
        signal,
        noise,
        params["Lambda"],
        invtt,
    )

    return

def MTmap(
    comm,
    params,
    sweeptstamps,
    nsweeps,
    az,
    az_min,
    az_max,
    hwp_angle,
    nces,
    data_size_proc,
    nb_blocks_loc,
    local_blocks_sizes,
    nnz,
    pixels,
    pixweights,
    signal,
    noise,
    invtt
):
    """
    Compute a map estimator through the marginalization over a set of time-domain templates.
    
    Parameters
    ----------
    * `comm`: Communicator over which data is distributed
    * `params`: Parameter dictionary
    * `sweeptstamps`: Index stamps indicating the borders of the compact supports of polynomial templates
    * `nsweeps`: Number of polynomial supports for each local CES
    * `az`: Boresight azimuth 
    * `az_min`: Minimum azimuth for each local CES
    * `az_max`: Maximum azimuth for each local CES
    * `hwp_angle`: HWP angle
    * `nces`: Number of local CES
    * `data_size_proc`: Data sizes in full communicator
    * `nb_blocks_loc`: Number of local observations
    * `local_blocks_sizes`: Local data sizes
    * `nnz`: Number of non-zero elements per row
    * `pixels`: Pixel indices of non-zero values
    * `pixweights`: Corresponding matrix values
    * `signal`: Signal buffer
    * `noise`: Noise buffer
    * `invtt`: Inverse noise weights

    """
    if not available:
        raise RuntimeError("No libmappraiser available, cannot reconstruct the map")
    
    outpath = params["path_output"].encode('ascii')
    ref = params["ref"].encode('ascii')
    
    fsamp = params["fsample"]
    # remove unit of fsamp
    try:
        f_unit = fsamp.unit
        fsamp = float(fsamp / (1.0 * f_unit))
    except AttributeError:
        pass
    
    # Format 2D arrays as arrays of uintp (pointers: void*)
    sweeptstamps_p = None
    az_p = None
    hwp_angle_p = None
    if sweeptstamps is not None:
        sweeptstamps_p = (sweeptstamps.__array_interface__['data'][0]
        + np.arange(sweeptstamps.shape[0])*sweeptstamps.strides[0]).astype(np.uintp)
    if az is not None:
        az_p = (az.__array_interface__['data'][0] + np.arange(az.shape[0])*az.strides[0]).astype(np.uintp)
    if hwp_angle is not None:
        hwp_angle_p = (hwp_angle.__array_interface__['data'][0] + np.arange(hwp_angle.shape[0])*hwp_angle.strides[0]).astype(np.uintp)
    
    # Call C-function
    _mappraiser.MTmap(
        encode_comm(comm),
        outpath,
        ref,
        params["solver"],
        params["precond"],
        params["Z_2lvl"],
        params["ptcomm_flag"],
        params["npoly"],
        params["nhwp"],
        params["hwpssbaseline_length"],
        params["sss"],
        params["sbins"],
        params["tol"],
        params["maxiter"],
        params["enlFac"],
        params["ortho_alg"],
        params["bs_red"],
        params["nside"],
        sweeptstamps_p,
        nsweeps,
        az_p,
        az_min,
        az_max,
        hwp_angle_p,
        nces,
        data_size_proc,
        nb_blocks_loc,
        local_blocks_sizes,
        nnz,
        pixels,
        pixweights,
        signal,
        noise,
        fsamp,
        invtt
    )
    
    return