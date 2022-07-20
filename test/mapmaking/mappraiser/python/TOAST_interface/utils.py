# This script contains a list of routines to set mappraiser parameters, and
# apply the Mappraiser operator during a TOD2MAP TOAST(3) pipeline

# @author: Hamza El Bouhargani
# @date: January 2020

import argparse
import copy
import os
import re

import numpy as np
from toast.timing import function_timer
from toast.utils import GlobalTimers, Logger, Timer, dtype_to_aligned, memreport
from toast.ops.memory_counter import MemoryCounter

from TOAST_interface import Mappraiser


def add_mappraiser_args(parser):
    """Add libmappraiser arguments"""

    parser.add_argument("--outpath", required=False, default="./", help="Output path")
    parser.add_argument(
        "--ref", required=False, default="run0", help="Output maps references"
    )
    parser.add_argument(
        "--Lambda",
        required=False,
        default=16384,
        type=np.int,
        help="Half bandwidth (lambda) of noise covariance",
    )
    parser.add_argument(
        "--uniform_w",
        required=False,
        default=0,
        type=np.int,
        help="Activate for uniform white noise model: 0->off, 1->on",
    )
    parser.add_argument(
        "--solver",
        required=False,
        default=0,
        type=np.int,
        help="Choose map-making solver: 0->PCG, 1->ECG",
    )
    parser.add_argument(
        "--precond",
        required=False,
        default=0,
        type=np.int,
        help="Choose map-making preconditioner: 0->BD, 1->2lvl a priori, 2->2lvl a posteriori",
    )
    parser.add_argument(
        "--Z_2lvl", required=False, default=0, type=np.int, help="2lvl deflation size"
    )
    parser.add_argument(
        "--ptcomm_flag",
        required=False,
        default=6,
        type=np.int,
        help="Choose collective communication scheme",
    )
    parser.add_argument(
        "--tol",
        required=False,
        default=1e-6,
        type=np.double,
        help="Tolerance parameter for convergence",
    )
    parser.add_argument(
        "--maxiter",
        required=False,
        default=500,
        type=np.int,
        help="Maximum number of iterations in Mappraiser",
    )
    parser.add_argument(
        "--enlFac",
        required=False,
        default=1,
        type=np.int,
        help="Enlargement factor for ECG",
    )
    parser.add_argument(
        "--ortho_alg",
        required=False,
        default=1,
        type=np.int,
        help="Orthogonalization scheme for ECG. O:odir, 1:omin",
    )
    parser.add_argument(
        "--bs_red",
        required=False,
        default=0,
        type=np.int,
        help="Use dynamic search reduction",
    )
    parser.add_argument(
        "--conserve-memory",
        dest="conserve_memory",
        required=False,
        action="store_true",
        help="Conserve memory when staging libMappraiser buffers [default]",
    )
    parser.add_argument(
        "--no-conserve-memory",
        dest="conserve_memory",
        required=False,
        action="store_false",
        help="Do not conserve memory when staging libMappraiser buffers",
    )
    parser.set_defaults(conserve_memory=True)

    # `nside` may already be added
    try:
        parser.add_argument(
            "--nside", required=False, default=512, type=np.int, help="Healpix NSIDE"
        )
    except argparse.ArgumentError:
        pass
    # Common flag mask may already be added
    try:
        parser.add_argument(
            "--common-flag-mask",
            required=False,
            default=1,
            type=np.uint8,
            help="Common flag mask",
        )
    except argparse.ArgumentError:
        pass
    # `sample-rate` may be already added
    try:
        parser.add_argument(
            "--sample-rate",
            required=False,
            default=100.0,
            type=np.float,
            help="Detector sample rate (Hz)",
        )
    except argparse.ArgumentError:
        pass
    return


@function_timer
def setup_mappraiser(args):
    """Create a Mappraiser parameter dictionary.

    Initialize the Mappraiser parameters from the command line arguments.

    """
    params = {}

    params["nside"] = args.nside
    params["Lambda"] = args.Lambda
    params["uniform_w"] = args.uniform_w
    params["samplerate"] = args.sample_rate
    params["output"] = args.outpath
    params["ref"] = args.ref
    params["solver"] = args.solver
    params["precond"] = args.precond
    params["Z_2lvl"] = args.Z_2lvl
    params["pointing_commflag"] = args.ptcomm_flag
    params["tol"] = args.tol
    params["maxiter"] = args.maxiter
    params["enlFac"] = args.enlFac
    params["ortho_alg"] = args.ortho_alg
    params["bs_red"] = args.bs_red

    return params


@function_timer
def apply_mappraiser(
    args,
    comm,
    data,
    params,
    signalname,
    noisename,
    time_comms=None,
    telescope_data=None,
    verbose=True,
):
    """Use libmappraiser to run the ML map-making

    Args:
        time_comms (iterable) :  Series of disjoint communicators that
            map, e.g., seasons and days.  Each entry is a tuple of
            the form (`name`, `communicator`)
        telescope_data (iterable) : series of disjoint TOAST data
            objects.  Each entry is tuple of the form (`name`, `data`).
    """
    if comm.comm_world is None:
        raise RuntimeError("Mappraiser requires MPI")

    log = Logger.get()
    total_timer = Timer()
    total_timer.start()
    if comm.world_rank == 0 and verbose:
        log.info("Making maps")

    mappraiser = OpMappraiser(
        params=params,
        purge=True,
        name=signalname,
        noise_name=noisename,
        conserve_memory=args.conserve_memory,
    )

    if time_comms is None:
        time_comms = [("all", comm.comm_world)]

    if telescope_data is None:
        telescope_data = [("all", data)]

    timer = Timer()
    for time_name, time_comm in time_comms:
        for tele_name, tele_data in telescope_data:
            if len(time_name.split("-")) == 3:
                # Special rules for daily maps
                if args.do_daymaps:
                    continue
                if len(telescope_data) > 1 and tele_name == "all":
                    # Skip daily maps over multiple telescopes
                    continue

            timer.start()
            # N.B: code below is for Madam but may be useful to copy in Mappraiser
            # once we start doing multiple maps in one run
            # madam.params["file_root"] = "{}_telescope_{}_time_{}".format(
            #     file_root, tele_name, time_name
            # )
            # if time_comm == comm.comm_world:
            #     madam.params["info"] = info
            # else:
            #     # Cannot have verbose output from concurrent mapmaking
            #     madam.params["info"] = 0
            # if (time_comm is None or time_comm.rank == 0) and verbose:
            #     log.info("Mapping {}".format(madam.params["file_root"]))
            mappraiser.exec(tele_data, time_comm)

            if time_comm is not None:
                time_comm.barrier()
            if comm.world_rank == 0 and verbose:
                timer.report_clear(
                    "Mapping {}_telescope_{}_time_{}".format(
                        args.outpath,
                        tele_name,
                        time_name,
                    )
                )

    if comm.comm_world is not None:
        comm.comm_world.barrier()
    total_timer.stop()
    if comm.world_rank == 0 and verbose:
        total_timer.report("Mappraiser total")

    return


# helper functions adapted from toast/src/ops/madam_utils.py


def log_time_memory(
    data, timer=None, timer_msg=None, mem_msg=None, full_mem=False, prefix=""
):
    log = Logger.get()
    data.comm.comm_world.barrier()
    restart = False

    if timer is not None:
        if timer.is_running():
            timer.stop()
            restart = True

        if data.comm.world_rank == 0:
            msg = "{} {}: {:0.1f} s".format(prefix, timer_msg, timer.seconds())
            log.debug(msg)

    if mem_msg is not None:
        # Dump toast memory use
        mem_count = MemoryCounter(silent=True)
        mem_count.total_bytes = 0
        toast_bytes = mem_count.apply(data)

        if data.comm.group_rank == 0:
            msg = "{} {} Group {} memory = {:0.2f} GB".format(
                prefix, mem_msg, data.comm.group, toast_bytes / 1024**2
            )
            log.debug(msg)
        if full_mem:
            _ = memreport(
                msg="{} {}".format(prefix, mem_msg), comm=data.comm.comm_world
            )
    if restart:
        timer.start()


def stage_local(
    data,
    nsamp,
    view,
    dets,
    detdata_name,
    mappraiser_buffer,
    interval_starts,
    nnz,
    nnz_stride,
    shared_flags,
    shared_mask,
    det_flags,
    det_mask,
    do_purge=False,
    operator=None,
    compute_blocksizes=False,
    b_buffer=None,
):
    """Helper function to fill a mappraiser buffer from a local detdata key."""
    n_det = len(dets)
    nb_obs_loc = len(data.obs)
    if compute_blocksizes:
        assert b_buffer is not None
    interval_offset = 0
    do_flags = False
    if shared_flags is not None or det_flags is not None:
        do_flags = True
        # Flagging should only be enabled when we are processing the pixel indices
        # (which is how madam effectively implements flagging).  So we will set
        # all flagged samples to "-1" below.
        # TODO: how do we handle flagging in mappraiser?
        if nnz != 1:
            raise RuntimeError(
                "Internal error on mappraiser copy.  Only pixel indices should be flagged."
            )
    for iobs, ob in enumerate(data.obs):
        views = ob.view[view]
        for idet, det in enumerate(dets):
            if det not in ob.local_detectors:
                continue
            if operator is not None:
                # Synthesize data for staging
                obs_data = data.select(obs_uid=ob.uid)
                operator.apply(obs_data, detectors=[det])
            # Loop over views
            for ivw, vw in enumerate(views):
                view_samples = None
                if vw.start is None:
                    # This is a view of the whole obs
                    view_samples = ob.n_local_samples
                else:
                    view_samples = vw.stop - vw.start
                if compute_blocksizes:
                    b_buffer[idet * nb_obs_loc + iobs] += view_samples
                offset = interval_starts[interval_offset + ivw]
                flags = None
                if do_flags:
                    # Using flags
                    flags = np.zeros(view_samples, dtype=np.uint8)
                if shared_flags is not None:
                    flags |= views.shared[shared_flags][ivw] & shared_mask

                slc = slice(
                    (idet * nsamp + offset) * nnz,
                    (idet * nsamp + offset + view_samples) * nnz,
                    1,
                )
                if nnz > 1:
                    mappraiser_buffer[slc] = views.detdata[detdata_name][ivw][
                        det
                    ].flatten()[::nnz_stride]
                else:
                    mappraiser_buffer[slc] = views.detdata[detdata_name][ivw][
                        det
                    ].flatten()
                detflags = None
                if do_flags:
                    if det_flags is None:
                        detflags = flags
                    else:
                        detflags = np.copy(flags)
                        detflags |= views.detdata[det_flags][ivw][det] & det_mask
                    mappraiser_buffer[slc][detflags != 0] = -1
        if do_purge:
            del ob.detdata[detdata_name]
        interval_offset += len(views)
    return


def stage_in_turns(
    data,
    nodecomm,
    n_copy_groups,
    nsamp,
    view,
    dets,
    detdata_name,
    mappraiser_dtype,
    interval_starts,
    nnz,
    nnz_stride,
    shared_flags,
    shared_mask,
    det_flags,
    det_mask,
    operator=None,
    compute_blocksizes=False,
):
    """When purging data, take turns staging it."""
    raw = None
    wrapped = None
    for copying in range(n_copy_groups):
        if nodecomm.rank % n_copy_groups == copying:
            # Our turn to copy data
            storage, _ = dtype_to_aligned(mappraiser_dtype)
            raw = storage.zeros(nsamp * len(dets) * nnz)
            wrapped = raw.array()
            if compute_blocksizes:
                # create buffer for local_block_sizes
                b_storage, _ = dtype_to_aligned(np.int32)
                # TODO: don't define this dtype explicitly
                b_raw = b_storage.zeros(len(data.obs) * len(dets))
                b_wrapped = b_raw.array()
            else:
                b_wrapped = None
            stage_local(
                data,
                nsamp,
                view,
                dets,
                detdata_name,
                wrapped,
                interval_starts,
                nnz,
                nnz_stride,
                shared_flags,
                shared_mask,
                det_flags,
                det_mask,
                do_purge=True,
                operator=operator,
                compute_blocksizes=compute_blocksizes,
                blocksizes_buffer=b_wrapped,
            )
        nodecomm.barrier()
    if compute_blocksizes:
        return raw, wrapped, b_raw, b_wrapped
    else:
        return raw, wrapped


def restore_local(
    data,
    nsamp,
    view,
    dets,
    detdata_name,
    detdata_dtype,
    mappraiser_buffer,
    interval_starts,
    nnz,
):
    """Helper function to create a detdata buffer from mappraiser data."""
    n_det = len(dets)
    interval = 0
    for ob in data.obs:
        # Create the detector data
        if nnz == 1:
            ob.detdata.create(detdata_name, dtype=detdata_dtype)
        else:
            ob.detdata.create(detdata_name, dtype=detdata_dtype, sample_shape=(nnz,))
        # Loop over views
        views = ob.view[view]
        for ivw, vw in enumerate(views):
            view_samples = None
            if vw.start is None:
                # This is a view of the whole obs
                view_samples = ob.n_local_samples
            else:
                view_samples = vw.stop - vw.start
            offset = interval_starts[interval]
            ldet = 0
            for det in dets:
                if det not in ob.local_detectors:
                    continue
                idet = ob.local_detectors.index(det)
                slc = slice(
                    (idet * nsamp + offset) * nnz,
                    (idet * nsamp + offset + view_samples) * nnz,
                    1,
                )
                if nnz > 1:
                    views.detdata[detdata_name][ivw][ldet] = mappraiser_buffer[
                        slc
                    ].reshape((-1, nnz))
                else:
                    views.detdata[detdata_name][ivw][ldet] = mappraiser_buffer[slc]
                ldet += 1
            interval += 1
    return


def restore_in_turns(
    data,
    nodecomm,
    n_copy_groups,
    nsamp,
    view,
    dets,
    detdata_name,
    detdata_dtype,
    mappraiser_buffer,
    mappraiser_buffer_raw,
    interval_starts,
    nnz,
):
    """When restoring data, take turns copying it."""
    for copying in range(n_copy_groups):
        if nodecomm.rank % n_copy_groups == copying:
            # Our turn to copy data
            restore_local(
                data,
                nsamp,
                view,
                dets,
                detdata_name,
                detdata_dtype,
                mappraiser_buffer,
                interval_starts,
                nnz,
            )
            mappraiser_buffer_raw.clear()
        nodecomm.barrier()
    return
