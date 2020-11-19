# This script contains a list of routines to set mappraiser parameters, and
# apply the OpMappraiser operator during a TOD2MAP TOAST pipeline

#@author: Hamza El Bouhargani
#@date: January 2020

import argparse
import copy
import os
import re

import numpy as np
from toast.timing import function_timer, Timer
from toast.utils import Logger, Environment

from TOAST_interface import OpMappraiser

def add_mappraiser_args(parser):
    """ Add libmappraiser arguments
    """

    parser.add_argument(
        "--outpath", required=False, default="./", help="Output path"
    )
    parser.add_argument("--ref", required=False, default="run0", help="Output maps references"
    )
    parser.add_argument("--map-maker",dest="map_maker", required=False, default="ML", help="Map-making procedure"
    )
    parser.add_argument("--npoly", required=False, default=0, type=np.int, help="Order of polynomial templates"
    )
    parser.add_argument("--nhwp", required=False, default=0, type=np.int, help="Order of HWPSS templates"
    )
    parser.add_argument("--hwpss-base", required=False, default=10.0, type=np.double, help="HWPSS baseline length in seconds (default = 10.0s)"
    )
    parser.add_argument("--sss", required=False, default=0, type=np.int, help="SSS template on/off (only effective in MT mapper): 0-> off, 1->on"
    )
    parser.add_argument("--sbins", required=False, default=20, type=np.int, help="Number of azimuth bins in the SSS template"
    )
    parser.add_argument("--Lambda", required=False, default=16384, type=np.int, help="Half bandwidth (lambda) of noise covariance"
    )
    parser.add_argument("--solver", required=False, default=0, type=np.int, help="Choose map-making solver: 0->PCG, 1->ECG"
    )
    parser.add_argument("--ptcomm_flag", required=False, default=6, type=np.int, help="Choose collective communication scheme"
    )
    parser.add_argument("--tol", required=False, default=1e-6, type=np.double, help="Tolerance parameter for convergence"
    )
    parser.add_argument("--maxiter", required=False, default=500, type=np.int, help="Maximum number of iterations in Mappraiser"
    )
    parser.add_argument("--enlFac", required=False, default=1, type=np.int, help="Enlargement factor for ECG"
    )
    parser.add_argument("--ortho_alg", required=False, default=1, type=np.int, help="Orthogonalization scheme for ECG. O:odir, 1:omin"
    )
    parser.add_argument("--bs_red", required=False, default = 0, type=np.int, help="Use dynamic search reduction"
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
    """ Create a Mappraiser parameter dictionary.

    Initialize the Mappraiser parameters from the command line arguments.

    """
    params = {}

    params["map-maker"] = args.map_maker
    params["npoly"] = args.npoly
    params["nhwp"] = args.nhwp
    params["hwpss-base"] = args.hwpss_base
    params["sss"] = args.sss
    params["sbins"] = args.sbins
    params["nside"] = args.nside
    params["Lambda"] = args.Lambda
    params["samplerate"] = args.sample_rate
    params["hwp_rpm"] = args.hwp_rpm
    params["output"] = args.outpath
    params["ref"] = args.ref
    params["solver"] = args.solver
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
    """ Use libmappraiser to run unbiased map-making

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
        params= params,
        purge=True,
        name=signalname,
        noise_name = noisename,
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
                timer.report_clear("Mapping {}_telescope_{}_time_{}".format(
                args.outpath,
                tele_name,
                time_name,
                ))

    if comm.comm_world is not None:
        comm.comm_world.barrier()
    total_timer.stop()
    if comm.world_rank == 0 and verbose:
        total_timer.report("Mappraiser total")

    return
