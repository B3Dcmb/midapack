import os
import re

import numpy as np
import traitlets
from astropy import units as u

from toast.mpi import use_mpi
from toast.observation import default_values as defaults
from toast.templates import Offset
from toast.timing import function_timer, Timer
from toast.traits import Bool, Dict, Instance, Int, Float, Unicode, trait_docs
from toast.utils import Environment, Logger, dtype_to_aligned
from toast.ops.delete import Delete
from toast.ops.mapmaker import MapMaker
from toast.ops.operator import Operator

import toml

# local imports
from .utils import (
    compute_invtt,
    compute_local_block_sizes,
    log_time_memory,
    restore_in_turns,
    restore_local,
    stage_in_turns,
    stage_local,
    stage_polymetadata,
    stage_azscan,
    stage_hwpangle,
)

mappraiser = None
if use_mpi:
    try:
        import mappraiser_wrapper as mappraiser
    except ImportError:
        mappraiser = None


def available():
    """(bool): True if libmappraiser is found in the library search path."""
    global mappraiser
    return (mappraiser is not None) and mappraiser.available


def madam_params_from_mapmaker(mapmaker):
    """Utility function that configures Madam to match the TOAST mapmaker"""

    if not isinstance(mapmaker, MapMaker):
        raise RuntimeError("Need an instance of MapMaker to configure from")

    destripe_pixels = mapmaker.binning.pixel_pointing
    map_pixels = mapmaker.map_binning.pixel_pointing

    params = {
        "nside_cross": destripe_pixels.nside,
        "nside_map": map_pixels.nside,
        "nside_submap": map_pixels.nside_submap,
        "path_output": mapmaker.output_dir,
        "write_hits": mapmaker.write_hits,
        "write_matrix": mapmaker.write_invcov,
        "write_wcov": mapmaker.write_cov,
        "write_mask": mapmaker.write_rcond,
        "info": 3,
        "iter_max": mapmaker.iter_max,
        "pixlim_cross": mapmaker.solve_rcond_threshold,
        "pixlim_map": mapmaker.map_rcond_threshold,
        "cglimit": mapmaker.convergence,
    }
    sync_type = mapmaker.map_binning.sync_type
    if sync_type == "allreduce":
        params["allreduce"] = True
    elif sync_type == "alltoallv":
        params["concatenate_messages"] = True
        params["reassign_submaps"] = True
    else:
        msg = f"Unknown sync_type: {sync_type}"
        raise RuntimeError(msg)

    # Destriping parameters

    for template in mapmaker.template_matrix.templates:
        if isinstance(template, Offset):
            baselines = template
            break
    else:
        baselines = None

    if baselines is None or not baselines.enabled:
        params.update(
            {
                "write_binmap": True,
                "write_map": False,
                "kfirst": False,
            }
        )
    else:
        params.update(
            {
                "write_binmap": False,
                "write_map": True,
                "kfilter": baselines.use_noise_prior,
                "kfirst": True,
                "base_first": baselines.step_time.to_value(u.s),
                "precond_width_min": baselines.precond_width,
                "precond_width_max": baselines.precond_width,
                "good_baseline_fraction": baselines.good_fraction,
            }
        )

    return params


@trait_docs
class Mappraiser(Operator):
    """Operator which passes data to libmappraiser for map-making."""

    # Class traits

    API = Int(0, help="Internal interface version for this operator")

    params = Dict(default_value={}, help="Parameters to pass to mappraiser")

    paramfile = Unicode(
        None, allow_none=True, help="Read mappraiser parameters from this file"
    )

    # N.B: timestamps are not currently used in MAPPRAISER.
    # However, that may change in the future.
    times = Unicode(defaults.times, help="Observation shared key for timestamps")

    det_data = Unicode(
        defaults.det_data, help="Observation detdata key for the timestream data"
    )

    noise_name = Unicode(
        "noise",
        allow_none=True,
        help="Observation detdata key for noise data (if None, triggers noiseless mode)",
    )

    az_name = Unicode(
        "azimuth",
        allow_none=True,
        help="Observation boresight azimuth key"
    )

    hwpangle_name = Unicode(
        "hwp_angle",
        allow_none=True,
        help="Observation hwp angle key"
    )

    noiseless = Bool(
        False,
        help="Activate noiseless mode"
    )

    det_flags = Unicode(
        defaults.det_flags,
        allow_none=True,
        help="Observation detdata key for flags to use",
    )

    det_flag_mask = Int(
        defaults.det_mask_invalid, help="Bit mask value for optional detector flagging"
    )

    shared_flags = Unicode(
        defaults.shared_flags,
        allow_none=True,
        help="Observation shared key for telescope flags to use",
    )

    shared_flag_mask = Int(
        defaults.shared_mask_invalid,
        help="Bit mask value for optional shared flagging",
    )

    # instance of PixelsHealpix, operator which generates healpix pixel numbers
    pixel_pointing = Instance(
        klass=Operator,
        allow_none=True,
        help="This must be an instance of a pixel pointing operator",
    )

    # instance of StokesWeights, operator which generates I/Q/U pointing weights
    stokes_weights = Instance(
        klass=Operator,
        allow_none=True,
        help="This must be an instance of a Stokes weights operator",
    )

    # ! N.B: this is not supported for the moment.
    view = Unicode(
        None, allow_none=True, help="Use this view of the data in all observations"
    )

    # There is not tod cleaning in MAPPRAISER.
    # det_out = Unicode(
    #     None,
    #     allow_none=True,
    #     help="Observation detdata key for output destriped timestreams",
    # )

    noise_model = Unicode(
        "noise_model", help="Observation key containing the noise model"
    )

    purge_det_data = Bool(
        False,
        help="If True, clear all observation detector data after copying to mappraiser buffers",
    )

    restore_det_data = Bool(
        False,
        help="If True, restore detector data to observations on completion",
    )

    # ! N.B: not supported for the moment
    mcmode = Bool(
        False,
        help="If true, Madam will store auxiliary information such as pixel matrices and noise filter",
    )

    copy_groups = Int(
        1,
        help="The processes on each node are split into this number of groups to copy data in turns",
    )

    translate_timestamps = Bool(
        False, help="Translate timestamps to enforce monotonity"
    )

    noise_scale = Unicode(
        "noise_scale",
        help="Observation key with optional scaling factor for noise PSDs",
    )

    mem_report = Bool(
        False, help="Print system memory use while staging/unstaging data"
    )

    # Some traits for noise estimation

    save_psd = Bool(False, help="Save noise PSD information during inv_tt computation")

    apod_window_type = Unicode(
        "chebwin", help="Type of apodisation window to use during noise PSD estimation"
    )

    nperseg = Int(
        0,
        help="If 0, set nperseg = timestream length to compute the noise periodograms. If > 0, nperseg = Lambda.",
    )

    bandwidth = Int(16384, help="Half-bandwidth for the noise model")

    # Some useful traits for debugging

    noiseless = Bool(False, help="Activate noiseless mode")

    fill_noise_zero = Bool(
        False, help="Fill the noise vector with zeros just before calling Mappraiser"
    )

    downscale = Int(
        1,
        help="Scale down the noise by the sqrt of this number to artifically increase S/N ratio",
    )

    signal_fraction = Float(1.0, help="Fraction of the sky signal to keep")

    limit_det = Int(
        None,
        allow_none=True,
        help="Limit the number of local detectors to this number.",
    )

    # Additional parameters for the C library

    # map-making method
    map_maker = Unicode("ML",
        allow_none=False,
        help="Map-making approach: 'ML' for maximum likelihood (downweigthing approach), 'MT' for marginalized templates.",
    )

    # solver
    solver = Int(0, help="Choose mapmaking solver (0->PCG, 1->ECG)")
    tol = Float(1e-6, help="Convergence threshold for the iterative solver")
    maxiter = Int(3000, help="Maximum number of iterations allowed for the solver")
    enlFac = Int(1, help="Enlargement factor for ECG")
    bs_red = Int(0, help="Use dynamic search reduction")

    # preconditioner
    precond = Int(
        0, help="Choose preconditioner (0->BJ, 1->2lvl a priori, 2->2lvl a posteriori"
    )
    z_2lvl = Int(0, help="Size of 2lvl deflation space")
    ortho_alg = Int(1, help="Orthogonalization scheme for ECG (O->odir, 1->omin)")

    # communication algorithm
    ptcomm_flag = Int(6, help="Choose collective communication scheme")

    realization = Int(0, help="Noise realization index (for gap filling)")
    
    # Template relevant additions
    npoly = Int(0, help="Order of polynomial templates")
    fixed_polybase = Int(0, help="Flag to set a specific polynomial baseline length. When False the baseline length is set to the sweep duration (ground only).")
    polybaseline_length = Float(10.0, help="Polynomial baseline length in seconds (default = 10.0s)")
    nhwp = Int(0, help="Order of HWPSS templates")
    hwpssbaseline_length = Float(10.0, help="HWPSS baseline length in seconds (default = 10.0s)")
    sss = Int(0, help="Scan Synchronous Signal template flag.")
    sbins = Int(0, help="Number of azimuth bins in the SSS template")
    uniform_w = Int(0, help="Flag to use a uniform white noise model in reduction.")

    @traitlets.validate("shared_flag_mask")
    def _check_shared_flag_mask(self, proposal):
        check = proposal["value"]
        if check < 0:
            raise traitlets.TraitError("Shared flag mask should be a positive integer")
        return check

    @traitlets.validate("det_flag_mask")
    def _check_det_flag_mask(self, proposal):
        check = proposal["value"]
        if check < 0:
            raise traitlets.TraitError("Det flag mask should be a positive integer")
        return check

    @traitlets.validate("restore_det_data")
    def _check_restore_det_data(self, proposal):
        check = proposal["value"]
        if check and not self.purge_det_data:
            raise traitlets.TraitError(
                "Cannot set restore_det_data since purge_det_data is False"
            )
        # if check and self.det_out is not None:
        #     raise traitlets.TraitError(
        #         "Cannot set restore_det_data since det_out is not None"
        #     )
        return check

    # @traitlets.validate("det_out")
    # def _check_det_out(self, proposal):
    #     check = proposal["value"]
    #     if check is not None and self.restore_det_data:
    #         raise traitlets.TraitError(
    #             "If det_out is not None, restore_det_data should be False"
    #         )
    #     return check

    @traitlets.validate("pixel_pointing")
    def _check_pixel_pointing(self, proposal):
        pixels = proposal["value"]
        if pixels is not None:
            if not isinstance(pixels, Operator):
                raise traitlets.TraitError(
                    "pixel_pointing should be an Operator instance"
                )
            # Check that this operator has the traits we expect
            for trt in ["pixels", "create_dist", "view"]:
                if not pixels.has_trait(trt):
                    msg = f"pixel_pointing operator should have a '{trt}' trait"
                    raise traitlets.TraitError(msg)
        return pixels

    @traitlets.validate("stokes_weights")
    def _check_stokes_weights(self, proposal):
        weights = proposal["value"]
        if weights is not None:
            if not isinstance(weights, Operator):
                raise traitlets.TraitError(
                    "stokes_weights should be an Operator instance"
                )
            # Check that this operator has the traits we expect
            for trt in ["weights", "view"]:
                if not weights.has_trait(trt):
                    msg = f"stokes_weights operator should have a '{trt}' trait"
                    raise traitlets.TraitError(msg)
        return weights

    @traitlets.validate("params")
    def _check_params(self, proposal):
        check = proposal["value"]
        if "info" not in check:
            # The user did not specify the info level- set it from the toast loglevel
            env = Environment.get()
            level = env.log_level()
            if level == "DEBUG":
                check["info"] = 2
            elif level == "VERBOSE":
                check["info"] = 3
            else:
                check["info"] = 1
        return check

    # Check the traits that are not yet supported
    @traitlets.validate("view")
    def _check_view(self, proposal):
        check = proposal["value"]
        if check is not None:
            raise traitlets.TraitError("Views of the data are currently not supported")
        return check

    @traitlets.validate("mcmode")
    def _check_mcmode(self, proposal):
        check = proposal["value"]
        if check:
            raise traitlets.TraitError("MC mode is not currently supported")
        return check

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._cached = False
        self._logprefix = "Mappraiser:"

    def clear(self):
        """Delete the underlying memory.
        This will forcibly delete the C-allocated memory and invalidate all python
        references to the buffers.
        """

        if self._cached:
            # MAPPRAISER does not have caches to clear (no 'cached' mode implemented)
            # madam.clear_caches()
            self._cached = False

        for atr in [
            "timestamps",
            "signal",
            "blocksizes",
            "noise",
            "invtt",
            "pixels",
            "pixweights",
        ]:
            atrname = "_mappraiser_{}".format(atr)
            rawname = "{}_raw".format(atrname)
            if hasattr(self, atrname):
                delattr(self, atrname)
                raw = getattr(self, rawname)
                if raw is not None:
                    raw.clear()
                setattr(self, rawname, None)
                setattr(self, atrname, None)

    def __del__(self):
        self.clear()

    @function_timer
    def _exec(self, data, detectors=None, **kwargs):
        log = Logger.get()
        timer = Timer()
        timer.start()

        if not available():
            raise RuntimeError("Mappraiser is either not installed or MPI is disabled")

        if len(data.obs) == 0:
            raise RuntimeError(
                "Mappraiser requires every supplied data object to "
                "contain at least one observation"
            )

        for trait in "det_data", "pixel_pointing", "stokes_weights":
            if getattr(self, trait) is None:
                msg = f"You must set the '{trait}' trait before calling exec()"
                raise RuntimeError(msg)

        # Combine parameters from an external file and other parameters passed in

        params = dict()

        # repeat_keys = ["detset", "detset_nopol", "survey"]

        if self.paramfile is not None:
            if data.comm.world_rank == 0:
                line_pat = re.compile(r"(\S+)\s+=\s+(\S+)")
                comment_pat = re.compile(r"^\s*\#.*")
                with open(self.paramfile, "r") as f:
                    for line in f:
                        if comment_pat.match(line) is None:
                            line_mat = line_pat.match(line)
                            if line_mat is not None:
                                k = line_mat.group(1)
                                v = line_mat.group(2)
                                # if k in repeat_keys:
                                #     if k not in params:
                                #         params[k] = [v]
                                #     else:
                                #         params[k].append(v)
                                # else:
                                #     params[k] = v
                                params[k] = v
            if data.comm.comm_world is not None:
                params = data.comm.comm_world.bcast(params, root=0)
            for k, v in self.params.items():
                # if k in repeat_keys:
                #     if k not in params:
                #         params[k] = [v]
                #     else:
                #         params[k].append(v)
                # else:
                #     params[k] = v
                params[k] = v

        params.update({
            "Lambda": self.bandwidth,
            "solver": self.solver,
            "precond": self.precond,
            "Z_2lvl": self.z_2lvl,
            "ptcomm_flag": self.ptcomm_flag,
            "tol": np.double(self.tol),
            "maxiter": self.maxiter,
            "enlFac": self.enlFac,
            "ortho_alg": self.ortho_alg,
            "bs_red": self.bs_red,
            "map_maker": self.map_maker,
            "npoly": self.npoly,
            "fixed_polybase": self.fixed_polybase,
            "polybaseline_length": np.double(self.polybaseline_length),
            "nhwp": self.nhwp,
            "hwpssbaseline_length": np.double(self.hwpssbaseline_length),
            "sss": self.sss,
            "sbins": self.sbins,
            "uniform_w": self.uniform_w,
        })

        # params dictionary overrides operator traits for mappraiser C library arguments
        if self.params is not None:
            params.update(self.params)

        for key in "path_output", "ref":
            if key not in params:
                msg = f"You must set the '{key}' key of mappraiser.params before calling exec()"
                raise RuntimeError(msg)

        if "fsample" not in params:
            params["fsample"] = data.obs[0].telescope.focalplane.sample_rate.to_value(
                u.Hz
            )

        # Set mappraiser parameters that depend on our traits
        if self.mcmode:
            params["mcmode"] = True
        else:
            params["mcmode"] = False

        # if self.det_out is not None:
        #     params["write_tod"] = True
        # else:
        #     params["write_tod"] = False

        # Check if noiseless mode is activated
        if self.noiseless or (self.noise_name is None):
            self.noiseless = True
            self.noise_name = None
            # diagonal noise covariance
            params["Lambda"] = 1

        # Set half bandwidth of noise model to 1 in uniform weighting mode or MT map-making
        if params["uniform_w"] == 1 or params["map_maker"] == "MT":
            params["Lambda"] = 1

        # Log the libmappraiser parameters that were used.
        if data.comm.world_rank == 0:
            with open(
                    os.path.join(params["path_output"], "mappraiser_args_log.toml"),
                    "w",
            ) as f:
                toml.dump(params, f)

        # Check input parameters and compute the sizes of Mappraiser data objects
        if data.comm.world_rank == 0:
            msg = "{} Computing data sizes".format(self._logprefix)
            log.info(msg)
        (
            all_dets,
            nsamp,
            nnz,
            nnz_full,
            nnz_stride,
            interval_starts,
            psd_freqs,
        ) = self._prepare(params, data, detectors)

        log.info_rank(
            f"{self._logprefix} Parsed parameters in",
            comm=data.comm.comm_world,
            timer=timer,
        )

        # Stage data
        if data.comm.world_rank == 0:
            msg = "{} Copying toast data to buffers".format(self._logprefix)
            log.info(msg)

        (
            signal_dtype,
            data_size_proc,
            sweeptstamps,
            nsweeps,
            az,
            az_min,
            az_max,
            hwp_angle,
            nces,
        ) = self._stage_data(
            params,
            data,
            all_dets,
            nsamp,
            nnz,
            nnz_full,
            nnz_stride,
            interval_starts,
            psd_freqs,
        )

        log.info_rank(
            f"{self._logprefix} Staged data in",
            comm=data.comm.comm_world,
            timer=timer,
        )

        # Compute the map
        if params["map_maker"] == "ML":
            if data.comm.world_rank == 0:
                msg = "{} Computing the ML (LF downweighting) map".format(self._logprefix)
                log.info(msg)
            self._MLmap(
                params,
                data,
                data_size_proc,
                len(data.obs) * len(all_dets),
                nnz,
            )
        elif params["map_maker"] == "MT":
            if data.comm.world_rank == 0:
                msg = "{} Computing the MT (marginalized templates) map".format(self._logprefix)
                log.info(msg)
            self._MTmap(
                params,
                data,
                data_size_proc,
                len(data.obs) * len(all_dets),
                nnz,
                sweeptstamps,
                nsweeps,
                az,
                az_min,
                az_max,
                hwp_angle,
                nces,
            )
        else:
            raise RuntimeError(
                "Unvalid Map-making approach please choose:"
                " -'ML' for Maximum Likelihood map-making (or low frequency down-weighting)"
                " -'MT' for Marginalized Templates map-making"
            )

        log.info_rank(
            f"{self._logprefix} Processed time data in",
            comm=data.comm.comm_world,
            timer=timer,
        )

        # Unstage data
        if data.comm.world_rank == 0:
            msg = "{} Copying buffers back to toast data".format(self._logprefix)
            log.info(msg)
        self._unstage_data(
            params,
            data,
            all_dets,
            nsamp,
            nnz,
            nnz_full,
            interval_starts,
            signal_dtype,
        )

        log.info_rank(
            f"{self._logprefix} Unstaged data in",
            comm=data.comm.comm_world,
            timer=timer,
        )

        return

    def _finalize(self, data, **kwargs):
        return

    def _requires(self):
        req = {
            "meta": [self.noise_model],
            "shared": [
                self.times,
            ],
            "detdata": [self.det_data],
            "intervals": list(),
        }
        if self.view is not None:
            req["intervals"].append(self.view)
        if self.shared_flags is not None:
            req["shared"].append(self.shared_flags)
        if self.det_flags is not None:
            req["detdata"].append(self.det_flags)
        return req

    def _provides(self):
        prov = {"detdata": list()}
        # if self.det_out is not None:
        #     prov["detdata"].append(self.det_out)
        return prov

    @function_timer
    def _prepare(self, params, data, detectors):
        """Examine the data and determine quantities needed to set up Mappraiser buffers"""
        log = Logger.get()
        timer = Timer()
        timer.start()

        params["nside"] = self.pixel_pointing.nside

        # MAPPRAISER requires a fixed set of detectors and pointing matrix non-zeros.
        # Here we find the superset of local detectors used, and also the number
        # of pointing matrix elements.

        nsamp = 0

        # MAPPRAISER uses monolithic data buffers and specifies contiguous data intervals
        # in that buffer.  The starting sample index is used to mark the transition
        # between data intervals.
        interval_starts = list()

        # This quantity is only used for printing the fraction of samples in valid
        # ranges specified by the View.  Only samples actually in the view are copied
        # to Mappraiser buffers.
        # N.B: For the moment this is useless since MAPPRAISER does not use data views
        nsamp_valid = 0

        all_dets = set()
        nnz_full = None
        psd_freqs = None

        for ob in data.obs:
            # Get the detectors we are using for this observation
            dets = ob.select_local_detectors(detectors)

            if self.limit_det is not None:
                # cut all detectors but one
                dets = set(list(dets)[: self.limit_det])

            all_dets.update(dets)

            # Check that the timestamps exist.
            if self.times not in ob.shared:
                msg = (
                    "Shared timestamps '{}' does not exist in observation '{}'".format(
                        self.times, ob.name
                    )
                )
                raise RuntimeError(msg)

            # Check that the detector data and pointing exists in the observation
            if self.det_data not in ob.detdata:
                msg = "Detector data '{}' does not exist in observation '{}'".format(
                    self.det_data, ob.name
                )
                raise RuntimeError(msg)

            # Check that the noise model exists, and that the PSD frequencies are the
            # same across all observations (required by Mappraiser).
            if not self.noiseless:
                if self.noise_model not in ob:
                    msg = "Noise model '{}' not in observation '{}'".format(
                        self.noise_model, ob.name
                    )
                    raise RuntimeError(msg)
                if psd_freqs is None:
                    psd_freqs = np.array(
                        ob[self.noise_model].freq(ob.local_detectors[0]).to_value(u.Hz),
                        dtype=np.float64,
                    )
                else:
                    check_freqs = (
                        ob[self.noise_model].freq(ob.local_detectors[0]).to_value(u.Hz)
                    )
                    if not np.allclose(psd_freqs, check_freqs):
                        raise RuntimeError(
                            "All PSDs passed to Mappraiser must have the same frequency binning."
                        )

            # Are we using a view of the data?  If so, we will only be copying data in
            # those valid intervals.
            if self.view is not None:
                if self.view not in ob.intervals:
                    msg = "View '{}' does not exist in observation {}".format(
                        self.view, ob.name
                    )
                    raise RuntimeError(msg)
                # Go through all the intervals that will be used for our data view
                # and accumulate the number of samples.
                for intvw in ob.intervals[self.view]:
                    interval_starts.append(nsamp_valid)
                    nsamp_valid += intvw.last - intvw.first + 1
            else:
                interval_starts.append(nsamp_valid)
                nsamp_valid += ob.n_local_samples
            nsamp += ob.n_local_samples

        if data.comm.world_rank == 0:
            log.info(
                "{} {:.2f} % of samples are included in valid intervals.".format(
                    self._logprefix, nsamp_valid * 100.0 / nsamp
                )
            )

        nsamp = nsamp_valid

        interval_starts = np.array(interval_starts, dtype=np.int64)
        all_dets = sorted(all_dets)

        nnz_full = len(self.stokes_weights.mode)
        nnz_stride = None  # N.B: not used, useful for temperature-only maps

        if nnz_full != 3:
            msg = f"Mappraiser cannot make maps with nnz = {nnz_full}"
            raise RuntimeError(msg)
        else:
            nnz = nnz_full
            nnz_stride = 1

        if data.comm.world_rank == 0 and "path_output" in params:
            os.makedirs(params["path_output"], exist_ok=True)
            if self.save_psd:
                psdout = os.path.join(params["path_output"], "psd")
                os.makedirs(psdout, exist_ok=True)

        data.comm.comm_world.barrier()
        timer.stop()
        if data.comm.world_rank == 0:
            msg = "{} Compute data dimensions: {:0.1f} s".format(
                self._logprefix, timer.seconds()
            )
            log.debug(msg)

        return (
            all_dets,
            nsamp,
            nnz,
            nnz_full,
            nnz_stride,
            interval_starts,
            psd_freqs,
        )

    @function_timer
    def _stage_data(
            self,
            params,
            data,
            all_dets,
            nsamp,
            nnz,
            nnz_full,
            nnz_stride,
            interval_starts,
            psd_freqs,
    ):
        """Create mappraiser-compatible buffers.
        Collect the data into Mappraiser buffers.  If we are purging TOAST data to save
        memory, then optionally limit the number of processes that are copying at once.
        """
        log = Logger.get()
        timer = Timer()

        nodecomm = data.comm.comm_group_node

        # Determine how many processes per node should copy at once.
        n_copy_groups = 1
        if self.purge_det_data:
            # We will be purging some data- see if we should reduce the number of
            # processes copying in parallel (if we are not purging data, there
            # is no benefit to staggering the copy).
            if self.copy_groups > 0:
                n_copy_groups = min(self.copy_groups, nodecomm.size)

        if not self._cached:
            # Only do this if we have not cached the data yet.
            log_time_memory(
                data,
                prefix=self._logprefix,
                mem_msg="Before staging",
                full_mem=self.mem_report,
            )

        # Copy timestamps and PSDs all at once, since they are never purged.

        psds = dict()

        timer.start()

        # if not self._cached:
        #     timestamp_storage, _ = dtype_to_aligned(mappraiser.TIMESTAMP_TYPE)
        #     self._mappraiser_timestamps_raw = timestamp_storage.zeros(nsamp)
        #     self._mappraiser_timestamps = self_mappraiser_timestamps_raw.array()

        #     interval = 0
        #     time_offset = 0.0

        #     for ob in data.obs:
        #         for vw in ob.view[self.view].shared[self.times]:
        #             offset = interval_starts[interval]
        #             slc = slice(offset, offset + len(vw), 1)
        #             self._mappraiser_timestamps[slc] = vw
        #             if self.translate_timestamps:
        #                 off = self._mappraiser_timestamps[offset] - time_offset
        #                 self._mappraiser_timestamps[slc] -= off
        #                 time_offset = self._mappraiser_timestamps[slc][-1] + 1.0
        #             interval += 1

        #         # Get the noise object for this observation and create new
        #         # entries in the dictionary when the PSD actually changes.  The detector
        #         # weights are obtained from the noise model.

        #         nse = ob[self.noise_model]
        #         nse_scale = 1.0
        #         if self.noise_scale is not None:
        #             if self.noise_scale in ob:
        #                 nse_scale = float(ob[self.noise_scale])

        #         for det in all_dets:
        #             if det not in ob.local_detectors:
        #                 continue
        #             psd = nse.psd(det).to_value(u.K**2 * u.second) * nse_scale**2
        #             detw = nse.detector_weight(det)
        #             if det not in psds:
        #                 psds[det] = [(0.0, psd, detw)]
        #             else:
        #                 if not np.allclose(psds[det][-1][1], psd):
        #                     psds[det] += [(ob.shared[self.times][0], psd, detw)]

        #     log_time_memory(
        #         data,
        #         timer=timer,
        #         timer_msg="Copy timestamps and PSDs",
        #         prefix=self._logprefix,
        #         mem_msg="After timestamp staging",
        #         full_mem=self.mem_report,
        #     )

        # Copy the signal.  We always need to do this, even if we are running MCs.

        signal_dtype = data.obs[0].detdata[self.det_data].dtype

        if self._cached:
            # We have previously created the mappraiser buffers.  We just need to fill
            # them from the toast data.  Since both already exist we just copy the
            # contents.
            stage_local(
                data,
                nsamp,
                self.view,
                all_dets,
                self.det_data,
                self._mappraiser_signal,
                interval_starts,
                1,
                1,
                None,
                None,
                None,
                None,
                do_purge=False,
            )
        else:
            # Signal buffers do not yet exist
            if self.purge_det_data:
                # Allocate in a staggered way.
                self._mappraiser_signal_raw, self._mappraiser_signal = stage_in_turns(
                    data,
                    nodecomm,
                    n_copy_groups,
                    nsamp,
                    self.view,
                    all_dets,
                    self.det_data,
                    mappraiser.SIGNAL_TYPE,
                    interval_starts,
                    1,
                    1,
                    None,
                    None,
                    None,
                    None,
                )
            else:
                # Allocate and copy all at once.
                storage, _ = dtype_to_aligned(mappraiser.SIGNAL_TYPE)
                self._mappraiser_signal_raw = storage.zeros(nsamp * len(all_dets))
                self._mappraiser_signal = self._mappraiser_signal_raw.array()

                stage_local(
                    data,
                    nsamp,
                    self.view,
                    all_dets,
                    self.det_data,
                    self._mappraiser_signal,
                    interval_starts,
                    1,
                    1,
                    None,
                    None,
                    None,
                    None,
                    do_purge=False,
                )
            # Create buffer for local_block_sizes
            b_storage, _ = dtype_to_aligned(mappraiser.PIXEL_TYPE)
            self._mappraiser_blocksizes_raw = b_storage.zeros(
                len(data.obs) * len(all_dets)
            )
            self._mappraiser_blocksizes = self._mappraiser_blocksizes_raw.array()

        # Compute sizes of local data blocks
        compute_local_block_sizes(
            data,
            self.view,
            all_dets,
            self._mappraiser_blocksizes,
        )

        # Gather data sizes of the full communicator in global array
        data_size_proc = np.array(
            data.comm.comm_world.allgather(len(self._mappraiser_signal)), dtype=np.int32
        )

        # debug purposes
        if self.signal_fraction < 1.0:
            self._mappraiser_signal *= self.signal_fraction

        log_time_memory(
            data,
            timer=timer,
            timer_msg="Copy signal",
            prefix=self._logprefix,
            mem_msg="After signal staging",
            full_mem=self.mem_report,
        )

        # Copy the noise.
        # For the moment, in the absence of a gap-filling procedure in MAPPRAISER,
        # we separate signal and noise in the simulations

        # Check if the simulation contains any noise at all.
        if self.noiseless:
            msg = "{} Noiseless mode -> noise buffer filled with zeros".format(
                self._logprefix
            )
            log.info_rank(
                msg,
                data.comm.comm_world,
            )

        if self._cached:
            # We have previously created the mappraiser buffers.  We just need to fill
            # them from the toast data.  Since both already exist we just copy the
            # contents.
            stage_local(
                data,
                nsamp,
                self.view,
                all_dets,
                self.noise_name,
                self._mappraiser_noise,
                interval_starts,
                1,
                1,
                None,
                None,
                None,
                None,
                do_purge=False,
            )
        else:
            # Signal buffers do not yet exist
            if self.purge_det_data:
                # Allocate in a staggered way.
                self._mappraiser_noise_raw, self._mappraiser_noise = stage_in_turns(
                    data,
                    nodecomm,
                    n_copy_groups,
                    nsamp,
                    self.view,
                    all_dets,
                    self.noise_name,
                    mappraiser.SIGNAL_TYPE,
                    interval_starts,
                    1,
                    1,
                    None,
                    None,
                    None,
                    None,
                )
            else:
                # Allocate and copy all at once.
                storage, _ = dtype_to_aligned(mappraiser.SIGNAL_TYPE)
                self._mappraiser_noise_raw = storage.zeros(nsamp * len(all_dets))
                self._mappraiser_noise = self._mappraiser_noise_raw.array()

                stage_local(
                    data,
                    nsamp,
                    self.view,
                    all_dets,
                    self.noise_name,
                    self._mappraiser_noise,
                    interval_starts,
                    1,
                    1,
                    None,
                    None,
                    None,
                    None,
                    do_purge=False,
                )
            # Create buffer for invtt
            storage, _ = dtype_to_aligned(mappraiser.INVTT_TYPE)
            self._mappraiser_invtt_raw = storage.zeros(
                len(data.obs) * len(all_dets) * params["Lambda"]
            )
            self._mappraiser_invtt = self._mappraiser_invtt_raw.array()

        # Compute invtt
        if not self.noiseless:
            # Scale down the noise if we want
            if self.downscale > 1:
                self._mappraiser_noise /= np.sqrt(self.downscale)
            
            if params["uniform_w"] == 1:
                self._mappraiser_invtt[:] = 1.0
            else:
                compute_invtt(
                    len(data.obs),
                    len(all_dets),
                    self._mappraiser_noise,
                    self._mappraiser_blocksizes,
                    params["Lambda"],
                    params["fsample"],
                    self._mappraiser_invtt,
                    mappraiser.INVTT_TYPE,
                    apod_window_type=self.apod_window_type,
                    print_info=(data.comm.world_rank == 0),
                    save_psd=(self.save_psd and data.comm.world_rank == 0),
                    save_dir=os.path.join(params["path_output"], "psd"),
                )

            if self.fill_noise_zero:
                # just set the noise to zero
                # that way we can make the mapmaker iterate but on signal-only data
                self._mappraiser_noise[:] = 0.0
        else:
            self._mappraiser_noise[:] = 0.0
            self._mappraiser_invtt[:] = 1.0

        log_time_memory(
            data,
            timer=timer,
            timer_msg="Copy noise",
            prefix=self._logprefix,
            mem_msg="After noise staging",
            full_mem=self.mem_report,
        )

        # Copy the pointing

        nested_pointing = self.pixel_pointing.nest
        if not nested_pointing:
            # Any existing pixel numbers are in the wrong ordering
            Delete(detdata=[self.pixel_pointing.pixels]).apply(data)
            self.pixel_pointing.nest = True

        if not self._cached:
            # We do not have the pointing yet.
            self._mappraiser_pixels_raw, self._mappraiser_pixels = stage_in_turns(
                data,
                nodecomm,
                n_copy_groups,
                nsamp,
                self.view,
                all_dets,
                self.pixel_pointing.pixels,
                mappraiser.PIXEL_TYPE,
                interval_starts,
                nnz,
                1,
                self.shared_flags,
                self.shared_flag_mask,
                self.det_flags,
                self.det_flag_mask,
                operator=self.pixel_pointing,
                n_repeat=nnz,
            )

            # Arrange pixel indices for MAPPRAISER
            self._mappraiser_pixels *= nnz
            for i in range(nnz):
                self._mappraiser_pixels[i::nnz] += i

            (
                self._mappraiser_pixweights_raw,
                self._mappraiser_pixweights,
            ) = stage_in_turns(
                data,
                nodecomm,
                n_copy_groups,
                nsamp,
                self.view,
                all_dets,
                self.stokes_weights.weights,
                mappraiser.WEIGHT_TYPE,
                interval_starts,
                nnz,
                nnz_stride,
                None,
                None,
                None,
                None,
                operator=self.stokes_weights,
            )

            log_time_memory(
                data,
                timer=timer,
                timer_msg="Copy pointing",
                prefix=self._logprefix,
                mem_msg="After pointing staging",
                full_mem=self.mem_report,
            )

        if not nested_pointing:
            # Any existing pixel numbers are in the wrong ordering
            Delete(detdata=[self.pixel_pointing.pixels]).apply(data)
            self.pixel_pointing.nest = False
        
        # Stage instrument data to build templates in case of template marginalization
        sweeptstamps = None
        nsweeps = None
        az = None
        az_min = None
        az_max = None
        hwp_angle = None
        nces = len(data.obs)
        if params["map_maker"] == "MT":
            
            (
                sweeptstamps,
                nsweeps,
            ) = stage_polymetadata(
                data,
                self.view,
                params,
                self.shared_flags,
            )

            if params["sss"]:
                (
                    az,
                    az_min,
                    az_max,
                ) = stage_azscan(
                    data,
                    self.view,
                    self.az_name,
                )
            if params["nhwp"] != 0:
                hwp_angle = stage_hwpangle(
                    data,
                    self.view,
                    self.hwpangle_name,                    
                )

        # The following is basically useless for Mappraiser.

        # psdinfo = None

        # if not self._cached:
        #     # Detectors weights.  Madam assumes a single noise weight for each detector
        #     # that is constant.  We set this based on the first observation or else use
        #     # uniform weighting.

        #     ndet = len(all_dets)
        #     detweights = np.ones(ndet, dtype=np.float64)

        #     if len(psds) > 0:
        #         npsdbin = len(psd_freqs)
        #         npsd = np.zeros(ndet, dtype=np.int64)
        #         psdstarts = []
        #         psdvals = []
        #         for idet, det in enumerate(all_dets):
        #             if det not in psds:
        #                 raise RuntimeError("Every detector must have at least one PSD")
        #             psdlist = psds[det]
        #             npsd[idet] = len(psdlist)
        #             for psdstart, psd, detw in psdlist:
        #                 psdstarts.append(psdstart)
        #                 psdvals.append(psd)
        #             detweights[idet] = psdlist[0][2]
        #         npsdtot = np.sum(npsd)
        #         psdstarts = np.array(psdstarts, dtype=np.float64)
        #         psdvals = np.hstack(psdvals).astype(mappraiser.PSD_TYPE)
        #         npsdval = psdvals.size
        #     else:
        #         # Uniform weighting
        #         npsd = np.ones(ndet, dtype=np.int64)
        #         npsdtot = np.sum(npsd)
        #         psdstarts = np.zeros(npsdtot)
        #         npsdbin = 10
        #         fsample = 10.0
        #         psd_freqs = np.arange(npsdbin) * fsample / npsdbin
        #         npsdval = npsdbin * npsdtot
        #         psdvals = np.ones(npsdval)

        #     psdinfo = (detweights, npsd, psdstarts, psd_freqs, psdvals)

        #     log_time_memory(
        #         data,
        #         timer=timer,
        #         timer_msg="Collect PSD info",
        #         prefix=self._logprefix,
        #     )
        timer.stop()

        # return signal_dtype and info on size of distributed data
        return (
            signal_dtype,
            data_size_proc,
            sweeptstamps,
            nsweeps,
            az,
            az_min,
            az_max,
            hwp_angle,
            nces,
        )

    @function_timer
    def _unstage_data(
            self,
            params,
            data,
            all_dets,
            nsamp,
            nnz,
            nnz_full,
            interval_starts,
            signal_dtype,
    ):
        """
        Restore data to TOAST observations.
        Optionally copy the signal and pointing back to TOAST if we previously
        purged it to save memory.  Also copy the destriped timestreams if desired.
        """
        log = Logger.get()
        timer = Timer()

        nodecomm = data.comm.comm_group_node

        # Determine how many processes per node should copy at once.
        n_copy_groups = 1
        if self.purge_det_data:
            # We MAY be restoring some data- see if we should reduce the number of
            # processes copying in parallel (if we are not purging data, there
            # is no benefit to staggering the copy).
            if self.copy_groups > 0:
                n_copy_groups = min(self.copy_groups, nodecomm.size)

        log_time_memory(
            data,
            prefix=self._logprefix,
            mem_msg="Before un-staging",
            full_mem=self.mem_report,
        )

        # Copy the signal

        timer.start()

        out_name = self.det_data
        # if self.det_out is not None:
        #     out_name = self.det_out

        # if self.det_out is not None or (self.purge_det_data and self.restore_det_data):
        if self.purge_det_data and self.restore_det_data:
            # We are copying some kind of signal back
            if not self.mcmode:
                # We are not running multiple realizations, so delete as we copy.
                # Restore signal
                restore_in_turns(
                    data,
                    nodecomm,
                    n_copy_groups,
                    nsamp,
                    self.view,
                    all_dets,
                    out_name,
                    signal_dtype,
                    self._mappraiser_signal,
                    self._mappraiser_signal_raw,
                    interval_starts,
                    1,
                )
                del self._mappraiser_signal
                del self._mappraiser_signal_raw
                del self._mappraiser_blocksizes
                del self._mappraiser_blocksizes_raw
                # Restore noise
                restore_in_turns(
                    data,
                    nodecomm,
                    n_copy_groups,
                    nsamp,
                    self.view,
                    all_dets,
                    out_name,
                    signal_dtype,
                    self._mappraiser_noise,
                    self._mappraiser_noise_raw,
                    interval_starts,
                    1,
                )
                del self._mappraiser_noise
                del self._mappraiser_noise_raw
                del self._mappraiser_invtt
                del self._mappraiser_invtt_raw
            else:
                # We want to re-use the signal buffer, just copy.
                restore_local(
                    data,
                    nsamp,
                    self.view,
                    all_dets,
                    out_name,
                    signal_dtype,
                    self._mappraiser_signal,
                    interval_starts,
                    1,
                )
                restore_local(
                    data,
                    nsamp,
                    self.view,
                    all_dets,
                    out_name,
                    signal_dtype,
                    self._mappraiser_noise,
                    interval_starts,
                    1,
                )

            log_time_memory(
                data,
                timer=timer,
                timer_msg="Copy signal",
                prefix=self._logprefix,
                mem_msg="After restoring signal",
                full_mem=self.mem_report,
            )

        if not self.mcmode:
            # We can clear the cached pointing
            del self._mappraiser_pixels
            del self._mappraiser_pixels_raw
            del self._mappraiser_pixweights
            del self._mappraiser_pixweights_raw
        return

    @function_timer
    def _MLmap(self, params, data, data_size_proc, nb_blocks_loc, nnz):
        """Compute the ML map from buffered data."""
        log_time_memory(
            data,
            prefix=self._logprefix,
            mem_msg="Just before libmappraiser.MLmap",
            full_mem=self.mem_report,
        )

        # ? how to use mcmode (not yet supported)
        # if self._cached:
        # -> call mappraiser in "cached" mode
        # else:
        # -> call mappraiser in normal mode
        # -> if self.mcmode: self._cached=True

        mappraiser.MLmap(
            data.comm.comm_world,
            params,
            data_size_proc,
            nb_blocks_loc,
            self._mappraiser_blocksizes,
            nnz,
            self._mappraiser_pixels,
            self._mappraiser_pixweights,
            self._mappraiser_signal,
            self._mappraiser_noise,
            self._mappraiser_invtt,
        )

        return
    
    @function_timer
    def _MTmap(self, params, data, data_size_proc, nb_blocks_loc, nnz, 
        sweeptstamps, nsweeps, az, az_min, az_max, hwp_angle, nces):
        """Compute the MT map from buffered data."""
        log_time_memory(
            data,
            prefix=self._logprefix,
            mem_msg="Just before libmappraiser.MTmap",
            full_mem=self.mem_report,
        )

        # ? how to use mcmode (not yet supported)
        # if self._cached:
        # -> call mappraiser in "cached" mode
        # else:
        # -> call mappraiser in normal mode
        # -> if self.mcmode: self._cached=True
        
        # Force number of OMP threads to 1
        os.environ["OMP_NUM_THREADS"] = "1"

        mappraiser.MTmap(
            data.comm.comm_world,
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
            self._mappraiser_blocksizes,
            nnz,
            self._mappraiser_pixels,
            self._mappraiser_pixweights,
            self._mappraiser_signal,
            self._mappraiser_noise,
            self._mappraiser_invtt,
        )

        return

    def _requires(self):
        req = self.pixel_pointing.requires()
        req.update(self.stokes_weights.requires())
        req["meta"].extend([self.noise_model])
        req["detdata"].extend([self.det_data])
        if self.shared_flags is not None:
            req["shared"].append(self.shared_flags)
        if self.det_flags is not None:
            req["detdata"].append(self.det_flags)
        return req

    def _provides(self):
        return dict()
