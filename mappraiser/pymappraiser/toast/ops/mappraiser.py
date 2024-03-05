import os
import re
import tomlkit

import numpy as np
import traitlets
from astropy import units as u

from toast.mpi import MPI, use_mpi
from toast.observation import default_values as defaults
from toast.timing import function_timer, Timer
from toast.traits import Bool, Dict, Instance, Int, Float, Unicode, trait_docs
from toast.utils import Environment, Logger, dtype_to_aligned
from toast.ops.delete import Delete
from toast.ops.operator import Operator
from toast.ops.arithmetic import Combine

from .mappraiser_utils import (
    compute_autocorrelations,
    log_time_memory,
    stage_in_turns,
    stage_local,
)

libmappraiser = None
if use_mpi:
    try:
        import pymappraiser.wrapper as libmappraiser
    except ImportError:
        libmappraiser = None


def available():
    """(bool): True if libmappraiser is found in the library search path."""
    global libmappraiser
    return (libmappraiser is not None) and libmappraiser.available


@trait_docs
class Mappraiser(Operator):
    """Operator which passes data to libmappraiser for map-making."""

    # Class traits

    API = Int(0, help="Internal interface version for this operator")

    params = Dict(default_value={}, help="Parameters to pass to mappraiser")

    paramfile = Unicode(
        None, allow_none=True, help="Read mappraiser parameters from this file"
    )

    noise_name = Unicode(
        "noise",
        allow_none=True,
        help="Observation detdata key for noise data (if None, triggers noiseless mode)",
    )

    atm_name = Unicode(
        None,
        allow_none=True,
        help="Observation detdata key for atm data to be added with the noise (if None, assume no atm)",
    )

    det_data = Unicode(
        defaults.det_data, help="Observation detdata key for the timestream data"
    )

    det_mask = Int(
        defaults.det_mask_nonscience,
        help="Bit mask value for per-detector flagging",
    )

    det_flags = Unicode(
        defaults.det_flags,
        allow_none=True,
        help="Observation detdata key for flags to use",
    )

    det_flag_mask = Int(
        defaults.det_mask_nonscience,
        help="Bit mask value for detector sample flagging",
    )

    shared_flags = Unicode(
        defaults.shared_flags,
        allow_none=True,
        help="Observation shared key for telescope flags to use",
    )

    shared_flag_mask = Int(
        defaults.shared_mask_nonscience,
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

    noise_model = Unicode(
        "noise_model", help="Observation key containing the noise model"
    )

    purge_det_data = Bool(
        False,
        help="If True, clear all observation detector data after copying to mappraiser buffers",
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

    noise_scale = Unicode(
        "noise_scale",
        help="Observation key with optional scaling factor for noise PSDs",
    )

    mem_report = Bool(
        False, help="Print system memory use while staging/unstaging data"
    )

    pair_diff = Bool(
        False,
        help="Process differenced timestreams between orthogonal detectors in pairs.",
    )

    estimate_spin_zero = Bool(
        False,
        help="When doing pair-diff, still estimate a spin-zero field (T) alongside Q and U.",
    )

    # Some traits for noise estimation

    save_psd = Bool(False, help="Save noise PSD information during inv_tt computation")

    apod_window_type = Unicode(
        "chebwin", help="Type of apodisation window to use during noise PSD estimation"
    )

    nperseg = Int(
        0,
        help="If 0, set nperseg = timestream length to compute the noise periodograms.",
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

    # solver
    solver = Int(0, help="Choose mapmaking solver (0->PCG, 1->ECG)")
    tol = Float(1e-12, help="Convergence threshold for the iterative solver")
    maxiter = Int(3000, help="Maximum number of iterations allowed for the solver")
    enlFac = Int(1, help="Enlargement factor for ECG")
    bs_red = Int(0, help="Use dynamic search reduction")

    # preconditioner
    precond = Int(
        0, help="Choose preconditioner (0->BJ, 1->2lvl a priori, 2->2lvl a posteriori)"
    )
    z_2lvl = Int(0, help="Size of 2lvl deflation space")
    ortho_alg = Int(1, help="Orthogonalization scheme for ECG (O->odir, 1->omin)")

    # communication algorithm
    ptcomm_flag = Int(6, help="Choose collective communication scheme")

    # gap treatment strategy
    # 0 -> Condition on gaps having zero signal
    # 1 -> Marginalize on gap contents using 1 extra pixel/scan/detector
    # 2 -> Iterative noise weighting (nested PCG)
    # 3 -> Iterative noise weighting without gaps
    # 4 -> Marginalize on gap contents using 1 extra pixel/proc
    gap_stgy = Int(0, help="Strategy for handling timestream gaps")

    # Gap filling
    # choosing False is only possible on simulations
    # Mappraiser will take the noise from the simulation
    # as if the noise reconstruction inside the gaps was perfect
    do_gap_filling = Bool(True, help="Perform gap filling on the data")

    realization = Int(0, help="Noise realization index (for gap filling)")

    @traitlets.validate("det_mask")
    def _check_det_mask(self, proposal):
        check = proposal["value"]
        if check < 0:
            raise traitlets.TraitError("Det mask should be a positive integer")
        return check

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

    @traitlets.validate("mcmode")
    def _check_mcmode(self, proposal):
        check = proposal["value"]
        if check:
            raise traitlets.TraitError("MC mode is not currently supported")
        return check

    # Checks for mapmaker parameters (solver, ...)
    @traitlets.validate("gap_stgy")
    def _check_gap_stgy(self, proposal):
        check = proposal["value"]
        if check not in (0, 1, 2, 3, 4):
            msg = "Invalid gap_stgy - accepted values are:\n"
            msg += "0 -> condition on gaps having zero signal\n"
            msg += (
                "1 -> marginalize on gap contents using 1 extra pixel/scan/detector\n"
            )
            msg += "2 -> iterative noise weighting 'nested PCG' (completely ignore the gaps)\n"
            msg += "3 -> iterative noise weighting without gaps\n"
            msg += "4 -> marginalize on gap contents using 1 extra pixel/proc"
            raise traitlets.TraitError(msg)
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
            # not implemented
            # libmappraiser.clear_caches()
            self._cached = False

        for atr in [
            "signal",
            "blocksizes",
            "detindxs",
            "obsindxs",
            "telescopes",
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
                                params[k] = v
            if data.comm.comm_world is not None:
                params = data.comm.comm_world.bcast(params, root=0)
            for k, v in self.params.items():
                params[k] = v

        params.update(
            {
                "lambda": self.bandwidth,
                "solver": self.solver,
                "precond": self.precond,
                "Z_2lvl": self.z_2lvl,
                "ptcomm_flag": self.ptcomm_flag,
                "tol": np.double(self.tol),
                "maxiter": self.maxiter,
                "enlFac": self.enlFac,
                "ortho_alg": self.ortho_alg,
                "bs_red": self.bs_red,
                "gap_stgy": self.gap_stgy,
                "do_gap_filling": self.do_gap_filling,
                "realization": self.realization,
            }
        )

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

        # Check if noiseless mode is activated
        if self.noiseless or (self.noise_name is None):
            self.noiseless = True
            self.noise_name = None
            # diagonal noise covariance
            params["lambda"] = 1

        # Pair-differencing checks
        if not self.pair_diff:
            self.estimate_spin_zero = True

        # Log the libmappraiser parameters that were used.
        if data.comm.world_rank == 0:
            with open(
                os.path.join(params["path_output"], "mappraiser_args_log.toml"),
                "w",
            ) as f:
                tomlkit.dump(params, f, sort_keys=True)

        # Check input parameters and compute the sizes of Mappraiser data objects
        if data.comm.world_rank == 0:
            msg = "{} Computing data sizes".format(self._logprefix)
            log.info(msg)
        (
            all_dets,
            nsamp,
            nnz,
            nnz_full,
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

        signal_dtype, data_size_proc, nblock_loc = self._stage_data(
            params,
            data,
            all_dets,
            nsamp,
            nnz,
            nnz_full,
            interval_starts,
            psd_freqs,
        )

        log.info_rank(
            f"{self._logprefix} Staged data in",
            comm=data.comm.comm_world,
            timer=timer,
        )

        # Compute the ML map
        if not kwargs.get("test_gap_fill"):
            if data.comm.world_rank == 0:
                msg = "{} Computing the ML map".format(self._logprefix)
                log.info(msg)
            self._MLmap(
                params,
                data,
                data_size_proc,
                nblock_loc,
                nnz,
            )
        else:
            # merge signal and noise
            self._mappraiser_signal += self._mappraiser_noise

            # save a copy of the original TOD (since it will be modified by gap-filling procedure)
            original_tod = np.copy(self._mappraiser_signal)

            # identify the gaps (1 -> gap, 0 -> valid)
            tgaps = np.zeros_like(self._mappraiser_pixels[::nnz], dtype=np.uint8)
            tgaps[self._mappraiser_pixels[::nnz] < 0] = 1

            # perform gap-filling
            self._gap_filling(
                params,
                data,
                data_size_proc,
                nblock_loc,
                nnz,
            )

            if data.comm.world_rank == 0:
                # save results
                np.savez_compressed(
                    os.path.join(params["path_output"], f"data/gf_{params['ref']}"),
                    tod=original_tod,
                    gaps=tgaps,
                    gf=self._mappraiser_signal,
                )

        log.info_rank(
            f"{self._logprefix} Processed time data in",
            comm=data.comm.comm_world,
            timer=timer,
        )

        return

    def _finalize(self, data, **kwargs):
        self.clear()

    def _requires(self):
        req = {
            "meta": [self.noise_model],
            "shared": [
                self.times,
            ],
            "detdata": [self.det_data],
            "intervals": list(),
        }
        if self.shared_flags is not None:
            req["shared"].append(self.shared_flags)
        if self.det_flags is not None:
            req["detdata"].append(self.det_flags)
        return req

    def _provides(self):
        prov = {"detdata": list()}
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

        all_dets = set()
        nnz_full = None
        psd_freqs = None

        for ob in data.obs:
            # Get the detectors we are using for this observation
            local_dets = ob.select_local_detectors(detectors, flagmask=self.det_mask)
            all_dets.update(set(local_dets))

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

            # Not using views, so there is only one interval per observation
            interval_starts.append(nsamp)
            nsamp += ob.n_local_samples

        interval_starts = np.array(interval_starts, dtype=np.int64)
        all_dets = list(sorted(all_dets))[: self.limit_det]

        nnz_full = len(self.stokes_weights.mode)

        # Check that Stokes weights operator has mode "iqu"
        if nnz_full != 3:
            msg = "Mappraiser assumes that I,Q,U weights are provided\n"
            msg += f"'mode' trait of stokes_weights operator has length {nnz_full} != 3"
            raise RuntimeError(msg)

        if self.pair_diff:
            if self.estimate_spin_zero:
                nnz = nnz_full
            else:
                nnz = nnz_full - 1
        else:
            nnz = nnz_full

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
        #     timestamp_storage, _ = dtype_to_aligned(libmappraiser.TIMESTAMP_TYPE)
        #     self._mappraiser_timestamps_raw = timestamp_storage.zeros(nsamp)
        #     self._mappraiser_timestamps = self_mappraiser_timestamps_raw.array()

        #     interval = 0
        #     time_offset = 0.0

        #     for ob in data.obs:

        #         # Get the noise object for this observation and create new
        #         # entries in the dictionary when the PSD actually changes.  The detector
        #         # weights are obtained from the noise model.

        #         nse = ob[self.noise_model]
        #         nse_scale = 1.0
        #         if self.noise_scale is not None:
        #             if self.noise_scale in ob:
        #                 nse_scale = float(ob[self.noise_scale])

        #         local_dets = set(ob.select_local_detectors(flagmask=self.det_mask))
        #         for det in all_dets:
        #             if det not in local_dets:
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

        # Are we doing pair differencing? If yes, number of dets will be 2 times smaller
        ndet = len(all_dets) if not self.pair_diff else (len(all_dets) // 2)

        # Number of observations
        nobs = len(data.obs)

        # Number of local blocks
        nblock_loc = ndet * nobs

        # Copy the signal.  We always need to do this, even if we are running MCs.

        signal_dtype = data.obs[0].detdata[self.det_data].dtype

        if self._cached:
            # We have previously created the mappraiser buffers.  We just need to fill
            # them from the toast data.  Since both already exist we just copy the
            # contents.
            stage_local(
                data,
                nsamp,
                all_dets,
                self.det_data,
                self._mappraiser_signal,
                interval_starts,
                1,
                self.det_mask,
                None,
                None,
                None,
                self.det_flag_mask,
                do_purge=False,
                pair_diff=self.pair_diff,
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
                    all_dets,
                    self.det_data,
                    libmappraiser.SIGNAL_TYPE,
                    interval_starts,
                    1,
                    self.det_mask,
                    None,
                    None,
                    None,
                    self.det_flag_mask,
                    pair_diff=self.pair_diff,
                )
            else:
                # Allocate and copy all at once.
                storage, _ = dtype_to_aligned(libmappraiser.SIGNAL_TYPE)
                self._mappraiser_signal_raw = storage.zeros(nsamp * ndet)
                self._mappraiser_signal = self._mappraiser_signal_raw.array()

                stage_local(
                    data,
                    nsamp,
                    all_dets,
                    self.det_data,
                    self._mappraiser_signal,
                    interval_starts,
                    1,
                    self.det_mask,
                    None,
                    None,
                    None,
                    self.det_flag_mask,
                    do_purge=False,
                    pair_diff=self.pair_diff,
                )

        # Create buffer for local_block_sizes
        storage, _ = dtype_to_aligned(libmappraiser.PIXEL_TYPE)
        self._mappraiser_blocksizes_raw = storage.zeros(nblock_loc)
        self._mappraiser_blocksizes = self._mappraiser_blocksizes_raw.array()

        # Create buffers for detindx, obindx, telescope
        storage, _ = dtype_to_aligned(np.uint64)
        self._mappraiser_detindxs_raw = storage.zeros(nblock_loc)
        self._mappraiser_obsindxs_raw = storage.zeros(nblock_loc)
        self._mappraiser_telescopes_raw = storage.zeros(nblock_loc)
        self._mappraiser_detindxs = self._mappraiser_detindxs_raw.array()
        self._mappraiser_obsindxs = self._mappraiser_obsindxs_raw.array()
        self._mappraiser_telescopes = self._mappraiser_telescopes_raw.array()

        # Compute sizes of local data blocks, store telescope, ob and det uids (for gap-filling procedure)
        for iobs, ob in enumerate(data.obs):
            dets = ob.select_local_detectors(all_dets, flagmask=self.det_mask)
            sindx = ob.session.uid
            telescope = ob.telescope.uid
            nse = ob[self.noise_model]
            for idet, key in enumerate(nse.all_keys_for_dets(dets)):
                iblock: int
                if self.pair_diff:
                    if idet % 2 == 1:
                        continue
                    else:
                        iblock = (idet // 2) * nobs + iobs
                else:
                    iblock = idet * nobs + iobs
                self._mappraiser_detindxs[iblock] = nse.index(key)
                self._mappraiser_obsindxs[iblock] = sindx
                self._mappraiser_telescopes[iblock] = telescope
                self._mappraiser_blocksizes[iblock] = ob.n_local_samples

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

        # Check if the simulation contains any noise at all.
        if self.noiseless:
            msg = "{} Noiseless mode -> noise buffer filled with zeros".format(
                self._logprefix
            )
            log.info_rank(
                msg,
                data.comm.comm_world,
            )

        # Is there any atmosphere to be added with the noise?
        if self.atm_name is not None:
            if self.atm_name != self.noise_name:
                msg = "{} Adding atmosphere '{}' to noise buffer '{}'".format(
                    self._logprefix, self.atm_name, self.noise_name
                )
                log.info_rank(
                    msg,
                    data.comm.comm_world,
                )
                Combine(
                    op="add",
                    first=self.noise_name,
                    second=self.atm_name,
                    result=self.noise_name,
                ).apply(data)

        if self._cached:
            # We have previously created the mappraiser buffers.  We just need to fill
            # them from the toast data.  Since both already exist we just copy the
            # contents.
            stage_local(
                data,
                nsamp,
                all_dets,
                self.noise_name,
                self._mappraiser_noise,
                interval_starts,
                1,
                self.det_mask,
                None,
                None,
                None,
                self.det_flag_mask,
                do_purge=False,
                pair_diff=self.pair_diff,
            )
        else:
            # Noise buffers do not yet exist
            if self.purge_det_data:
                # Allocate in a staggered way.
                self._mappraiser_noise_raw, self._mappraiser_noise = stage_in_turns(
                    data,
                    nodecomm,
                    n_copy_groups,
                    nsamp,
                    all_dets,
                    self.noise_name,
                    libmappraiser.SIGNAL_TYPE,
                    interval_starts,
                    1,
                    self.det_mask,
                    None,
                    None,
                    None,
                    self.det_flag_mask,
                    pair_diff=self.pair_diff,
                )
            else:
                # Allocate and copy all at once.
                storage, _ = dtype_to_aligned(libmappraiser.SIGNAL_TYPE)
                self._mappraiser_noise_raw = storage.zeros(nsamp * ndet)
                self._mappraiser_noise = self._mappraiser_noise_raw.array()

                stage_local(
                    data,
                    nsamp,
                    all_dets,
                    self.noise_name,
                    self._mappraiser_noise,
                    interval_starts,
                    1,
                    self.det_mask,
                    None,
                    None,
                    None,
                    self.det_flag_mask,
                    do_purge=False,
                    pair_diff=self.pair_diff,
                )
            # Create buffer for invtt and tt
            storage, _ = dtype_to_aligned(libmappraiser.INVTT_TYPE)
            self._mappraiser_invtt_raw = storage.zeros(nblock_loc * params["lambda"])
            self._mappraiser_tt_raw = storage.zeros(nblock_loc * params["lambda"])
            self._mappraiser_invtt = self._mappraiser_invtt_raw.array()
            self._mappraiser_tt = self._mappraiser_tt_raw.array()

        if not self.noiseless:
            # Scale down the noise if we want
            if self.downscale > 1:
                self._mappraiser_noise /= np.sqrt(self.downscale)

            # Compute noise autocorrelation and inverse autocorrelation functions
            compute_autocorrelations(
                nobs,
                ndet,
                self._mappraiser_noise,
                self._mappraiser_blocksizes,
                params["lambda"],
                params["fsample"],
                self._mappraiser_invtt,
                self._mappraiser_tt,
                libmappraiser.INVTT_TYPE,
                apod_window_type=self.apod_window_type,
                print_info=(data.comm.world_rank == 0),
                save_psd=(self.save_psd and data.comm.world_rank == 0),
                save_dir=os.path.join(params["path_output"], "psd_fits"),
            )

            if self.fill_noise_zero:
                # just set the noise to zero
                # that way we can make the mapmaker iterate but on signal-only data
                self._mappraiser_noise[:] = 0.0
        else:
            self._mappraiser_noise[:] = 0.0
            self._mappraiser_invtt[:] = 1.0
            self._mappraiser_tt[:] = 1.0

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
                all_dets,
                self.pixel_pointing.pixels,
                libmappraiser.PIXEL_TYPE,
                interval_starts,
                nnz,
                self.det_mask,
                self.shared_flags,
                self.shared_flag_mask,
                self.det_flags,
                self.det_flag_mask,
                operator=self.pixel_pointing,
                n_repeat=nnz,
                pair_skip=self.pair_diff,
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
                all_dets,
                self.stokes_weights.weights,
                libmappraiser.WEIGHT_TYPE,
                interval_starts,
                nnz,
                self.det_mask,
                None,
                None,
                None,
                self.det_flag_mask,
                operator=self.stokes_weights,
                pair_skip=self.pair_diff,
                select_qu=not self.estimate_spin_zero,  # FIXME is this correct ??
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
        #         psdvals = np.hstack(psdvals).astype(libmappraiser.PSD_TYPE)
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
            nblock_loc,
        )

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

        libmappraiser.MLmap(
            data.comm.comm_world,
            params,
            data_size_proc,
            nb_blocks_loc,
            self._mappraiser_blocksizes,
            self._mappraiser_detindxs,
            self._mappraiser_obsindxs,
            self._mappraiser_telescopes,
            nnz,
            self._mappraiser_pixels,
            self._mappraiser_pixweights,
            self._mappraiser_signal,
            self._mappraiser_noise,
            self._mappraiser_invtt,
            self._mappraiser_tt,
        )

        return

    @function_timer
    def _gap_filling(self, params, data, data_size_proc, nb_blocks_loc, nnz):
        """Perform gap-filling on the data (signal + noise)"""
        libmappraiser.gap_filling(
            data.comm.comm_world,
            data_size_proc,
            nb_blocks_loc,
            self._mappraiser_blocksizes,
            params,
            self._mappraiser_detindxs,
            self._mappraiser_obsindxs,
            self._mappraiser_telescopes,
            nnz,
            self._mappraiser_pixels,
            self._mappraiser_signal,
            self._mappraiser_invtt,
            self._mappraiser_tt,
        )

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
