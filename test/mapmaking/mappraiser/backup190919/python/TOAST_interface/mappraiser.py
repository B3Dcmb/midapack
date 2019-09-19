# This is the interface code of MAPPRAISER with TOAST.
# It contains all the routines to read TOAST data objects and run on the fly
# map-making from TOAST simulations. The code is for the most part copied from the OpMadam model.

#@author: Hamza El Bouhargani
#@date: May 2019


from ctypes.util import find_library
import os

from numpy.fft import fft, fftfreq, fftshift
from scipy import interpolate
import scipy.signal
import math

import ctypes as ct
import healpy as hp
import numpy as np
import numpy.ctypeslib as npc
from toast.cache import Cache
from toast.mpi import MPI
from toast.op import Operator
import toast.timing as timing

try:
    import mappraiser_wrapper as mappraiser
except:
    mappraiser = None

import gc

try:
    import psutil

    def memreport(comm=None, msg=""):
        """ Gather and report the amount of allocated, free and swapped system memory
        """
        if psutil is None:
            return
        vmem = psutil.virtual_memory()._asdict()
        gc.collect()
        vmem2 = psutil.virtual_memory()._asdict()
        memstr = "Memory usage {}\n".format(msg)
        for key, value in vmem.items():
            value2 = vmem2[key]
            if comm is None:
                vlist = [value]
                vlist2 = [value2]
            else:
                vlist = comm.gather(value)
                vlist2 = comm.gather(value2)
            if comm is None or comm.rank == 0:
                vlist = np.array(vlist, dtype=np.float64)
                vlist2 = np.array(vlist2, dtype=np.float64)
                if key != "percent":
                    # From bytes to better units
                    if np.amax(vlist) < 2 ** 20:
                        vlist /= 2 ** 10
                        vlist2 /= 2 ** 10
                        unit = "kB"
                    elif np.amax(vlist) < 2 ** 30:
                        vlist /= 2 ** 20
                        vlist2 /= 2 ** 20
                        unit = "MB"
                    else:
                        vlist /= 2 ** 30
                        vlist2 /= 2 ** 30
                        unit = "GB"
                else:
                    unit = "% "
                if comm is None or comm.size == 1:
                    memstr += "{:>12} : {:8.3f} {}\n".format(key, vlist[0], unit)
                    if np.abs(vlist2[0] - vlist[0]) / vlist[0] > 1e-3:
                        memstr += "{:>12} : {:8.3f} {} (after GC)\n".format(
                            key, vlist2[0], unit
                        )
                else:
                    med1 = np.median(vlist)
                    memstr += (
                        "{:>12} : {:8.3f} {}  < {:8.3f} +- {:8.3f} {}  "
                        "< {:8.3f} {}\n".format(
                            key,
                            np.amin(vlist),
                            unit,
                            med1,
                            np.std(vlist),
                            unit,
                            np.amax(vlist),
                            unit,
                        )
                    )
                    med2 = np.median(vlist2)
                    if np.abs(med2 - med1) / med1 > 1e-3:
                        memstr += (
                            "{:>12} : {:8.3f} {}  < {:8.3f} +- {:8.3f} {}  "
                            "< {:8.3f} {} (after GC)\n".format(
                                key,
                                np.amin(vlist2),
                                unit,
                                med2,
                                np.std(vlist2),
                                unit,
                                np.amax(vlist2),
                                unit,
                            )
                        )
        if comm is None or comm.rank == 0:
            print(memstr, flush=True)
        if comm is not None:
            comm.Barrier()
        return


except:

    def memreport(comm=None, msg=""):
        return

class OpMappraiser(Operator):
    """
    Operator which passes data to libmappraiser for map-making.
    Args:
        params (dictionary): parameters to mappraiser
        detweights (dictionary): individual noise weights to use for each
            detector.
        pixels (str): the name of the cache object (<pixels>_<detector>)
            containing the pixel indices to use.
        pixels_nested (bool): Set to False if the pixel numbers are in
            ring ordering. Default is True.
        weights (str): the name of the cache object (<weights>_<detector>)
            containing the pointing weights to use.
        name (str): the name of the cache object (<name>_<detector>) to
            use for the detector timestream.  If None, use the TOD.
        noise_name (str) : the name of the cache object (<name>_<detector>) to
            use for the noise timestream. If None, skip.
        flag_name (str): the name of the cache object
            (<flag_name>_<detector>) to use for the detector flags.
            If None, use the TOD.
        flag_mask (int): the integer bit mask (0-255) that should be
            used with the detector flags in a bitwise AND.
        common_flag_name (str): the name of the cache object
            to use for the common flags.  If None, use the TOD.
        common_flag_mask (int): the integer bit mask (0-255) that should
            be used with the common flags in a bitwise AND.
        apply_flags (bool): whether to apply flags to the pixel numbers.
        purge (bool): if True, clear any cached data that is copied into
            the Mappraiser buffers.
        purge_tod (bool): if True, clear any cached signal that is
            copied into the Mappraiser buffers.
        purge_pixels (bool): if True, clear any cached pixels that are
            copied into the Mappraiser buffers.
        purge_weights (bool): if True, clear any cached weights that are
            copied into the Mappraiser buffers.
        purge_flags (bool): if True, clear any cached flags that are
            copied into the Mappraiser buffers.
        dets (iterable):  List of detectors to map. If left as None, all
            available detectors are mapped.
        noise (str): Keyword to use when retrieving the noise object
            from the observation.
        conserve_memory(bool/int): Stagger the Mappraiser buffer staging on node.
        translate_timestamps(bool): Translate timestamps to enforce
            monotonity.
    """

    def __init__(
        self,
        params={},
        detweights=None,
        pixels="pixels",
        pixels_nested=True,
        weights="weights",
        name="signal",
        noise_name=None,
        flag_name=None,
        flag_mask=255,
        common_flag_name=None,
        common_flag_mask=255,
        apply_flags=False,
        purge=False,
        dets=None,
        purge_tod=False,
        purge_pixels=False,
        purge_weights=False,
        purge_flags=False,
        noise="noise",
        intervals="intervals",
        conserve_memory=False,
        translate_timestamps=True,
    ):
        # We call the parent class constructor, which currently does nothing
        super().__init__()
        # MAPPRAISER uses time-based distribution
        self._name = name
        self._noise_name = noise_name
        self._flag_name = flag_name
        self._flag_mask = flag_mask
        self._common_flag_name = common_flag_name
        self._common_flag_mask = common_flag_mask
        self._pixels = pixels
        self._pixels_nested = pixels_nested
        self._weights = weights
        self._detw = detweights
        self._purge = purge
        if self._purge:
            self._purge_tod = True
            self._purge_pixels = True
            self._purge_weights = True
            self._purge_flags = True
        else:
            self._purge_tod = purge_tod
            self._purge_pixels = purge_pixels
            self._purge_weights = purge_weights
            self._purge_flags = purge_flags
        self._apply_flags = apply_flags
        self._params = params
        if dets is not None:
            self._dets = set(dets)
        else:
            self._dets = None
        self._noisekey = noise
        self._intervals = intervals
        self._cache = Cache()
        self._mappraiser_timestamps = None
        self._mappraiser_noise = None
        self._mappraiser_pixels = None
        self._mappraiser_pixweights = None
        self._mappraiser_signal = None
        self._mappraiser_invtt = None
        self._conserve_memory = int(conserve_memory)
        self._translate_timestamps = translate_timestamps
        self._verbose = True

    def __del__(self):
        self._cache.clear()

    @property
    def available(self):
        """
        (bool): True if libmappraiser is found in the library search path.
        """
        return mappraiser is not None and mappraiser.available

    def exec(self, data, comm=None):
        """
        Copy data to Mappraiser-compatible buffers and make a map.
        Args:
            data (toast.Data): The distributed data.
        """
        if not self.available:
            raise RuntimeError("libmappraiser is not available")

        if len(data.obs) == 0:
            raise RuntimeError(
                "OpMappraiser requires every supplied data object to "
                "contain at least one observation"
            )

        auto_timer = timing.auto_timer(type(self).__name__)

        if comm is None:
            # Just use COMM_WORLD
            comm = data.comm.comm_world

        (
            dets,
            nsamp,
            ndet,
            nnz,
            nnz_full,
            nnz_stride,
            psdfreqs,
            nside,
        ) = self._prepare(data, comm)

        Lambda = self._params["Lambda"]

        data_size_proc, nobsloc, local_blocks_sizes, signal_type, noise_type, pixels_dtype, weight_dtype = self._stage_data(
            data,
            comm,
            nsamp,
            ndet,
            nnz,
            nnz_full,
            nnz_stride,
            psdfreqs,
            Lambda,
            dets,
            nside,
        )

        self._MLmap(comm, data_size_proc, nobsloc*ndet, local_blocks_sizes, nnz, Lambda)

        self._unstage_data(
            comm,
            data,
            nsamp,
            nnz,
            nnz_full,
            dets,
            signal_type,
            noise_type,
            pixels_dtype,
            nside,
            weight_dtype,
        )

        return

    def _MLmap(self, comm, data_size_proc, nb_blocks_loc, local_blocks_sizes, nnz, Lambda):
        """ Compute the ML map
        """
        auto_timer = timing.auto_timer(type(self).__name__)
        if self._verbose:
            memreport(comm, "just before calling libmappraiser.MLmap")

        # Compute the Maximum Likelihood map
        os.environ["OMP_NUM_THREADS"] = "1"
        mappraiser.MLmap(
            comm,
            self._params,
            data_size_proc,
            nb_blocks_loc,
            local_blocks_sizes,
            nnz,
            self._mappraiser_pixels,
            self._mappraiser_pixweights,
            self._mappraiser_signal,
            self._mappraiser_noise,
            Lambda,
            self._mappraiser_invtt,
            )
        os.environ["OMP_NUM_THREADS"] = "4"

        return

    def _count_samples(self, data):
        """ Loop over the observations and count the number of samples.
        """
        if len(data.obs) != 1:
            nsamp = 0
            tod0 = data.obs[0]["tod"]
            detectors0 = tod0.local_dets
            for obs in data.obs:
                tod = obs["tod"]
                # For the moment, we require that all observations have
                # the same set of detectors
                detectors = tod.local_dets
                dets_are_same = True
                if len(detectors0) != len(detectors):
                    dets_are_same = False
                else:
                    for det1, det2 in zip(detectors0, detectors):
                        if det1 != det2:
                            dets_are_same = False
                            break
                if not dets_are_same:
                    raise RuntimeError(
                        "When calling Mappraiser, all TOD assigned to a process "
                        "must have the same local detectors."
                    )
                nsamp += tod.local_samples[1]
        else:
            tod = data.obs[0]["tod"]
            nsamp = tod.local_samples[1]
        return nsamp

    def _get_period_ranges(self, data, detectors):
        """ Collect the ranges of every observation.
        """
        # # Discard intervals that are too short to fit a baseline
        # if "basis_order" in self.params:
        #     norder = int(self.params["basis_order"]) + 1
        # else:
        #     norder = 1
        # norder = 1

        psdfreqs = None
        # period_lengths = []
        # obs_period_ranges = []

        for obs in data.obs:
            tod = obs["tod"]
            # Check that all noise objects have the same binning
            if self._noisekey in obs.keys():
                nse = obs[self._noisekey]
                if nse is not None:
                    if psdfreqs is None:
                        psdfreqs = nse.freq(detectors[0]).astype(np.float64).copy()
                    for det in detectors:
                        check_psdfreqs = nse.freq(det)
                        if not np.allclose(psdfreqs, check_psdfreqs):
                            raise RuntimeError(
                                "All PSDs passed to Mappraiser must have"
                                " the same frequency binning."
                            )
        #     # Collect the valid intervals for this observation
        #     period_ranges = []
        #     if self._intervals in obs:
        #         intervals = obs[self._intervals]
        #     else:
        #         intervals = None
        #     local_intervals = tod.local_intervals(intervals)
        #
        #     for ival in local_intervals:
        #         local_start = ival.first
        #         local_stop = ival.last + 1
        #         if local_stop - local_start < norder:
        #             continue
        #         period_lengths.append(local_stop - local_start)
        #         period_ranges.append((local_start, local_stop))
        #     obs_period_ranges.append(period_ranges)
        #
        # # Update the number of samples based on the valid intervals
        #
        # nsamp_tot_full = comm.allreduce(nsamp, op=MPI.SUM)
        # nperiod = len(period_lengths)
        # period_lengths = np.array(period_lengths, dtype=np.int64)
        # nsamp = np.sum(period_lengths, dtype=np.int64)
        # nsamp_tot = comm.allreduce(nsamp, op=MPI.SUM)
        # if nsamp_tot == 0:
        #     raise RuntimeError(
        #         "No samples in valid intervals: nsamp_tot_full = {}, "
        #         "nsamp_tot = {}".format(nsamp_tot_full, nsamp_tot)
        #     )
        # if comm.rank == 0:
        #     print(
        #         "OpMappraiser: {:.2f} % of samples are included in valid "
        #         "intervals.".format(nsamp_tot * 100.0 / nsamp_tot_full), flush = True
        #     )

        # # Mappraiser expects starting indices, not period lengths
        # periods = np.zeros(nperiod, dtype=np.int64)
        # for i, n in enumerate(period_lengths[:-1]):
        #     periods[i + 1] = periods[i] + n

        return psdfreqs

    def _psd2invtt(self, psdfreqs, psd, Lambda):
        """ Generate the first row of the Toeplitz blocks from the psds
        """
        # parameters
        sampling_freq = self._params["samplerate"]
        f_defl = sampling_freq/(np.pi*Lambda)
        df = f_defl/2
        block_size = 2**(int(math.log(sampling_freq*1./df,2))+1)

        # Extracting psd info
        psd_sim_m1 = np.reciprocal(psd)

        # Initialize full size inverse PSD in frequency domain
        fs = fftfreq(block_size, 1./sampling_freq)
        psdm1 = np.zeros_like(fs)

        # Perform interpolation to get full size PSD from TOAST provided PSD
        tck = interpolate.splrep(psdfreqs, psd_sim_m1, s=0) #s=0 : no smoothing
        psdfit = interpolate.splev(np.abs(fs[:int(block_size/2)+1]), tck, der=0)
        psdfit[0] = 0 #set offset noise contribution to zero
        psdm1[:int(block_size/2)] = psdfit[:int(block_size/2)]
        psdm1[int(block_size/2):] = np.flip(psdfit[1:],0)

        # Compute noise autocorrelation and inverse noise autocorrelation functions
        inv_tt = np.real(np.fft.ifft(psdm1, n=block_size))

        # Define apodization window
        window = scipy.signal.gaussian(2*Lambda, 1./2*Lambda)
        window = np.fft.ifftshift(window)
        window = window[:Lambda]
        window = np.pad(window,(0,int(block_size/2-(Lambda))),'constant')
        symw = np.zeros(block_size)
        symw[:int(block_size/2)] = window
        symw[int(block_size/2):] = np.flip(window,0)

        inv_tt_w = np.multiply(symw, inv_tt, dtype = mappraiser.INVTT_TYPE)

        return inv_tt_w[:Lambda]

    def _prepare(self, data, comm):
        """ Examine the data object.
        """
        auto_timer = timing.auto_timer(type(self).__name__)

        nsamp = self._count_samples(data)

        # Determine the detectors and the pointing matrix non-zeros
        # from the first observation. Mappraiser will expect these to remain
        # unchanged across observations.

        tod = data.obs[0]["tod"]

        if self._dets is None:
            dets = tod.local_dets
        else:
            dets = [det for det in tod.local_dets if det in self._dets]
        ndet = len(dets)

        # We get the number of Non-zero pointing weights per pixel, from the
        # shape of the data from the first detector

        nnzname = "{}_{}".format(self._weights, dets[0])
        nnz_full = tod.cache.reference(nnzname).shape[1]

        if nnz_full != 3:
            raise RuntimeError(
                    "OpMappraiser: Don't know how to make a map "
                    "with nnz={}".format(nnz_full)
                )
            nnz = 3
            nnz_stride = 1
        else:
            nnz = nnz_full
            nnz_stride = 1

        if "nside" not in self._params:
            raise RuntimeError(
                'OpMappraiser: "nside" must be set in the parameter dictionary'
            )
        nside = int(self._params["nside"])


        # Inspect the valid intervals across all observations to
        # determine the number of samples per detector
        psdfreqs = self._get_period_ranges(data, dets)

        return (
            dets,
            nsamp,
            ndet,
            nnz,
            nnz_full,
            nnz_stride,
            psdfreqs,
            nside,
        )

    def _stage_time(self, data, detectors, nsamp, psdfreqs, Lambda):
        """ Stage the timestamps and use them to build PSD inputs.
        """
        auto_timer = timing.auto_timer(type(self).__name__)
        # self._mappraiser_timestamps = self._cache.create(
        #     "timestamps", mappraiser.TIMESTAMP_TYPE, (nsamp,)
        # )

        # offset = 0
        # time_offset = 0
        # psds = {}
        invtt_list = []
        for iobs, obs in enumerate(data.obs):
            tod = obs["tod"]
            # period_ranges = obs_period_ranges[iobs]

            # # Collect the timestamps for the valid intervals
            # timestamps = tod.local_times().copy()
            # if self._translate_timestamps:
            #     # Translate the time stamps to be monotonous
            #     timestamps -= timestamps[0] - time_offset
            #     time_offset = timestamps[-1] + 1
            #
            # for istart, istop in period_ranges:
            #     nn = istop - istart
            #     ind = slice(offset, offset + nn)
            #     self._mappraiser_timestamps[ind] = timestamps[istart:istop]
            #     offset += nn

            # get the noise object for this observation and create new
            # entries in the dictionary when the PSD actually changes
            if self._noisekey in obs.keys():
                nse = obs[self._noisekey]
                if "noise_scale" in obs:
                    noise_scale = obs["noise_scale"]
                else:
                    noise_scale = 1
                if nse is not None:
                    for det in detectors:
                        psd = nse.psd(det) * noise_scale ** 2
                        invtt = self._psd2invtt(psdfreqs, psd, Lambda)
                        invtt_list.append(invtt)
                        # if det not in psds:
                        #     psds[det] = [(0, psd)]
                        # else:
                        #     if not np.allclose(psds[det][-1][1], psd):
                        #         psds[det] += [(timestamps[0], psd)]

        return invtt_list

    def _stage_signal(self, data, detectors, nsamp, ndet):
        """ Stage signal
        """
        auto_timer = timing.auto_timer(type(self).__name__)
        self._mappraiser_signal = self._cache.create(
            "signal", mappraiser.SIGNAL_TYPE, (nsamp * ndet,)
        )
        self._mappraiser_signal[:] = np.nan

        global_offset = 0
        local_blocks_sizes = []
        for iobs, obs in enumerate(data.obs):
            tod = obs["tod"]
            # period_ranges = obs_period_ranges[iobs]

            for idet, det in enumerate(detectors):
                # Get the signal.
                signal = tod.local_signal(det, self._name)
                signal_dtype = signal.dtype
                offset = global_offset
                local_V_size = len(signal)
                dslice = slice(idet * nsamp + offset, idet * nsamp + offset + local_V_size)
                self._mappraiser_signal[dslice] = signal
                offset += local_V_size
                local_blocks_sizes.append(local_V_size)


                del signal

            for idet, det in enumerate(detectors):
                if self._name is not None and (
                    self._purge_tod
                ):
                    cachename = "{}_{}".format(self._name, det)
                    tod.cache.clear(pattern=cachename)

            global_offset = offset

        local_blocks_sizes = np.array(local_blocks_sizes, dtype=np.int32)

        return signal_dtype, local_blocks_sizes

    def _stage_noise(self, data, detectors, nsamp, ndet):
        """ Stage noise timestream (detector noise + atmosphere)
        """
        auto_timer = timing.auto_timer(type(self).__name__)
        self._mappraiser_noise = self._cache.create(
            "noise", mappraiser.SIGNAL_TYPE, (nsamp * ndet,)
        )
        if self._noise_name == None:
            self._mappraiser_noise = np.zeros_like(self._mappraiser_noise)
            return self._mappraiser_noise.dtype

        self._mappraiser_noise[:] = np.nan

        global_offset = 0
        for iobs, obs in enumerate(data.obs):
            tod = obs["tod"]
            # period_ranges = obs_period_ranges[iobs]

            for idet, det in enumerate(detectors):
                # Get the signal.
                noise = tod.local_signal(det, self._noise_name)
                noise_dtype = noise.dtype
                offset = global_offset
                nn = len(noise)
                dslice = slice(idet * nsamp + offset, idet * nsamp + offset + nn)
                self._mappraiser_noise[dslice] = noise
                offset += nn


                del noise

            for idet, det in enumerate(detectors):
                if self._noise_name is not None and (
                    self._purge_tod
                ):
                    cachename = "{}_{}".format(self._noise_name, det)
                    tod.cache.clear(pattern=cachename)

            global_offset = offset


        return noise_dtype

    def _stage_pixels(self, data, detectors, nsamp, ndet, nnz, nside):
        """ Stage pixels
        """
        auto_timer = timing.auto_timer(type(self).__name__)
        self._mappraiser_pixels = self._cache.create(
            "pixels", mappraiser.PIXEL_TYPE, (nsamp * ndet * nnz,)
        )
        self._mappraiser_pixels[:] = -1

        global_offset = 0
        for iobs, obs in enumerate(data.obs):
            tod = obs["tod"]
            # period_ranges = obs_period_ranges[iobs]

            commonflags = None
            for idet, det in enumerate(detectors):
                # Optionally get the flags, otherwise they are
                # assumed to have been applied to the pixel numbers.
                # Mappraiser doesn't use flags for now but might be useful for
                # future updates.

                if self._apply_flags:
                    detflags = tod.local_flags(det, self._flag_name)
                    commonflags = tod.local_common_flags(self._common_flag_name)
                    flags = np.logical_or(
                        (detflags & self._flag_mask) != 0,
                        (commonflags & self._common_flag_mask) != 0,
                    )
                    del detflags

                # get the pixels for the valid intervals from the cache

                pixelsname = "{}_{}".format(self._pixels, det)
                pixels = tod.cache.reference(pixelsname)
                pixels_dtype = pixels.dtype

                if not self._pixels_nested:
                    # Madam expects the pixels to be in nested ordering.
                    # This is not the case for Mappraiser but keeping it for now
                    pixels = pixels.copy()
                    good = pixels >= 0
                    pixels[good] = hp.ring2nest(nside, pixels[good])

                if self._apply_flags:
                    pixels = pixels.copy()
                    pixels[flags] = -1

                offset = global_offset
                nn = len(pixels)
                dslice = slice(
                    (idet * nsamp + offset) * nnz,
                    (idet * nsamp + offset + nn) * nnz,
                )
                # nnz = 3 is a mandatory assumption here (could easily be generalized ...)
                self._mappraiser_pixels[dslice] = nnz * np.repeat(pixels,nnz)
                self._mappraiser_pixels[dslice][1::nnz] += 1
                self._mappraiser_pixels[dslice][2::nnz] += 2
                offset += nn

                del pixels
                if self._apply_flags:
                    del flags

            # Always purge the pixels but restore them from the Madam
            # buffers when purge_pixels=False
            for idet, det in enumerate(detectors):
                pixelsname = "{}_{}".format(self._pixels, det)
                tod.cache.clear(pattern=pixelsname)
                if self._name is not None and (
                    self._purge_tod
                ):
                    cachename = "{}_{}".format(self._name, det)
                    tod.cache.clear(pattern=cachename)
                if self._purge_flags and self._flag_name is not None:
                    cacheflagname = "{}_{}".format(self._flag_name, det)
                    tod.cache.clear(pattern=cacheflagname)

            del commonflags
            if self._purge_flags and self._common_flag_name is not None:
                tod.cache.clear(pattern=self._common_flag_name)
            global_offset = offset

        return pixels_dtype

    def _stage_pixweights(
        self, data, detectors, nsamp, ndet, nnz, nnz_full, nnz_stride
    ):
        """Now collect the pixel weights
        """
        auto_timer = timing.auto_timer(type(self).__name__)

        self._mappraiser_pixweights = self._cache.create(
            "pixweights", mappraiser.WEIGHT_TYPE, (nsamp * ndet * nnz,)
        )
        self._mappraiser_pixweights[:] = 0

        global_offset = 0
        for iobs, obs in enumerate(data.obs):
            tod = obs["tod"]
            # period_ranges = obs_period_ranges[iobs]
            for idet, det in enumerate(detectors):
                # get the pixels and weights for the valid intervals
                # from the cache
                weightsname = "{}_{}".format(self._weights, det)
                weights = tod.cache.reference(weightsname)
                weight_dtype = weights.dtype
                offset = global_offset
                nn = len(weights)
                dwslice = slice(
                    (idet * nsamp + offset) * nnz,
                    (idet * nsamp + offset + nn) * nnz,
                )
                self._mappraiser_pixweights[dwslice] = weights.flatten()[
                    ::nnz_stride
                ]
                offset += nn
                del weights
            # Purge the weights but restore them from the Madam
            # buffers when purge_weights=False.
            # Handle special case when Madam only stores a subset of
            # the weights.
            if not self._purge_weights and (nnz != nnz_full):
                pass
            else:
                for idet, det in enumerate(detectors):
                    # get the pixels and weights for the valid intervals
                    # from the cache
                    weightsname = "{}_{}".format(self._weights, det)
                    tod.cache.clear(pattern=weightsname)

            global_offset = offset

        return weight_dtype

    def _stage_data(
        self,
        data,
        comm,
        nsamp,
        ndet,
        nnz,
        nnz_full,
        nnz_stride,
        psdfreqs,
        Lambda,
        detectors,
        nside,
    ):
        """ create Mappraiser-compatible buffers
        Collect the TOD into Mappraiser buffers. Process pixel weights
        Separate from the rest to reduce the memory high water mark
        When the user has set purge=True
        Moving data between toast and Mappraiser buffers has an overhead.
        We perform the operation in a staggered fashion to have the
        overhead only once per node.
        """
        auto_timer = timing.auto_timer(type(self).__name__)

        if self._conserve_memory:
            # The user has elected to stagger staging the data on each
            # node to avoid exhausting memory
            nodecomm = comm.Split_type(MPI.COMM_TYPE_SHARED, comm.rank)
            if self._conserve_memory == 1:
                nread = nodecomm.size
            else:
                nread = min(self._conserve_memory, nodecomm.size)
        else:
            nodecomm = MPI.COMM_SELF
            nread = 1

        for iread in range(nread):
            nodecomm.Barrier()
            if nodecomm.rank % nread != iread:
                continue
            invtt_list = self._stage_time(data, detectors, nsamp, psdfreqs, Lambda)
            self._mappraiser_invtt = np.array([np.array(invtt_i, dtype= mappraiser.INVTT_TYPE) for invtt_i in invtt_list])
            del invtt_list
            self._mappraiser_invtt = np.concatenate(self._mappraiser_invtt)
            signal_dtype, local_blocks_sizes = self._stage_signal(
                data, detectors, nsamp, ndet
            )
            noise_dtype = self._stage_noise(
                data, detectors, nsamp, ndet
            )
            pixels_dtype = self._stage_pixels(
                data, detectors, nsamp, ndet, nnz, nside
            )
            weight_dtype = self._stage_pixweights(
                data,
                detectors,
                nsamp,
                ndet,
                nnz,
                nnz_full,
                nnz_stride,
            )
        del nodecomm

        # detweights is either a dictionary of weights specified at
        # construction time, or else we use uniform weighting.
        detw = {}
        if self._detw is None:
            for idet, det in enumerate(detectors):
                detw[det] = 1.0
        else:
            detw = self._detw

        detweights = np.zeros(ndet, dtype=np.float64)
        for idet, det in enumerate(detectors):
            detweights[idet] = detw[det]

        # if len(psds) > 0:
        #     npsdbin = len(psdfreqs)
        #
        #     npsd = np.zeros(ndet, dtype=np.int64)
        #     psdstarts = []
        #     psdvals = []
        #     for idet, det in enumerate(detectors):
        #         if det not in psds:
        #             raise RuntimeError("Every detector must have at least " "one PSD")
        #         psdlist = psds[det]
        #         npsd[idet] = len(psdlist)
        #         for psdstart, psd in psdlist:
        #             psdstarts.append(psdstart)
        #             psdvals.append(psd)
        #     npsdtot = np.sum(npsd)
        #     psdstarts = np.array(psdstarts, dtype=np.float64)
        #     psdvals = np.hstack(psdvals).astype(mappraiser.PSD_TYPE)
        #     npsdval = psdvals.size
        # else:
        #     print("white \n")
        #     npsd = np.ones(ndet, dtype=np.int64)
        #     npsdtot = np.sum(npsd)
        #     psdstarts = np.zeros(npsdtot)
        #     npsdbin = 10
        #     fsample = 10.0
        #     psdfreqs = np.arange(npsdbin) * fsample / npsdbin
        #     npsdval = npsdbin * npsdtot
        #     psdvals = np.ones(npsdval)
        # psdinfo = (detweights, npsd, psdstarts, psdfreqs, psdvals)

        # Get global array of data sizes of the full communicator
        data_size_proc = np.array(comm.allgather(len(self._mappraiser_signal)), dtype=np.int32)
        # Get number of local observations
        nobsloc = len(data.obs)

        return data_size_proc, nobsloc, local_blocks_sizes, signal_dtype, noise_dtype, pixels_dtype, weight_dtype

    def _unstage_data(
        self,
        comm,
        data,
        nsamp,
        nnz,
        nnz_full,
        detectors,
        signal_type,
        noise_type,
        pixels_dtype,
        nside,
        weight_dtype,
    ):
        """ Clear Mappraiser buffers, restore pointing into TOAST caches.
        """
        auto_timer = timing.auto_timer(type(self).__name__)
        # self._mappraiser_timestamps = None
        # self._cache.destroy("timestamps")

        if self._conserve_memory:
            nodecomm = comm.Split_type(MPI.COMM_TYPE_SHARED, comm.rank)
            nread = nodecomm.size
        else:
            nodecomm = MPI.COMM_SELF
            nread = 1

        for iread in range(nread):
            nodecomm.Barrier()
            if nodecomm.rank % nread != iread:
                continue
            self._mappraiser_signal = None
            self._cache.destroy("signal")
            self._mappraiser_noise = None
            self._cache.destroy("noise")

            # if not self._purge_pixels:
            #     # restore the pixels from the Mappraiser buffers
            #     global_offset = 0
            #     for obs, period_ranges in zip(data.obs, obs_period_ranges):
            #         tod = obs["tod"]
            #         nlocal = tod.local_samples[1]
            #         for idet, det in enumerate(detectors):
            #             pixels = -np.ones(nlocal, dtype=pixels_dtype)
            #             offset = global_offset
            #             for istart, istop in period_ranges:
            #                 nn = istop - istart
            #                 dslice = slice(
            #                     idet * nsamp + offset, idet * nsamp + offset + nn
            #                 )
            #                 pixels[istart:istop] = self._mappraiser_pixels[dslice]
            #                 offset += nn
            #             npix = 12 * nside ** 2
            #             good = np.logical_and(pixels >= 0, pixels < npix)
            #             if not self._pixels_nested:
            #                 pixels[good] = hp.nest2ring(nside, pixels[good])
            #             pixels[np.logical_not(good)] = -1
            #             cachename = "{}_{}".format(self._pixels, det)
            #             tod.cache.put(cachename, pixels, replace=True)
            #         global_offset = offset
            self._mappraiser_pixels = None
            self._cache.destroy("pixels")

            # if not self._purge_weights and nnz == nnz_full:
            #     # restore the weights from the Madam buffers
            #     global_offset = 0
            #     for obs, period_ranges in zip(data.obs, obs_period_ranges):
            #         tod = obs["tod"]
            #         nlocal = tod.local_samples[1]
            #         for idet, det in enumerate(detectors):
            #             weights = np.zeros([nlocal, nnz], dtype=weight_dtype)
            #             offset = global_offset
            #             for istart, istop in period_ranges:
            #                 nn = istop - istart
            #                 dwslice = slice(
            #                     (idet * nsamp + offset) * nnz,
            #                     (idet * nsamp + offset + nn) * nnz,
            #                 )
            #                 weights[istart:istop] = self._mappraiser_pixweights[
            #                     dwslice
            #                 ].reshape([-1, nnz])
            #                 offset += nn
            #             cachename = "{}_{}".format(self._weights, det)
            #             tod.cache.put(cachename, weights, replace=True)
            #         global_offset = offset
            self._mappraiser_pixweights = None
            self._cache.destroy("pixweights")
        del nodecomm
        return
