# This is the interface code of MAPPRAISER with TOAST.
# It contains all the routines to read TOAST data objects and run on the fly
# map-making from TOAST simulations. The code is for the most part copied from the OpMadam model.

#@author: Hamza El Bouhargani
#@date: May 2019
#@latest_update: January 2020

from toast.mpi import MPI, use_mpi

import os

import healpy as hp
import numpy as np

from toast.cache import Cache
from toast.op import Operator
from toast.timing import function_timer, Timer
from toast.utils import Logger, memreport
import toast.qarray as qa

from numpy.fft import fft, fftfreq, fftshift
from scipy import interpolate
import scipy.signal
import math


mappraiser = None
if use_mpi:
    try:
        import mappraiser_wrapper as mappraiser
    except ImportError:
        mappraiser = None

def count_caches(data, comm, nodecomm, mappraisercache, msg=""):
    """ Count the amount of memory in the TOD caches
    """
    my_todsize = 0
    for obs in data.obs:
        tod = obs["tod"]
        my_todsize += tod.cache.report(silent=True)
    my_cachesize = mappraisercache.report(silent=True)
    node_todsize = nodecomm.allreduce(my_todsize, MPI.SUM)
    node_cachesize = nodecomm.allreduce(my_cachesize, MPI.SUM)
    if comm.rank == 0:
        print(
            "Node has {:.3f} GB allocated in TOAST TOD caches and "
            "{:.3f} GB in Mappraiser caches ({:.3f} GB total) {}".format(
                node_todsize / 2 ** 30,
                node_cachesize / 2 ** 30,
                (node_todsize + node_cachesize) / 2 ** 30,
                msg,
            )
        )
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
        # Call the parent class constructor
        super().__init__()

        # mappraiser uses time-based distribution
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

    @function_timer
    def exec(self, data, comm=None):
        """
        Copy data to Mappraiser-compatible buffers and make a map.

        Args:
            data (toast.Data): The distributed data.

        Returns:
            None
        """
        if not self.available:
            raise RuntimeError("libmappraiser is not available")

        if len(data.obs) == 0:
            raise RuntimeError(
                "OpMappraiser requires every supplied data object to "
                "contain at least one observation"
            )


        if comm is None:
            # Use the word communicator from the distributed data.
            comm = data.comm.comm_world
        self._data = data
        self._comm = comm
        self._rank = comm.rank

        (
            dets,
            nsamp,
            ndet,
            nnz,
            nnz_full,
            nnz_stride,
            psdfreqs,
            nside,
        ) = self._prepare()

        data_size_proc, nobsloc, local_blocks_sizes, signal_type, noise_type, pixels_dtype, sweeptstamps, nsweeps, az, az_min, az_max, nces, weight_dtype = self._stage_data(
            nsamp,
            ndet,
            nnz,
            nnz_full,
            nnz_stride,
            psdfreqs,
            dets,
            nside,
        )
        if self._params["map-maker"] == 'ML':
            self._MLmap(data_size_proc, nobsloc*ndet, local_blocks_sizes, nnz)
        elif self._params["map-maker"] == 'MT':
            self._MTmap(sweeptstamps, nsweeps, az, az_min, az_max, nces, data_size_proc, nobsloc*ndet, local_blocks_sizes, nnz)
        else:
            raise RuntimeError(
                "Unvalid Map-making technique please choose:"
                " -'ML' for Maximum Likelihood map-Making"
                " -'MT' for Marginalized Templates map-making"
            )

        self._unstage_data(
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

    @function_timer
    def _MLmap(self, data_size_proc, nb_blocks_loc, local_blocks_sizes, nnz):
        """ Compute the ML map
        """
        if self._verbose:
            memreport("just before calling libmappraiser.MLmap", self._comm)

        # Compute the Maximum Likelihood map
        os.environ["OMP_NUM_THREADS"] = "1"
        mappraiser.MLmap(
            self._comm,
            self._params,
            data_size_proc,
            nb_blocks_loc,
            local_blocks_sizes,
            nnz,
            self._mappraiser_pixels,
            self._mappraiser_pixweights,
            self._mappraiser_signal,
            self._mappraiser_noise,
            self._params["Lambda"],
            self._mappraiser_invtt,
            )
        # os.environ["OMP_NUM_THREADS"] = "4"

        return

    @function_timer
    def _MTmap(self, sweeptstamps, nsweeps, az, az_min, az_max, nces, data_size_proc, nb_blocks_loc, local_blocks_sizes, nnz):
        """ Compute the Marginalized templates map
        """
        if self._verbose:
            memreport("just before calling libmappraiser.MTmap", self._comm)

        # Compute the marginalized templates map
        os.environ["OMP_NUM_THREADS"] = "1"
        mappraiser.MTmap(
            self._comm,
            self._params,
            sweeptstamps,
            nsweeps,
            az,
            az_min,
            az_max,
            nces,
            data_size_proc,
            nb_blocks_loc,
            local_blocks_sizes,
            nnz,
            self._mappraiser_pixels,
            self._mappraiser_pixweights,
            self._mappraiser_signal,
            self._mappraiser_noise,
            self._mappraiser_invtt,
            )
        # os.environ["OMP_NUM_THREADS"] = "4"

        return

    def _count_samples(self):
        """ Loop over the observations and count the number of samples.

        """
        if len(self._data.obs) != 1:
            nsamp = 0
            tod0 = self._data.obs[0]["tod"]
            detectors0 = tod0.local_dets
            for obs in self._data.obs:
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
            tod = self._data.obs[0]["tod"]
            nsamp = tod.local_samples[1]
        return nsamp

    def _get_period_ranges(self, detectors):
        """ Collect the ranges of every observation.
        (This routine taken as is from Madam has been truncated, for now it is
        only extracting the frequency binning of the PSDs)
        """
        psdfreqs = None

        for obs in self._data.obs:
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

        return psdfreqs

    def _psd2invtt(self, psdfreqs, psd):
        """ Generate the first rows of the Toeplitz blocks from the PSDs
        """
        # parameters
        sampling_freq = self._params["samplerate"]
        f_defl = sampling_freq/(np.pi*self._params["Lambda"])
        df = f_defl/2
        block_size = 2**(int(math.log(sampling_freq*1./df,2))+1)

        # Invert PSD
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

        # Compute inverse noise autocorrelation functions
        inv_tt = np.real(np.fft.ifft(psdm1, n=block_size))

        # Define apodization window
        window = scipy.signal.gaussian(2*self._params["Lambda"], 1./2*self._params["Lambda"])
        window = np.fft.ifftshift(window)
        window = window[:self._params["Lambda"]]
        window = np.pad(window,(0,int(block_size/2-(self._params["Lambda"]))),'constant')
        symw = np.zeros(block_size)
        symw[:int(block_size/2)] = window
        symw[int(block_size/2):] = np.flip(window,0)

        inv_tt_w = np.multiply(symw, inv_tt, dtype = mappraiser.INVTT_TYPE)

        return inv_tt_w[:self._params["Lambda"]]

    @function_timer
    def _prepare(self):
        """ Examine the data object.

        """
        log = Logger.get()
        timer = Timer()
        timer.start()

        nsamp = self._count_samples()

        # Determine the detectors and the pointing matrix non-zeros
        # from the first observation. Mappraiser will expect these to remain
        # unchanged across observations.

        tod = self._data.obs[0]["tod"]

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
        # N.B: Above comment is from OpMadam, for now
        # it only gives frequency binning of the PSDs
        psdfreqs = self._get_period_ranges(dets)

        self._comm.Barrier()
        if self._rank == 0 and self._verbose:
            timer.report_clear("Collect dataset dimensions")

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
    @function_timer
    def _stage_time(self, detectors, nsamp, psdfreqs):
        """ Stage the timestamps and use them to build PSD inputs.
        N.B: timestamps are not currently used in MAPPRAISER, however, this may
        change in the future. At this stage, the routine builds the time-domain
        Toeplitz blocks inputs (first rows).
        """
        # self._mappraiser_timestamps = self._cache.create(
        #     "timestamps", mappraiser.TIMESTAMP_TYPE, (nsamp,)
        # )

        # offset = 0
        # time_offset = 0
        # psds = {}
        invtt_list = []
        for iobs, obs in enumerate(self._data.obs):
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
                        invtt = self._psd2invtt(psdfreqs, psd)
                        if self._params["map-maker"] == 'ML':
                            invtt_list.append(invtt)
                        else:
                            invtt_list.append(invtt[:1])
                        # if det not in psds:
                        #     psds[det] = [(0, psd)]
                        # else:
                        #     if not np.allclose(psds[det][-1][1], psd):
                        #         psds[det] += [(timestamps[0], psd)]

        return invtt_list

    @function_timer
    def _stage_signal(self, detectors, nsamp, ndet, nodecomm, nread):
        """ Stage signal
        """
        log = Logger.get()
        timer = Timer()
        # Determine if we can purge the signal and avoid keeping two
        # copies in memory
        purge = self._name is not None and self._purge_tod
        if not purge:
            nread = 1
            nodecomm = MPI.COMM_SELF

        for iread in range(nread):
            nodecomm.Barrier()
            timer.start()
            if nodecomm.rank % nread == iread:
                self._mappraiser_signal = self._cache.create(
                "signal", mappraiser.SIGNAL_TYPE, (nsamp * ndet,)
                )
                self._mappraiser_signal[:] = np.nan

                global_offset = 0
                local_blocks_sizes = []
                for iobs, obs in enumerate(self._data.obs):
                    tod = obs["tod"]

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
                    # Purge only after all detectors are staged in case some are aliased
                    # cache.clear() will not fail if the object was already
                    # deleted as an alias
                    if purge:
                        for det in detectors:
                            cachename = "{}_{}".format(self._name, det)
                            tod.cache.clear(cachename)
                    global_offset = offset

                local_blocks_sizes = np.array(local_blocks_sizes, dtype=np.int32)
            if self._verbose and nread > 1:
                nodecomm.Barrier()
                if self._rank == 0:
                    timer.report_clear("Stage signal {} / {}".format(iread + 1, nread))

        return signal_dtype, local_blocks_sizes

    @function_timer
    def _stage_noise(self, detectors, nsamp, ndet, nodecomm, nread):
        """ Stage noise timestream (detector noise + atmosphere)
        """
        log = Logger.get()
        timer = Timer()
        # Determine if we can purge the signal and avoid keeping two
        # copies in memory
        purge = self._noise_name is not None and self._purge_tod
        if not purge:
            nread = 1
            nodecomm = MPI.COMM_SELF

        for iread in range(nread):
            nodecomm.Barrier()
            timer.start()
            if nodecomm.rank % nread == iread:
                self._mappraiser_noise = self._cache.create(
                "noise", mappraiser.SIGNAL_TYPE, (nsamp * ndet,)
                )
                if self._noise_name == None:
                    self._mappraiser_noise = np.zeros_like(self._mappraiser_noise)
                    return self._mappraiser_noise.dtype

                self._mappraiser_noise[:] = np.nan

                global_offset = 0
                for iobs, obs in enumerate(self._data.obs):
                    tod = obs["tod"]

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
                    # Purge only after all detectors are staged in case some are aliased
                    # cache.clear() will not fail if the object was already
                    # deleted as an alias
                    if purge:
                        for det in detectors:
                            cachename = "{}_{}".format(self._noise_name, det)
                            tod.cache.clear(cachename)
                    global_offset = offset
            if self._verbose and nread > 1:
                nodecomm.Barrier()
                if self._rank == 0:
                    timer.report_clear("Stage noise {} / {}".format(iread + 1, nread))

        return noise_dtype

    @function_timer
    def _stage_pixels(self, detectors, nsamp, ndet, nnz, nside):
        """ Stage pixels
        """
        self._mappraiser_pixels = self._cache.create(
            "pixels", mappraiser.PIXEL_TYPE, (nsamp * ndet * nnz,)
        )
        self._mappraiser_pixels[:] = -1

        global_offset = 0
        nces = 0
        sweeptstamps_list = []
        nsweeps_list = []
        az_list = []
        az_min_list = []
        az_max_list = []
        for iobs, obs in enumerate(self._data.obs): #assume only one obs per process for now
            tod = obs["tod"]

            # commonflags = None
            nces += 1
            sweeptstamps = [0]
            commonflags = tod.local_common_flags(self._common_flag_name)
            for iflg, flg in enumerate(commonflags[:-1]):
                if (flg & commonflags[iflg+1] <=1): #sweep direction changes
                    sweeptstamps.append(iflg+1)
            sweeptstamps.append(len(commonflags))
            nsweeps = len(sweeptstamps)-1
            sweeptstamps = np.array(sweeptstamps, dtype=np.int32)

            sweeptstamps_list.append(sweeptstamps)
            nsweeps_list.append(nsweeps)

            qazel = tod.read_boresight_azel()
            az = 180/np.pi *(2*np.pi-qa.to_position(qazel)[1])
            az_list.append(az)
            az_min_list.append(az.min())
            az_max_list.append(az.max())



            for idet, det in enumerate(detectors):
                # Optionally get the flags, otherwise they are
                # assumed to have been applied to the pixel numbers.
                # N.B: MAPPRAISER doesn't use flags for now but might be useful for
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
                    # This is also the case for Mappraiser, could be changed but keeping it for now
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

            # Always purge the pixels but restore them from the Mappraiser
            # buffers when purge_pixels = False
            # Purging MUST happen after all detectors are staged because
            # some the pixel numbers may be aliased between detectors
            for det in detectors:
                pixelsname = "{}_{}".format(self._pixels, det)
                # cache.clear() will not fail if the object was already
                # deleted as an alias
                tod.cache.clear(pixelsname)
                if self._purge_flags and self._flag_name is not None:
                    cacheflagname = "{}_{}".format(self._flag_name, det)
                    tod.cache.clear(cacheflagname)

            del commonflags
            if self._purge_flags and self._common_flag_name is not None:
                tod.cache.clear(self._common_flag_name)
            global_offset = offset

        sweeptstamps_list = np.array(sweeptstamps_list, dtype=np.int32)
        nsweeps_list = np.array(nsweeps_list, dtype=np.int32)
        az_list = np.array(az_list)
        az_min_list = np.array(az_min_list)
        az_max_list = np.array(az_max_list)

        return pixels_dtype, sweeptstamps_list, nsweeps_list, az_list, az_min_list, az_max_list, nces

    @function_timer
    def _stage_pixweights(
        self,
        detectors,
        nsamp,
        ndet,
        nnz,
        nnz_full,
        nnz_stride,
        nodecomm,
        nread,
    ):
        """Now collect the pixel weights
        """
        log = Logger.get()
        timer = Timer()
        # Determine if we can purge the pixel weights and avoid keeping two
        # copies of the weights in memory
        purge = self._purge_weights or (nnz == nnz_full)
        if not purge:
            nread = 1
            nodecomm = MPI.COMM_SELF
        for iread in range(nread):
            nodecomm.Barrier()
            timer.start()
            if nodecomm.rank % nread == iread:
                self._mappraiser_pixweights = self._cache.create(
                "pixweights", mappraiser.WEIGHT_TYPE, (nsamp * ndet * nnz,)
                )
                self._mappraiser_pixweights[:] = 0

                global_offset = 0
                for iobs, obs in enumerate(self._data.obs):
                    tod = obs["tod"]

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
                    # Purge the weights but restore them from the Mappraiser
                    # buffers when purge_weights=False.
                    if purge:
                        for idet, det in enumerate(detectors):
                            weightsname = "{}_{}".format(self._weights, det)
                            tod.cache.clear(pattern=weightsname)

                    global_offset = offset
            if self._verbose and nread > 1:
                nodecomm.Barrier()
                if self._rank == 0:
                    timer.report_clear(
                        "Stage pixel weights {} / {}".format(iread + 1, nread)
                    )
        return weight_dtype

    @function_timer
    def _stage_data(
        self,
        nsamp,
        ndet,
        nnz,
        nnz_full,
        nnz_stride,
        psdfreqs,
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
        log = Logger.get()
        nodecomm = self._comm.Split_type(MPI.COMM_TYPE_SHARED, self._rank)
        # Check if the user has elected to stagger staging the data on each
        # node to avoid exhausting memory
        if self._conserve_memory:
            if self._conserve_memory == 1:
                nread = nodecomm.size
            else:
                nread = min(self._conserve_memory, nodecomm.size)
        else:
            nread = 1

        self._comm.Barrier()
        timer_tot = Timer()
        timer_tot.start()

        # Stage time (Tpltz blocks in Mappraiser), it is never purged
        # so the staging is never stepped
        timer = Timer()
        timer.start()
        invtt_list = self._stage_time(detectors, nsamp, psdfreqs)
        self._mappraiser_invtt = np.array([np.array(invtt_i, dtype= mappraiser.INVTT_TYPE) for invtt_i in invtt_list])
        del invtt_list
        self._mappraiser_invtt = np.concatenate(self._mappraiser_invtt)
        if self._verbose:
            nodecomm.Barrier()
            if self._rank == 0:
                timer.report_clear("Stage time")
        memreport("after staging time", self._comm)  # DEBUG
        count_caches(
            self._data, self._comm, nodecomm, self._cache, "after staging time"
        )  # DEBUG

         # Stage signal.  If signal is not being purged, staging is not stepped
        timer.start()
        signal_dtype, local_blocks_sizes = self._stage_signal(
            detectors, nsamp, ndet, nodecomm, nread
        )
        if self._verbose:
            nodecomm.Barrier()
            if self._rank == 0:
                timer.report_clear("Stage signal")
        memreport("after staging signal", self._comm)  # DEBUG
        count_caches(
            self._data, self._comm, nodecomm, self._cache, "after staging signal"
        )  # DEBUG

        # Stage noise.  If noise is not being purged, staging is not stepped
        timer.start()
        noise_dtype = self._stage_noise(
           detectors, nsamp, ndet, nodecomm, nread
        )
        if self._verbose:
            nodecomm.Barrier()
            if self._rank == 0:
                timer.report_clear("Stage noise")
        memreport("after staging noise", self._comm)  # DEBUG
        count_caches(
            self._data, self._comm, nodecomm, self._cache, "after staging noise"
        )  # DEBUG

        # Stage pixels
        timer_step = Timer()
        timer_step.start()
        for iread in range(nread):
            nodecomm.Barrier()
            timer.start()
            if nodecomm.rank % nread == iread:
                pixels_dtype, sweeptstamps, nsweeps, az, az_min, az_max, nces = self._stage_pixels(
                    detectors, nsamp, ndet, nnz, nside
                )
            if self._verbose and nread > 1:
                nodecomm.Barrier()
                if self._rank == 0:
                    timer.report_clear("Stage pixels {} / {}".format(iread + 1, nread))
        if self._verbose:
            nodecomm.Barrier()
            if self._rank == 0:
                timer_step.report_clear("Stage pixels")
        memreport("after staging pixels", self._comm)  # DEBUG
        count_caches(
            self._data, self._comm, nodecomm, self._cache, "after staging pixels"
        )  # DEBUG

        # Stage pixel weights
        timer_step.start()
        weight_dtype = self._stage_pixweights(
            detectors,
            nsamp,
            ndet,
            nnz,
            nnz_full,
            nnz_stride,
            nodecomm,
            nread,
        )
        if self._verbose:
            nodecomm.Barrier()
            if self._rank == 0:
                timer_step.report_clear("Stage pixel weights")
        memreport("after staging pixel weights", self._comm)  # DEBUG
        count_caches(
            self._data, self._comm, nodecomm, self._cache, "after staging pixel weights"
        )  # DEBUG

        del nodecomm
        if self._rank == 0 and self._verbose:
            timer_tot.report_clear("Stage all data")

        # detweights is either a dictionary of weights specified at
        # construction time, or else we use uniform weighting.
        # N.B: This is essentially useless in current implementation
        detw = {}
        if self._detw is None:
            for idet, det in enumerate(detectors):
                detw[det] = 1.0
        else:
            detw = self._detw

        detweights = np.zeros(ndet, dtype=np.float64)
        for idet, det in enumerate(detectors):
            detweights[idet] = detw[det]

        # Get global array of data sizes of the full communicator
        data_size_proc = np.array(self._comm.allgather(len(self._mappraiser_signal)), dtype=np.int32)
        # Get number of local observations
        nobsloc = len(self._data.obs)

        return data_size_proc, nobsloc, local_blocks_sizes, signal_dtype, noise_dtype, pixels_dtype, sweeptstamps, nsweeps, az, az_min, az_max, nces, weight_dtype

    @function_timer
    def _unstage_signal(self, detectors, nsamp, signal_type):
        # N.B: useful when we want to get back data after mapmaking, not allowed for now
        # if self._name_out is not None:
        #     global_offset = 0
        #     for obs, period_ranges in zip(self._data.obs, obs_period_ranges):
        #         tod = obs["tod"]
        #         nlocal = tod.local_samples[1]
        #         for idet, det in enumerate(detectors):
        #             signal = np.ones(nlocal, dtype=signal_type) * np.nan
        #             offset = global_offset
        #             for istart, istop in period_ranges:
        #                 nn = istop - istart
        #                 dslice = slice(
        #                     idet * nsamp + offset, idet * nsamp + offset + nn
        #                 )
        #                 signal[istart:istop] = self._madam_signal[dslice]
        #                 offset += nn
        #             cachename = "{}_{}".format(self._name_out, det)
        #             tod.cache.put(cachename, signal, replace=True)
        #         global_offset = offset
        self._mappraiser_signal = None
        self._cache.destroy("signal")
        return

    @function_timer
    def _unstage_noise(self, detectors, nsamp, noise_type):
        # N.B: useful when we want to get back data after mapmaking, not allowed for now
        # if self._name_out is not None:
        #     global_offset = 0
        #     for obs, period_ranges in zip(self._data.obs, obs_period_ranges):
        #         tod = obs["tod"]
        #         nlocal = tod.local_samples[1]
        #         for idet, det in enumerate(detectors):
        #             signal = np.ones(nlocal, dtype=signal_type) * np.nan
        #             offset = global_offset
        #             for istart, istop in period_ranges:
        #                 nn = istop - istart
        #                 dslice = slice(
        #                     idet * nsamp + offset, idet * nsamp + offset + nn
        #                 )
        #                 signal[istart:istop] = self._madam_signal[dslice]
        #                 offset += nn
        #             cachename = "{}_{}".format(self._name_out, det)
        #             tod.cache.put(cachename, signal, replace=True)
        #         global_offset = offset
        self._mappraiser_noise = None
        self._cache.destroy("noise")
        return

    @function_timer
    def _unstage_pixels(self, detectors, nsamp, pixels_dtype, nside):
        # N.B: useful when we want to get back data after mapmaking, not allowed for now
        # if not self._purge_pixels:
        #     # restore the pixels from the Madam buffers
        #     global_offset = 0
        #     for obs, period_ranges in zip(self._data.obs, obs_period_ranges):
        #         tod = obs["tod"]
        #         nlocal = tod.local_samples[1]
        #         for idet, det in enumerate(detectors):
        #             pixels = -(np.ones(nlocal, dtype=pixels_dtype))
        #             offset = global_offset
        #             for istart, istop in period_ranges:
        #                 nn = istop - istart
        #                 dslice = slice(
        #                     idet * nsamp + offset, idet * nsamp + offset + nn
        #                 )
        #                 pixels[istart:istop] = self._madam_pixels[dslice]
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
        return

    @function_timer
    def _unstage_pixweights(
        self, detectors, nsamp, weight_dtype, nnz, nnz_full
    ):
        # N.B: useful when we want to get back data after mapmaking, not allowed for now
        # if not self._purge_weights and nnz == nnz_full:
        #     # restore the weights from the Madam buffers
        #     global_offset = 0
        #     for obs, period_ranges in zip(self._data.obs, obs_period_ranges):
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
        #                 weights[istart:istop] = self._madam_pixweights[dwslice].reshape(
        #                     [-1, nnz]
        #                 )
        #                 offset += nn
        #             cachename = "{}_{}".format(self._weights, det)
        #             tod.cache.put(cachename, weights, replace=True)
        #         global_offset = offset
        self._mappraiser_pixweights = None
        self._cache.destroy("pixweights")
        return

    def _unstage_data(
        self,
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
        """ Clear Mappraiser buffers, [restore pointing into TOAST caches-> not done currently].
        """
        log = Logger.get()
        # self._mappraiser_timestamps = None
        # self._cache.destroy("timestamps")

        if self._conserve_memory:
            nodecomm = self._comm.Split_type(MPI.COMM_TYPE_SHARED, self._rank)
            nread = nodecomm.size
        else:
            nodecomm = MPI.COMM_SELF
            nread = 1

        self._comm.Barrier()
        timer_tot = Timer()
        timer_tot.start()
        for iread in range(nread):
            timer_step = Timer()
            timer_step.start()
            timer = Timer()
            timer.start()
            if nodecomm.rank % nread == iread:
                self._unstage_signal(detectors, nsamp, signal_type)
            if self._verbose:
                nodecomm.Barrier()
                if self._rank == 0:
                    timer.report_clear(
                        "Unstage signal {} / {}".format(iread + 1, nread)
                    )
            if nodecomm.rank % nread == iread:
                self._unstage_noise(detectors, nsamp, noise_type)
            if self._verbose:
                nodecomm.Barrier()
                if self._rank == 0:
                    timer.report_clear(
                        "Unstage noise {} / {}".format(iread + 1, nread)
                    )
            if nodecomm.rank % nread == iread:
                self._unstage_pixels(
                    detectors, nsamp, pixels_dtype, nside
                )
            if self._verbose:
                nodecomm.Barrier()
                if self._rank == 0:
                    timer.report_clear(
                        "Unstage pixels {} / {}".format(iread + 1, nread)
                    )
            if nodecomm.rank % nread == iread:
                self._unstage_pixweights(
                    detectors, nsamp, weight_dtype, nnz, nnz_full
                )
            nodecomm.Barrier()
            if self._verbose and self._rank == 0:
                timer.report_clear(
                    "Unstage pixel weights {} / {}".format(iread + 1, nread)
                )
            if self._rank == 0 and self._verbose and nread > 1:
                timer_step.report_clear("Unstage data {} / {}".format(iread + 1, nread))
        self._comm.Barrier()
        if self._rank == 0 and self._verbose:
            timer_tot.report_clear("Unstage all data")

        del nodecomm
        return
