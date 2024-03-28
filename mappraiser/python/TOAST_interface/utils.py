# This script contains a list of routines to set mappraiser parameters, and
# apply the Mappraiser operator during a TOD2MAP TOAST(3) pipeline

# @author: Hamza El Bouhargani
# @date: January 2020

import os
import numpy as np
import scipy.signal
from scipy.optimize import curve_fit
from astropy import units as u

from toast.utils import Logger, dtype_to_aligned, memreport
from toast.ops.memory_counter import MemoryCounter


# Here are some helper functions adapted from toast/src/ops/madam_utils.py
def log_time_memory(
        data, timer=None, timer_msg=None, mem_msg=None, full_mem=False, prefix=""
):
    """(This function is taken from madam_utils.py)"""
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
            timer.clear()

    if mem_msg is not None:
        # Dump toast memory use
        mem_count = MemoryCounter(silent=True)
        mem_count.total_bytes = 0
        toast_bytes = mem_count.apply(data)

        if data.comm.group_rank == 0:
            msg = "{} {} Group {} memory = {:0.2f} GB".format(
                prefix, mem_msg, data.comm.group, toast_bytes / 1024 ** 2
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
        n_repeat=1,
):
    """Helper function to fill a mappraiser buffer from a local detdata key.
    (This function is taken from madam_utils.py)
    """
    # interval_offset = 0
    do_flags = False
    if shared_flags is not None or det_flags is not None:
        do_flags = True
        # Flagging should only be enabled when we are processing the pixel indices
        # (which is how mappraiser effectively implements flagging).  So we will set
        # all flagged samples to "-1" below.
    
    # The line below redefines the offset variable to have dets of a single 
    # observation contiguous in memory
    offset = 0
    
    for ob in data.obs:
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
                # This was used in a deprecated data arrangement:
                # offset = interval_starts[interval_offset + ivw]

                flags = None
                if do_flags:
                    # Using flags
                    flags = np.zeros(view_samples, dtype=np.uint8)
                if shared_flags is not None:
                    flags |= views.shared[shared_flags][ivw] & shared_mask

                # This was used in a deprecated data arrangement:
                # slc = slice(
                #     (idet * nsamp + offset) * nnz,
                #     (idet * nsamp + offset + view_samples) * nnz,
                #     1,
                # )
                slc = slice(
                    (offset) * nnz,
                    (offset + view_samples) * nnz,
                    1,
                )                
                if detdata_name is not None:
                    if nnz > 1:
                        mappraiser_buffer[slc] = np.repeat(
                            views.detdata[detdata_name][ivw][det].flatten()[::nnz_stride],
                            n_repeat,
                        )
                    else:
                        mappraiser_buffer[slc] = np.repeat(
                            views.detdata[detdata_name][ivw][det].flatten(),
                            n_repeat,
                        )
                else:
                    # Noiseless cases (noise_name=None).
                    mappraiser_buffer[slc] = 0.0

                detflags = None
                if do_flags:
                    if det_flags is None:
                        detflags = flags
                    else:
                        detflags = np.copy(flags)
                        detflags |= views.detdata[det_flags][ivw][det] & det_mask
                    # mappraiser's pixels buffer has nnz=3, not nnz=1
                    repeated_flags = np.repeat(detflags, n_repeat)
                    mappraiser_buffer[slc][repeated_flags != 0] = -1
                
                # The line below goes with the new redifinition of the offset variable:
                offset += view_samples
        if do_purge:
            del ob.detdata[detdata_name]
        # interval_offset += len(views)
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
        n_repeat=1,
):
    """When purging data, take turns staging it.
    (This function is taken from madam_utils.py)
    """
    raw = None
    wrapped = None
    for copying in range(n_copy_groups):
        if nodecomm.rank % n_copy_groups == copying:
            # Our turn to copy data
            storage, _ = dtype_to_aligned(mappraiser_dtype)
            raw = storage.zeros(nsamp * len(dets) * nnz)
            wrapped = raw.array()
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
                n_repeat=n_repeat,
            )
        nodecomm.barrier()
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
    """Helper function to create a detdata buffer from mappraiser data.
    (This function is taken from madam_utils.py)
    """
    # The line below redefines the offset variable to have dets of a single 
    # observation contiguous in memory
    offset = 0
    
    # interval = 0
    for ob in data.obs:
        # Create the detector data
        if nnz == 1:
            ob.detdata.create(detdata_name, dtype=detdata_dtype)
        else:
            ob.detdata.create(detdata_name, dtype=detdata_dtype, sample_shape=(nnz,))
        # Loop over detectors
        views = ob.view[view]
        ldet = 0
        for det in dets:
            if det not in ob.local_detectors:
                continue
            # idet = ob.local_detectors.index(det)
            # Loop over views
            for ivw, vw in enumerate(views):
                view_samples = None
                if vw.start is None:
                    # This is a view of the whole obs
                    view_samples = ob.n_local_samples
                else:
                    view_samples = vw.stop - vw.start
                # This was used in a deprecated data arrangement:
                # offset = interval_starts[interval]
                # slc = slice(
                #     (idet * nsamp + offset) * nnz,
                #     (idet * nsamp + offset + view_samples) * nnz,
                #     1,
                # )
                slc = slice(
                    (offset) * nnz,
                    (offset + view_samples) * nnz,
                    1,
                )
                if nnz > 1:
                    views.detdata[detdata_name][ivw][ldet] = mappraiser_buffer[
                        slc
                    ].reshape((-1, nnz))
                else:
                    views.detdata[detdata_name][ivw][ldet] = mappraiser_buffer[slc]
                offset += view_samples
                # This was used in a deprecated data arrangement:
                # interval += 1
            ldet += 1
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
    """When restoring data, take turns copying it.
    (This function is taken from madam_utils.py)
    """
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

def stack_padding(it):
    """Turns a list of arrays of different sizes to a matrix which rows
    are the different arrays padded with zeros such that they all have the same length."""
    # stack padding may be memory inefficient but was found as a solution
    # for the interfacing of double pointers in the C backend. Should remain until we find a better solution.
    def resize(row, size):
        new = np.array(row)
        new.resize(size)
        return new

    # find longest row length
    row_length = max(it, key=len).__len__()
    mat = np.array( [resize(row, row_length) for row in it] )

    return mat

def stage_polymetadata(
    data,
    view,
    params,
    shared_flags,
):
    """Builds a set of arrays of indices that mark the changes in scan direction and/or mark
    the switch from one polynomial baseline to the other (one array per local CES), also computes the number of baselines per CES
    and the local number of CES"""
    fsamp = params["fsample"]
    # remove unit of fsamp to avoid problems when computing periodogram
    try:
        f_unit = fsamp.unit
        fsamp = float(fsamp / (1.0 * f_unit))
    except AttributeError:
        pass
    sweeptstamps_list = []
    nsweeps_list = []
    for ob in data.obs:
        views = ob.view[view]
        for ivw, vw in enumerate(views):
            commonflags = views.shared[shared_flags][ivw]
            sweeptstamps = [0]
            if params["fixed_polybase"]:
                base_sup_id = int(params["polybaseline_length"] * fsamp)
                while base_sup_id < len(commonflags):
                    sweeptstamps.append(base_sup_id)
                    base_sup_id += int(params["polybaseline_length"] * fsamp)
            else:
                for iflg, flg in enumerate(commonflags[:-1]):
                    # detect switch LR/RL scan and viceversa (turnarounds included in sweeps)
                    # second condition was added because I noticed commonflags end with isolated turnaround flags, aka "3"
                    # which in the current implementation will be ignored
                    if ((flg & commonflags[iflg+1] < 8) and commonflags[iflg+1] >= 8):
                        sweeptstamps.append(iflg+1)
            sweeptstamps.append(len(commonflags))
            nsweeps = len(sweeptstamps)-1
            sweeptstamps = np.array(sweeptstamps, dtype=np.int32)

            sweeptstamps_list.append(sweeptstamps)
            nsweeps_list.append(nsweeps)
    
    sweeptstamps_list = stack_padding(sweeptstamps_list)
    nsweeps_list = np.array(nsweeps_list, dtype=np.int32)

    return sweeptstamps_list, nsweeps_list

def stage_azscan(
    data,
    view,
    az_name,
):
    """Stages the boresight azimuth scan for each local CES as well as the min and max values."""
    az_list = []
    az_min_list = []
    az_max_list = []
    for ob in data.obs:
        views = ob.view[view]
        for ivw, vw in enumerate(views):
            az = views.shared[az_name][ivw]
            az_list.append(az)
            az_min_list.append(az.min())
            az_max_list.append(az.max())
    
    az_list = stack_padding(az_list)
    az_min_list = np.array(az_min_list)
    az_max_list = np.array(az_max_list)

    return az_list, az_min_list, az_max_list

def stage_hwpangle(
    data,
    view,
    hwpangle_name,
):
    """Stage hwp angle for each local CES."""
    hwp_angle_list = []
    for ob in data.obs:
        views = ob.view[view]
        for ivw, vw in enumerate(views):
            hwp_angle = views.shared[hwpangle_name][ivw]
            if hwp_angle is not None:
                hwp_angle_list.append(hwp_angle)
    
    hwp_angle_list = stack_padding(hwp_angle_list)

    return hwp_angle_list

def compute_local_block_sizes(data, view, dets, buffer):
    """Compute the sizes of the local data blocks and store them in the provided buffer."""
    for iobs, ob in enumerate(data.obs):
        views = ob.view[view]
        for idet, det in enumerate(dets):
            if det not in ob.local_detectors:
                continue
            # Loop over views
            for vw in views:
                view_samples = None
                if vw.start is None:
                    # This is a view of the whole obs
                    view_samples = ob.n_local_samples
                else:
                    view_samples = vw.stop - vw.start
                buffer[iobs * len(dets) + idet] += view_samples
    return


def compute_invtt(
        nobs,
        ndet,
        mappraiser_noise,
        local_block_sizes,
        Lambda,
        fsamp,
        buffer,
        invtt_dtype,
        print_info=False,
        save_psd=False,
        save_dir=None,
        apod_window_type="chebwin",
):
    """Compute the first lines of the blocks of the banded noise covariance and store them in the provided buffer."""
    offset = 0
    for iobs in range(nobs):
        for idet in range(ndet):
            # # Old data arrangement:
            # blocksize = local_block_sizes[idet * nobs + iobs]
            # nsetod = mappraiser_noise[offset: offset + blocksize]
            # slc = slice(
            #     (idet * nobs + iobs) * Lambda,
            #     (idet * nobs + iobs) * Lambda + Lambda,
            #     1,
            # )    
            blocksize = local_block_sizes[iobs * ndet + idet]            
            nsetod = mappraiser_noise[offset: offset + blocksize]
            slc = slice(
                (iobs * ndet + idet) * Lambda,
                (iobs * ndet + idet) * Lambda + Lambda,
                1,
            )
            buffer[slc] = noise2invtt(
                nsetod,
                blocksize,
                Lambda,
                fsamp,
                idet,
                invtt_dtype,
                apod_window_type,
                verbose=(print_info and (idet == 0) and (iobs == 0)),
                save_psd=(save_psd and (idet == 0) and (iobs == 0)),
                save_dir=save_dir,
            )
            offset += blocksize
    return


def psd_model(f, sigma, alpha, fknee, fmin):
    return sigma * (1 + ((f + fmin) / fknee) ** alpha)


def logpsd_model(f, a, alpha, fknee, fmin):
    return a + np.log10(1 + ((f + fmin) / fknee) ** alpha)


def inversepsd_model(f, sigma, alpha, fknee, fmin):
    return sigma * 1.0 / (1 + ((f + fmin) / fknee) ** alpha)


def inverselogpsd_model(f, a, alpha, fknee, fmin):
    return a - np.log10(1 + ((f + fmin) / fknee) ** alpha)


def noise2invtt(
        nsetod,
        nn,
        Lambda,
        fsamp,
        idet,
        invtt_dtype,
        apod_window_type,
        verbose=False,
        save_psd=False,
        save_dir=None,
):
    """Computes a periodogram from a noise timestream, and fits a PSD model
    to it, which is then used to build the first row of a Toeplitz block.
    """
    # remove unit of fsamp to avoid problems when computing periodogram
    try:
        f_unit = fsamp.unit
        fsamp = float(fsamp / (1.0 * f_unit))
    except AttributeError:
        pass

    # Estimate psd from noise timestream

    # Length of segments used to estimate PSD (defines the lowest frequency we can estimate)
    nperseg = nn

    # Compute a periodogram with Welch's method
    f, psd = scipy.signal.welch(nsetod, fsamp, window='hann', nperseg=nperseg, detrend='linear')

    # Fit the psd model to the periodogram (in log scale)
    popt, pcov = curve_fit(
        logpsd_model,
        f[1:],
        np.log10(psd[1:]),
        p0=np.array([-7, -1.0, 0.1, 1e-6]),
        bounds=([-20, -10, 0.0, 0.0], [0.0, 0.0, 20, 1]),
        maxfev=1000,
    )

    if verbose:
        print(
            "\n[det "
            + str(idet)
            + "]: PSD fit log(sigma2) = %1.2f, alpha = %1.2f, fknee = %1.2f, fmin = %1.2e\n"
            % tuple(popt),
            flush=True,
        )
        print("[det {}]: PSD fit covariance: \n{}\n".format(idet, pcov), flush=True)

    # Initialize full size inverse PSD in frequency domain
    f_full = np.fft.rfftfreq(nn, 1.0 / fsamp)

    # Compute inverse noise psd from fit and extrapolate (if needed) to lowest frequencies
    psdm1 = inversepsd_model(f_full, 10 ** (-popt[0]), *popt[1:])

    # Compute inverse noise autocorrelation functions
    inv_tt = np.fft.irfft(psdm1, n=nn)

    # Define apodization window
    # Only allow max lambda = nn//2
    if Lambda > nn // 2:
        raise RuntimeError("Bandwidth cannot be larger than timestream.")

    if apod_window_type == "gaussian":
        q_apo = 3  # Apodization factor: cut happens at q_apo * sigma in the Gaussian window
        window = scipy.signal.get_window(
            ("general_gaussian", 1, 1 / q_apo * Lambda), 2 * Lambda
        )
    elif apod_window_type == "chebwin":
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.chebwin.html#scipy.signal.windows.chebwin
        at = 150
        window = scipy.signal.get_window(("chebwin", at), 2 * Lambda)
    else:
        raise RuntimeError(f"Apodisation window '{apod_window_type}' is not supported.")

    window = np.fft.ifftshift(window)
    window = window[:Lambda]

    # Apply window
    inv_tt_w = np.multiply(window, inv_tt[:Lambda], dtype=invtt_dtype)

    # Optionally save some PSDs for plots
    if save_psd:
        # Save TOD
        np.save(os.path.join(save_dir, "tod.npy"), nsetod)

        # save simulated PSD
        np.save(os.path.join(save_dir, "f_psd.npy"), f)
        np.save(os.path.join(save_dir, "psd_sim.npy"), psd)
        # ipsd = np.reciprocal(psd)
        # np.save(os.path.join(save_dir, "ipsd_sim.npy"), ipsd)

        # save fit of the PSD
        np.save(os.path.join(save_dir, "f_full.npy"), f_full[: nn // 2])
        # np.save(os.path.join(save_dir, "ipsd_fit.npy"), ipsd_fit[: nn // 2])

        # save effective PSDs
        ipsd_eff = compute_psd_eff(inv_tt_w, nn)
        np.save(os.path.join(save_dir, "ipsd_eff" + str(Lambda) + ".npy"), ipsd_eff)

    return inv_tt_w


def compute_psd_eff(t, n):
    # Compute "effective" PSD from truncated autocorrelation function
    lag = t.size
    circ_t = np.pad(t, (0, n - lag), "constant")
    if lag > 1:
        circ_t[-lag + 1:] = np.flip(t[1:], 0)
    psd_eff = np.abs(np.fft.fft(circ_t, n=n))
    return psd_eff[: n // 2]
