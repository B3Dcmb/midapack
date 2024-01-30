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
        pair_diff=False,
        pair_skip=False,
        select_qu=False,
):
    """Helper function to fill a mappraiser buffer from a local detdata key.
    (This function is taken from madam_utils.py)
    """
    interval_offset = 0
    do_flags = False
    if shared_flags is not None or det_flags is not None:
        do_flags = True
        # Flagging should only be enabled when we are processing the pixel indices
        # (which is how mappraiser effectively implements flagging).  So we will set
        # all flagged samples to "-1" below.
    if pair_diff and pair_skip:
        raise RuntimeError("pair_diff and pair_skip in stage_local are incompatible.")

    for ob in data.obs:
        views = ob.view[view]
        if pair_diff or pair_skip:
            for idet in range(0, len(dets) - 1, 2):
                det_a = dets[idet]
                det_b = dets[idet + 1]
                det_0_local = det_a in ob.local_detectors
                det_1_local = det_b in ob.local_detectors
                if not (det_0_local and det_1_local):
                    if det_0_local:
                        msg = "det_b is not in local detectors\n"
                    else:
                        msg = "det_a is not in local detectors\n"
                    msg += f"ob = {ob}\n"
                    msg += f"idet = {idet}\n"
                    msg += f"det_a = {det_a}\n"
                    msg += f"det_b = {det_b}"
                    raise RuntimeError(msg)
                if operator is not None:
                    # Synthesize data for staging
                    obs_data = data.select(obs_uid=ob.uid)
                    operator.apply(obs_data, detectors=[det_a, det_b])
                # Loop over views
                for ivw, vw in enumerate(views):
                    view_samples = None
                    if vw.start is None:
                        # This is a view of the whole obs
                        view_samples = ob.n_local_samples
                    else:
                        view_samples = vw.stop - vw.start
                    offset = interval_starts[interval_offset + ivw]
                    flags = None
                    if do_flags:
                        # Using flags
                        flags = np.zeros(view_samples, dtype=np.uint8)
                    if shared_flags is not None:
                        flags |= views.shared[shared_flags][ivw] & shared_mask

                    slc = slice(
                        ((idet // 2) * nsamp + offset) * nnz,
                        ((idet // 2) * nsamp + offset + view_samples) * nnz,
                        1,
                    )
                    if detdata_name is not None:
                        if nnz > 1:
                            if select_qu:
                                mappraiser_buffer[slc] = np.repeat(
                                    views.detdata[detdata_name][ivw][det_a][..., 1:].flatten()[::nnz_stride],
                                    n_repeat,
                                )
                            else:
                                mappraiser_buffer[slc] = np.repeat(
                                    views.detdata[detdata_name][ivw][det_a].flatten()[::nnz_stride],
                                    n_repeat,
                                )
                            if pair_diff:
                                # We are staging signal or noise
                                # Take the half difference
                                mappraiser_buffer[slc] = 0.5 * (mappraiser_buffer[slc] - np.repeat(
                                    views.detdata[detdata_name][ivw][det_b].flatten()[::nnz_stride],
                                    n_repeat,
                                ))
                        else:
                            mappraiser_buffer[slc] = np.repeat(
                                views.detdata[detdata_name][ivw][det_a].flatten(),
                                n_repeat,
                            )
                            if pair_diff:
                                mappraiser_buffer[slc] = 0.5 * (mappraiser_buffer[slc] - np.repeat(
                                    views.detdata[detdata_name][ivw][det_b].flatten(),
                                    n_repeat,
                                ))
                    else:
                        # Noiseless cases (noise_name=None).
                        mappraiser_buffer[slc] = 0.0

                    if do_flags:
                        if det_flags is None:
                            detflags = flags
                        else:
                            detflags = np.copy(flags)
                            detflags |= views.detdata[det_flags][ivw][det_a] & det_mask
                            detflags |= views.detdata[det_flags][ivw][det_b] & det_mask
                        # mappraiser's pixels buffer has nnz=3, not nnz=1
                        repeated_flags = np.repeat(detflags, n_repeat)
                        mappraiser_buffer[slc][repeated_flags != 0] = -1
        else:
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

                    if do_flags:
                        if det_flags is None:
                            detflags = flags
                        else:
                            detflags = np.copy(flags)
                            detflags |= views.detdata[det_flags][ivw][det] & det_mask
                        # mappraiser's pixels buffer has nnz=3, not nnz=1
                        repeated_flags = np.repeat(detflags, n_repeat)
                        mappraiser_buffer[slc][repeated_flags != 0] = -1
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
        n_repeat=1,
        pair_diff=False,
        pair_skip=False,
        select_qu=False,
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
            if pair_diff:
                raw = storage.zeros(nsamp * (len(dets) // 2) * nnz)
            else:
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
                pair_diff=pair_diff,
                pair_skip=pair_skip,
                select_qu=select_qu,
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


def apo_window(lambd: int, kind="chebwin") -> np.ndarray:
    if kind == "gaussian":
        q_apo = 3  # Apodization factor: cut happens at q_apo * sigma in the Gaussian window
        window = scipy.signal.get_window(
            ("general_gaussian", 1, 1 / q_apo * lambd), 2 * lambd
        )
    elif kind == "chebwin":
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.windows.chebwin.html#scipy.signal.windows.chebwin
        at = 150
        window = scipy.signal.get_window(("chebwin", at), 2 * lambd)
    else:
        raise RuntimeError(f"Apodisation window '{kind}' is not supported.")

    return np.fft.ifftshift(window)[:lambd]


def compute_autocorrelations(
        nobs,
        ndet,
        mappraiser_noise,
        local_block_sizes,
        lambda_,
        fsamp,
        buffer_inv_tt,
        buffer_tt,
        invtt_dtype,
        print_info=False,
        save_psd=False,
        save_dir="",
        apod_window_type="chebwin",
):
    """Compute the first lines of the blocks of the banded noise covariance and store them in the provided buffer."""
    offset = 0
    for iobs in range(nobs):
        for idet in range(ndet):
            blocksize = local_block_sizes[idet * nobs + iobs]
            nsetod = mappraiser_noise[offset: offset + blocksize]
            slc = slice(
                (idet * nobs + iobs) * lambda_,
                (idet * nobs + iobs) * lambda_ + lambda_,
                1,
            )
            buffer_inv_tt[slc], _ = noise_autocorrelation(
                nsetod,
                blocksize,
                lambda_,
                fsamp,
                idet,
                invtt_dtype,
                apod_window_type,
                verbose=(print_info and (idet == 0) and (iobs == 0)),
                save_psd=(save_psd and (idet == 0) and (iobs == 0)),
                save_dir=save_dir,
            )
            buffer_tt[slc] = compute_autocorr(1 / compute_psd_eff(buffer_inv_tt[slc], blocksize), lambda_)
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


def noise_autocorrelation(
        nsetod,
        nn,
        lambda_,
        fsamp,
        idet,
        invtt_dtype,
        apod_window_type,
        nperseg=0,
        verbose=False,
        save_psd=False,
        save_dir="",
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
    if nperseg == 0:
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
    ipsd_fit = inversepsd_model(f_full, 10 ** (-popt[0]), *popt[1:])
    psd_fit = psd_model(f_full, 10 ** (popt[0]), *popt[1:])
    psd_fit[0] = 0.0

    # Compute inverse noise auto-correlation functions
    inv_tt = np.fft.irfft(ipsd_fit, n=nn)
    tt = np.fft.irfft(psd_fit, n=nn)

    # Define apodization window
    # Only allow max lambda = nn//2
    if lambda_ > nn // 2:
        raise RuntimeError("Bandwidth cannot be larger than timestream.")

    window = apo_window(lambda_, kind=apod_window_type)

    # Apply window
    inv_tt_w = np.multiply(window, inv_tt[:lambda_], dtype=invtt_dtype)
    tt_w = np.multiply(window, tt[:lambda_], dtype=invtt_dtype)

    # Keep the same norm
    #
    # rescale_tt = np.sqrt(np.sum(np.square(tt)) / np.sum(np.square(tt_w)))
    # rescale_itt = np.sqrt(np.sum(np.square(inv_tt)) / np.sum(np.square(inv_tt_w)))
    #
    # print("correction for     tt = ", rescale_tt)
    # print("correction for inv_tt = ", rescale_itt)
    #
    # tt_w *= rescale_tt
    # inv_tt_w *= rescale_itt

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
        np.save(os.path.join(save_dir, "psd_fit.npy"), psd_fit[: nn // 2])
        # np.save(os.path.join(save_dir, "ipsd_fit.npy"), ipsd_fit[: nn // 2])

        # save effective PSDs
        ipsd_eff = compute_psd_eff(inv_tt_w, nn)
        psd_eff = compute_psd_eff(tt_w, nn)
        np.save(os.path.join(save_dir, "ipsd_eff" + str(lambda_) + ".npy"), ipsd_eff)
        np.save(os.path.join(save_dir, "psd_eff" + str(lambda_) + ".npy"), psd_eff)

    return inv_tt_w, tt_w


def compute_psd_eff(tt, m) -> np.ndarray:
    """
    Computes the power spectral density from a given autocorrelation function.

    :param tt: Input autocorrelation
    :param m: FFT size

    :return: The PSD (size = m // 2 + 1 beacuse we are using np.fft.rfft)
    """
    # Form the cylic autocorrelation
    lag = len(tt)
    circ_t = np.pad(tt, (0, m - lag), "constant")
    if lag > 1:
        circ_t[-lag + 1:] = np.flip(tt[1:], 0)

    # FFT
    psd = np.fft.rfft(circ_t, n=m)

    return np.real(psd)


def compute_autocorr(psd, lambd: int, apo=True) -> np.ndarray:
    """
    Computes the autocorrelation function from a given power spectral density.

    :param psd: Input PSD
    :param lambd: Assumed noise correlation length
    :param apo: if True, apodize the autocorrelation function

    :return: The autocorrelation function apodised and cut after `lambd` terms
    """

    # Compute the inverse FFT
    autocorr = np.fft.irfft(psd)[:lambd]

    # Apodisation
    if apo:
        return np.multiply(apo_window(lambd), autocorr)
    else:
        return autocorr
