/**
 * @file gap_filling.cpp
 * @brief Implementation of gap_filling routines
 * @author Simon Biquard
 * @date March 2023
 */

#include <algorithm>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <numeric>

#include <mappraiser/RunningSum.h>
#include <mappraiser/create_toeplitz.h>
#include <mappraiser/gap_filling.h>
#include <mappraiser/mapping.h>
#include <mappraiser/rng.h>
#include <mappraiser/stopwatch.h>
#include <mappraiser/weight.h>

// ____________________________________________________________
// GapFillRecap

mappraiser::GapFillRecap::GapFillRecap(const mappraiser::GapFillInfo &info) {
    n_blocks = info.n_blocks;

    mean_iterations = info.get_mean_iterations();
    mean_time = info.get_mean_seconds();

    const auto minmax_iter = std::minmax_element(info.nb_iterations.begin(),
                                                 info.nb_iterations.end());
    const auto minmax_seconds =
        std::minmax_element(info.block_times.begin(), info.block_times.end());

    min_iter = *minmax_iter.first;
    max_iter = *minmax_iter.second;
    min_iter_idx = static_cast<int>(
        std::distance(info.nb_iterations.begin(), minmax_iter.first));
    max_iter_idx = static_cast<int>(
        std::distance(info.nb_iterations.begin(), minmax_iter.second));

    min_time = *minmax_seconds.first;
    max_time = *minmax_seconds.second;
    min_time_idx = static_cast<int>(
        std::distance(info.block_times.begin(), minmax_seconds.first));
    max_time_idx = static_cast<int>(
        std::distance(info.block_times.begin(), minmax_seconds.second));
}

void mappraiser::GapFillRecap::send(MPI_Comm comm, int dest) const {
    MPI_Send(&n_blocks, 1, MPI_INT, dest, 10, comm);
    MPI_Send(&mean_iterations, 1, MPI_DOUBLE, dest, 0, comm);
    MPI_Send(&mean_time, 1, MPI_DOUBLE, dest, 1, comm);
    MPI_Send(&min_iter, 1, MPI_INT, dest, 2, comm);
    MPI_Send(&max_iter, 1, MPI_INT, dest, 3, comm);
    MPI_Send(&min_iter_idx, 1, MPI_INT, dest, 4, comm);
    MPI_Send(&max_iter_idx, 1, MPI_INT, dest, 5, comm);
    MPI_Send(&min_time, 1, MPI_DOUBLE, dest, 6, comm);
    MPI_Send(&max_time, 1, MPI_DOUBLE, dest, 7, comm);
    MPI_Send(&min_time_idx, 1, MPI_DOUBLE, dest, 8, comm);
    MPI_Send(&max_time_idx, 1, MPI_DOUBLE, dest, 9, comm);
}

void mappraiser::GapFillRecap::receive(MPI_Comm comm, int src) {
    MPI_Recv(&n_blocks, 1, MPI_INT, src, 10, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&mean_iterations, 1, MPI_DOUBLE, src, 0, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&mean_time, 1, MPI_DOUBLE, src, 1, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&min_iter, 1, MPI_INT, src, 2, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&max_iter, 1, MPI_INT, src, 3, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&min_iter_idx, 1, MPI_INT, src, 4, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&max_iter_idx, 1, MPI_INT, src, 5, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&min_time, 1, MPI_DOUBLE, src, 6, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&max_time, 1, MPI_DOUBLE, src, 7, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&min_time_idx, 1, MPI_DOUBLE, src, 8, comm, MPI_STATUS_IGNORE);
    MPI_Recv(&max_time_idx, 1, MPI_DOUBLE, src, 9, comm, MPI_STATUS_IGNORE);
}

void mappraiser::GapFillRecap::print(std::ostream &out) const {
    out << "number of blocks = " << n_blocks << "\n"
        << "  -> iter/block = " << mean_iterations << "\n"
        << "            min = " << min_iter << " (block " << min_iter_idx
        << ")\n"
        << "            max = " << max_iter << " (block " << max_iter_idx
        << ")\n"
        << "  -> time/block = " << mean_time << " s\n"
        << "            min = " << min_time << " s (block " << min_time_idx
        << ")\n"
        << "            max = " << max_time << " s (block " << max_time_idx
        << ")\n";
}

std::ostream &mappraiser::operator<<(std::ostream &out,
                                     const mappraiser::GapFillRecap &recap) {
    recap.print(out);
    return out;
}

// ____________________________________________________________
// GapFillInfo

mappraiser::GapFillInfo::GapFillInfo(int blocks, int gaps, int rank)
    : n_gaps(gaps), n_blocks(blocks), current_block(0), current_size(0),
      id(rank) {
    nb_iterations.resize(n_blocks);
    block_times.resize(n_blocks);
    pcg_times.resize(n_blocks);
    valid_fracs.resize(n_blocks);
}

double mappraiser::GapFillInfo::get_mean_iterations() const {
    double sum =
        std::accumulate(nb_iterations.begin(), nb_iterations.end(), 0.0);
    return sum / static_cast<double>(nb_iterations.size());
}

double mappraiser::GapFillInfo::get_mean_seconds() const {
    double sum = std::accumulate(block_times.begin(), block_times.end(), 0.0);
    return sum / static_cast<double>(block_times.size());
}

void mappraiser::GapFillInfo::print_curr_block() const {
    std::cout << "[id " << id << "] block = " << current_block << "/"
              << n_blocks - 1 << "; ";
    std::cout << "size = " << current_size << "; ";
    std::cout << "time = " << block_times[current_block] << " s; ";
    std::cout << "iterations = " << nb_iterations[current_block];
    std::cout << " (" << pcg_times[current_block] / nb_iterations[current_block]
              << " s / iter); ";
    std::cout << "valid fraction = " << valid_fracs[current_block] << " %\n";
}

// ____________________________________________________________
// Routines

void psd_from_tt(int fftlen, int lambda, int psdlen, const double *tt,
                 double *psd) {
    // FFTW variables
    double *circ_t;
    fftw_complex *cplx_psd;
    fftw_plan p;

    // Allocation
    circ_t = fftw_alloc_real(fftlen);
    cplx_psd = fftw_alloc_complex(psdlen);

    // FFTW plan
    p = fftw_plan_dft_r2c_1d(fftlen, circ_t, cplx_psd, FFTW_ESTIMATE);

    // Initialization of input
    std::copy(tt, (tt + lambda), circ_t);
    std::fill((circ_t + lambda), (circ_t + fftlen - lambda + 1), 0);
    std::reverse_copy(tt + 1, (tt + lambda), (circ_t + fftlen - lambda + 1));

    // Execute FFT plan
    fftw_execute(p);

    // Compute the PSD values
#pragma omp parallel for default(none) shared(psdlen, psd, cplx_psd)           \
    schedule(static)
    for (int i = 0; i < psdlen; ++i) {
        psd[i] = std::abs(std::complex<double>(cplx_psd[i][0], cplx_psd[i][1]));
    }

    // Zero out DC value
    psd[0] = 0;

    // Free allocated memory
    fftw_free(circ_t);
    fftw_free(cplx_psd);
    fftw_destroy_plan(p);
}

int mappraiser::find_valid_samples(Gap *gaps, size_t id0,
                                   std::vector<uint8_t> &valid) {
    int i_gap = -1;
    int n_good = 0;

    for (size_t i = 0; i < valid.size(); ++i) {
        // find closest gap before (or at) current sample
        while (i_gap + 1 < gaps->ngap &&
               gaps->id0gap[i_gap + 1] <= static_cast<long>(id0 + i)) {
            ++i_gap;
        }

        if (i_gap == -1)
            // the sample is valid because all gaps are strictly after
            valid[i] = 1;
        else
            // the sample is valid if the gap is too short
            valid[i] = (gaps->id0gap[i_gap] + gaps->lgap[i_gap] <=
                        static_cast<long>(id0 + i))
                           ? 1
                           : 0;

        if (valid[i])
            ++n_good;
    }
    return n_good;
}

void flagged_running_average(int samples, const double *buf, double *baseline,
                             const uint8_t *valid, int w0) {
    RunningSum<double, uint8_t> sum(w0);
    for (int i = 0; i < samples; ++i) {
        sum.process_index(samples, i, buf, valid);
        baseline[i] = sum.baseline_value();
    }
}

__attribute__((unused)) void running_average(int samples, const double *buf,
                                             double *baseline,
                                             const uint8_t *valid,
                                             int bandwidth) {
    // half size of window
    int w0 = bandwidth / 2;

#pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < samples; ++i) {
        double avg = 0;
        int count = 0;
        int start = i - w0;
        int end = i + w0 - 1;

        // compute the moving average of (valid) samples over [i-w0, i+w0]
#pragma omp simd
        for (int j = start; j < end + 1; ++j) {
            if (0 <= j && j < samples && valid[j]) {
                avg += buf[j];
                ++count;
            }
        }

        // handle case where we don't have enough valid samples to average over
        if (count < bandwidth / 10) {
            int w1 = w0 / 10;
            for (int j = -w1; j < -w0; ++j) {
                if (-1 < j + i && j + i < samples && valid[j + i]) {
                    avg += buf[j + i];
                    ++count;
                }
            }
            for (int j = w0; j < w1; ++j) {
                if (-1 < j + i && j + i < samples && valid[j + i]) {
                    avg += buf[j + i];
                    ++count;
                }
            }
        }

        baseline[i] = avg / static_cast<double>(count);
    }
}

void remove_baseline(int samples, double *buf, double *baseline,
                     const uint8_t *valid, int w0, bool rm) {
    // compute the baseline
    flagged_running_average(samples, buf, baseline, valid, w0);

    // remove the baseline in valid intervals
    if (rm) {
#pragma omp parallel for default(none) shared(samples, baseline, buf, valid)   \
    schedule(static)
        for (int i = 0; i < samples; ++i) {
            // if (valid[i]) { buf[i] -= baseline[i]; }
            buf[i] -= baseline[i];
        }
    }
}

#if 0
template <typename T> void printVector(const std::vector<T> &vec) {
    for (const auto &element : vec) {
        std::cout << +element << ", ";
    }
    std::cout << "\n";
}

template <typename T>
void printVector(const std::vector<T> &vec, size_t n_to_print,
                 size_t offset = 0) {
    if (offset + n_to_print > vec.size()) {
        offset = 0;
        n_to_print = vec.size();
    }
    for (size_t i = offset; i < offset + n_to_print; ++i) {
        std::cout << vec[i] << ", ";
    }
    std::cout << "\n";
}

template <typename T> void printCArray(const T *array, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        std::cout << +array[i] << ", ";
    }
    std::cout << "\n";
}
#endif

void sim_noise_tod(int samples, int lambda, const double *tt, double *buf,
                   uint64_t realization, uint64_t detindx, uint64_t obsindx,
                   uint64_t telescope, double sample_rate) {
    // Logical size of the fft
    int fftlen = 2;
    while (fftlen <= samples)
        fftlen *= 2;

    int npsd = (fftlen / 2) + 1;

    // Normalization
    double norm = sample_rate * static_cast<double>(npsd - 1);

    // FFTW variables
    double *tdata;
    fftw_complex *fdata;
    fftw_plan p;

    // Allocate the input/output buffer
    tdata = fftw_alloc_real(fftlen);
    fdata = fftw_alloc_complex(npsd);

    // Create the FFTW plan
    p = fftw_plan_dft_c2r_1d(fftlen, fdata, tdata, FFTW_ESTIMATE);

    // Generate Re/Im gaussian randoms in a complex array
    uint64_t key1 = realization * 4294967296 + telescope * 65536;
    uint64_t key2 = obsindx * 4294967296 + detindx;
    uint64_t counter1 = static_cast<uint64_t>(lambda) * 4294967296;
    uint64_t counter2 = 0; // incremented in loop during generation

    std::vector<double> rngdata(fftlen);
    mappraiser::rng_dist_normal(fftlen, key1, key2, counter1, counter2,
                                rngdata.data());

    // DC and Nyquist frequency are purely real
    fdata[0][0] = rngdata[0];
    fdata[0][1] = 0.0;
    fdata[npsd - 1][0] = rngdata[npsd - 1];
    fdata[npsd - 1][1] = 0.0;

    // Repack the other values
#pragma omp parallel for default(none) shared(fftlen, fdata, rngdata)          \
    schedule(static)
    for (int i = 1; i < (fftlen / 2); ++i) {
        fdata[i][0] = rngdata[i];
        fdata[i][1] = rngdata[fftlen - i];
    }

    // Compute PSD values (interpolated at the correct frequencies) from the
    // autocorrelation function
    std::vector<double> psd(npsd);
    psd_from_tt(fftlen, lambda, npsd, tt, psd.data());

    // Multiply our random data by the PSD values
#pragma omp parallel for default(none) shared(psd, norm, fdata) schedule(static)
    for (size_t i = 0; i < psd.size(); ++i) {
        double scale = std::sqrt(psd[i] * norm);
        fdata[i][0] *= scale;
        fdata[i][1] *= scale;
    }

    // Execute the FFT plan
    fftw_execute(p);

    // Backward FFT: 1/N factor not included by FFTW
#pragma omp parallel for default(none) shared(fftlen, tdata) schedule(static)
    for (int i = 0; i < fftlen; ++i) {
        tdata[i] /= fftlen;
    }

    // Copy as many samples as we need into our noise vector
    int offset = (fftlen - samples) / 2;
    std::copy((tdata + offset), (tdata + offset + samples), buf);

    // Subtract the DC level
    double DC = 0.0;
#pragma omp parallel for default(none) shared(samples, buf) schedule(static)   \
    reduction(+ : DC)
    for (int i = 0; i < samples; ++i) {
        DC += buf[i];
    }

    DC /= static_cast<double>(samples);

#pragma omp parallel for default(none) shared(samples, DC, buf) schedule(static)
    for (int i = 0; i < samples; ++i) {
        buf[i] -= DC;
    }

    // Free allocated memory
    fftw_free(fdata);
    fftw_free(tdata);
    fftw_destroy_plan(p);
}

void mappraiser::sim_constrained_noise_block(
    mappraiser::GapFillInfo &gfi, Tpltz *N_block, Tpltz *Nm1_block,
    double *noise, Gap *gaps, int w0, uint64_t realization, uint64_t detindx,
    uint64_t obsindx, uint64_t telescope, double sample_rate) {
    // get the number of samples, the global first index and the bandwidth
    const int samples = N_block->tpltzblocks[0].n;
    const auto id0 = static_cast<size_t>(N_block->tpltzblocks[0].idv);
    const int lambda = N_block->tpltzblocks[0].lambda;

    // copy the original noise vector
    std::vector<double> rhs(samples);
    std::copy(noise, (noise + samples), rhs.begin());

    // locate the valid samples
    std::vector<uint8_t> valid(samples);
    int n_good = mappraiser::find_valid_samples(gaps, id0, valid);

    gfi.store_valid_frac(std::ceil(100 * n_good / samples));

    // remove baseline (moving average)
#if 0
    mappraiser::system_stopwatch baseline_watch;
#endif
    std::vector<double> baseline(samples);
    if (w0 < 0)
        w0 = lambda / 2;
    remove_baseline(samples, rhs.data(), baseline.data(), valid.data(), w0,
                    true);
#if 0
    if (gfi.id == 0) {
        double milliseconds =
            baseline_watch.elapsed_time<double, std::chrono::milliseconds>();
        std::cout << "Baseline time (samples: " << samples
                  << ", lambda: " << lambda << ") -> " << milliseconds * 0.001
                  << " s\n";
    }
#endif

    // generate random noise realization "xi" with correlations
    std::vector<double> xi(samples);
    sim_noise_tod(samples, lambda, N_block->tpltzblocks[0].T_block, xi.data(),
                  realization, detindx, obsindx, telescope, sample_rate);

    // rhs = noise - xi
#pragma omp parallel for default(none) shared(samples, rhs, xi) schedule(static)
    for (int i = 0; i < samples; ++i) {
        rhs[i] -= xi[i];
    }

    // invert the system N x = (noise - xi)
    mappraiser::system_stopwatch stopwatch;

    WeightMatrix W_block = createWeightMatrix(Nm1_block, N_block, gaps, ITER);
    int nb_iterations = applyWeightMatrix(&W_block, rhs.data());

    gfi.store_pcg_time(
        stopwatch.elapsed_time<double, std::chrono::milliseconds>());
    gfi.store_nb_iterations(nb_iterations);

    // leaving this piece of code there for the moment (you never know)
    if (nb_iterations == 0) {
        std::cout << "[" << gfi.id << "] BLOCK WITH 0 ITER\n";
        gfi.print_curr_block();
        std::cout << "samples = " << samples << "    id0 = " << id0
                  << "    lambda = " << lambda << "\n";
        std::cout << "first gap: id0 = "
                  << gaps->id0gap[N_block->tpltzblocks[0].first_gap]
                  << "    lgap = "
                  << gaps->lgap[N_block->tpltzblocks[0].first_gap] << "\n";
        std::cout << "last gap: id0 = "
                  << gaps->id0gap[N_block->tpltzblocks[0].last_gap]
                  << "    lgap = "
                  << gaps->lgap[N_block->tpltzblocks[0].last_gap] << "\n";
    }

    // compute the constrained realization
    std::vector<double> constrained(samples);
    std::copy(rhs.begin(), rhs.end(), constrained.begin());
    stbmmProd(N_block, constrained.data());

#pragma omp parallel for default(none)                                         \
    shared(samples, valid, baseline, xi, constrained) schedule(static)
    for (int i = 0; i < samples; ++i) {
        if (valid[i]) {
            // add the baseline back
            constrained[i] += xi[i] + baseline[i];
        } else {
            constrained[i] += xi[i];
            // constrained[i] += xi[i] + baseline[i];
        }
    }

    // TODO check quality of constrained realization against original vector?
    // compute max(|constrained - noise|²)
    double max_err = -1.0;
    for (int i = 0; i < samples; ++i) {
        if (valid[i]) {
            double err =
                (constrained[i] - noise[i]) * (constrained[i] - noise[i]);
            if (err > max_err) {
                max_err = err;
            }
        }
    }

    double ratio = max_err / N_block->tpltzblocks[0].T_block[0];

    if (ratio > 1e-6)
        std::cout << "max sq err / sigma² = " << max_err << " / "
                  << N_block->tpltzblocks[0].T_block[0] << " = " << ratio
                  << "\n";

    // copy final result into noise vector
    // std::copy(constrained.begin(), constrained.end(), noise);
    // don't modify the valid samples
    for (int i = 0; i < samples; ++i) {
        if (!valid[i]) {
            noise[i] = constrained[i];
        }
    }
}

void mappraiser::sim_constrained_noise(
    mappraiser::GapFillInfo &gfi, Tpltz *N, Tpltz *Nm1, double *noise,
    Gap *gaps, uint64_t realization, const uint64_t *detindxs,
    const uint64_t *obsindxs, const uint64_t *telescopes, double sample_rate) {
    if (gfi.n_gaps == 0) {
        // nothing to do
        return;
    } else {
        // Loop through toeplitz blocks
        int t_id = 0;
        double *noise_block;
        Tpltz N_block, Nm1_block;

        for (int i = 0; i < N->nb_blocks_loc; ++i) {
            // define single-block Tpltz structures
            set_tpltz_struct(&N_block, N, &N->tpltzblocks[i]);
            set_tpltz_struct(&Nm1_block, Nm1, &Nm1->tpltzblocks[i]);

            // pointer to current block in the tod
            noise_block = (noise + t_id);

            // compute a constrained noise realization for this block
            gfi.set_current_block(i);
            gfi.set_current_size(N_block.local_V_size);
            mappraiser::system_stopwatch stopwatch;
            mappraiser::sim_constrained_noise_block(
                gfi, &N_block, &Nm1_block, noise_block, gaps, -1, realization,
                detindxs[i], obsindxs[i], telescopes[i], sample_rate);
            gfi.store_block_time(
                stopwatch.elapsed_time<double, std::chrono::milliseconds>());

            t_id += N->tpltzblocks[i].n;

#ifdef DEBUG
            if (i % 50 == 0) {
                gfi.print_curr_block();
            }
#endif
        }
    }
}

void perform_gap_filling(MPI_Comm comm, Tpltz *N, Tpltz *Nm1, double *noise,
                         Gap *gaps, uint64_t realization,
                         const uint64_t *detindxs, const uint64_t *obsindxs,
                         const uint64_t *telescopes, double sample_rate,
                         bool verbose) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (rank == 0) {
        std::cout << "[Gap filling] begin routine\n";
    }

    double st = MPI_Wtime();

    mappraiser::GapFillInfo gfi(N->nb_blocks_loc, gaps->ngap, rank);
    mappraiser::sim_constrained_noise(gfi, N, Nm1, noise, gaps, realization,
                                      detindxs, obsindxs, telescopes,
                                      sample_rate);

    MPI_Barrier(comm);
    double t = MPI_Wtime();

    if (rank == 0) {
        std::cout << "[Gap filling] finished in " << t - st
                  << " s (all procs)\n";
    }

    if (verbose) {
        // all procs compute their own recap
        mappraiser::GapFillRecap gf_recap(gfi);
        if (rank != 0) {
            // send recap to proc 0
            gf_recap.send(comm, 0);
        } else {
            // let proc 0 do the printing
            std::cout << "[Gap filling] recap from all processes\n";
            for (int i = 0; i < size; ++i) {
                if (i != 0) {
                    // reveice informations from proc i
                    gf_recap.receive(comm, i);
                }
                // print proc i recap
                std::cout << "[proc " << i << "] " << gf_recap;
            }
            std::cout << "\n";
        }
    }
}

void gap_filling(MPI_Comm comm, const int *data_size_proc, int nb_blocks_loc,
                 int *local_blocks_sizes, int nnz, double *tt, double *inv_tt,
                 int lambda, double *noise, int *indices, uint64_t realization,
                 const uint64_t *detindxs, const uint64_t *obsindxs,
                 const uint64_t *telescopes, double sample_rate) {
    //____________________________________________________________
    // MPI information

    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (rank == 0) {
        printf("--- Gap-filling routine ---\n");
        fflush(stdout);
    }

    //____________________________________________________________
    // Data distribution

    int64_t M = 0;
    int64_t gif = 0;
    int m = data_size_proc[rank];
    int nb_blocks_tot = 0;

    for (int i = 0; i < size; i++) {
        M += data_size_proc[i];
    }

    for (int i = 0; i < rank; i++) {
        gif += data_size_proc[i];
    }

    MPI_Allreduce(&nb_blocks_loc, &nb_blocks_tot, 1, MPI_INT, MPI_SUM, comm);

    if (rank == 0) {
        printf("  total tod size  = %ld \n", M);
        printf("  local tod size  = %d \n", m);
        printf("  total intervals = %d \n", nb_blocks_tot);
        printf("  local intervals = %d \n", nb_blocks_loc);
        fflush(stdout);
    }

    //____________________________________________________________
    // Pointing matrix initialization

    Mat A;
    Gap G;

    // fake pointing weights
    auto *weights = (double *)(calloc(m * nnz, sizeof(double)));

    A.trash_pix = 0;
    MatInit(&A, m, nnz, indices, weights, 6, comm);

    // Build pixel-to-time-domain mapping

    G.ngap = build_pixel_to_time_domain_mapping(&A);

    build_gap_struct(gif, &G, &A);

    free(A.id_last_pix);
    free(A.ll);

    if (rank == 0) {
        printf("Pixel-to-time-domain mapping -> detected %d timestream gaps \n",
               G.ngap);
        fflush(stdout);
    }

    //____________________________________________________________
    // Build Toeplitz matrices

    MPI_Barrier(comm);

    Block *tpltzblocks_N;
    Block *tpltzblocks_Nm1;
    Tpltz N;
    Tpltz Nm1;

    // flags for Toeplitz product strategy
    Flag flag_stgy;
    flag_stgy_init_auto(&flag_stgy);
    flag_stgy.flag_skip_build_gappy_blocks = 1;

    // Block definition
    tpltzblocks_N = (Block *)malloc(nb_blocks_loc * sizeof(Block));
    tpltzblocks_Nm1 = (Block *)malloc(nb_blocks_loc * sizeof(Block));
    defineBlocks_avg(tpltzblocks_N, tt, nb_blocks_loc, local_blocks_sizes,
                     lambda, gif);
    defineBlocks_avg(tpltzblocks_Nm1, inv_tt, nb_blocks_loc, local_blocks_sizes,
                     lambda, gif);

    // Matrix definition
    defineTpltz_avg(&N, M, 1, 1, tpltzblocks_N, nb_blocks_loc, nb_blocks_tot,
                    gif, m, flag_stgy, comm);

    defineTpltz_avg(&Nm1, M, 1, 1, tpltzblocks_Nm1, nb_blocks_loc,
                    nb_blocks_tot, gif, m, flag_stgy, comm);

    //____________________________________________________________
    // Gap filling

    MPI_Barrier(comm);

    compute_gaps_per_block(&G, Nm1.nb_blocks_loc, Nm1.tpltzblocks);
    copy_gap_info(Nm1.nb_blocks_loc, Nm1.tpltzblocks, N.tpltzblocks);

    // set signal in all gaps to zero
    // reset_relevant_gaps(signal, &Nm1, &G);
    // here we only have noise

    // set pointing weights to zero for gaps
    condition_extra_pix_zero(&A);

    // call routine to generate constrained realization
    MPI_Barrier(comm);
    double start = MPI_Wtime();

    perform_gap_filling(comm, &N, &Nm1, noise, &G, realization, detindxs,
                        obsindxs, telescopes, sample_rate, true);

    MPI_Barrier(comm);
    double end = MPI_Wtime();

    if (rank == 0) {
        printf("Performed gap-filling in %lf s\n\n", end - start);
        fflush(stdout);
    }

    //____________________________________________________________

    MatFree(&A);
    A.indices = nullptr;
    A.values = nullptr;

    free(tpltzblocks_N);
    free(tpltzblocks_Nm1);
    free(G.id0gap);
    free(G.lgap);
}

void sim_constrained_block(bool init, bool finalize, int samples, int lambda,
                           int w0, double *tt, double *inv_tt, double *noise,
                           int *pix, uint64_t realization, uint64_t detindx,
                           uint64_t obsindx, uint64_t telescope,
                           double sample_rate) {
    // Initialize MPI
    if (init)
        MPI_Init(nullptr, nullptr);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int nb_blocks_loc = 1;
    const int nb_blocks_tot = 1;
    const int64_t gif = 0;
    const int nnz = 1;

    const int local_blocks_sizes[nb_blocks_loc] = {samples};

    //____________________________________________________________
    // Pointing matrix initialization

    Mat A;
    Gap G;

    A.trash_pix = 0;
    MatSetIndices(&A, samples, nnz, pix);
    MatLocalShape(&A, 3);

    // Build pixel-to-time-domain mapping

    G.ngap = build_pixel_to_time_domain_mapping(&A);

    build_gap_struct(gif, &G, &A);

    free(A.id_last_pix);
    free(A.ll);

    //____________________________________________________________
    // Build Toeplitz matrices

    Block tpltzblocks_N[nb_blocks_loc];
    Block tpltzblocks_Nm1[nb_blocks_loc];
    Tpltz N;
    Tpltz Nm1;

    // flags for Toeplitz product strategy
    Flag flag_stgy;
    flag_stgy_init_auto(&flag_stgy);
    flag_stgy.flag_skip_build_gappy_blocks = 1;

    // Block definition
    defineBlocks_avg(tpltzblocks_N, tt, nb_blocks_loc,
                     (void *)local_blocks_sizes, lambda, gif);
    defineBlocks_avg(tpltzblocks_Nm1, inv_tt, nb_blocks_loc,
                     (void *)local_blocks_sizes, lambda, gif);

    // Matrix definition
    defineTpltz_avg(&N, samples, 1, 1, tpltzblocks_N, nb_blocks_loc,
                    nb_blocks_tot, gif, samples, flag_stgy, MPI_COMM_WORLD);

    defineTpltz_avg(&Nm1, samples, 1, 1, tpltzblocks_Nm1, nb_blocks_loc,
                    nb_blocks_tot, gif, samples, flag_stgy, MPI_COMM_WORLD);

    //____________________________________________________________
    // Gap filling

    compute_gaps_per_block(&G, Nm1.nb_blocks_loc, Nm1.tpltzblocks);
    copy_gap_info(Nm1.nb_blocks_loc, Nm1.tpltzblocks, N.tpltzblocks);

    mappraiser::GapFillInfo gfi(1, G.ngap, 0);
    mappraiser::sim_constrained_noise_block(gfi, &N, &Nm1, noise, &G, w0,
                                            realization, detindx, obsindx,
                                            telescope, sample_rate);

    //____________________________________________________________

    // MatFree(&A);
    // A.indices = nullptr;
    // A.values  = nullptr;

    free(G.id0gap);
    free(G.lgap);

    if (finalize)
        MPI_Finalize();
}
