/**
 * @file gap_filling.cpp
 * @brief Implementation of gap_filling routines
 * @author Simon Biquard
 * @date March 2023
 */

#include <iostream>
#include <fftw3.h>
#include <algorithm>
#include <cmath>
#include <complex>

#include "mappraiser/gap_filling.h"
#include "mappraiser/rng.h"
#include "mappraiser/noise_weighting.h"
#include "mappraiser/mapping.h"
#include "mappraiser/create_toeplitz.h"

void mappraiser::psd_from_tt ( int fftlen,
                               int lambda,
                               int psdlen,
                               const double *tt,
                               double *psd,
                               double rate ) {
    // Normalization
    double norm = rate * static_cast <double> (psdlen - 1);

    // FFTW variables
    double       *circ_t;
    fftw_complex *cplx_psd;
    fftw_plan    p;

    // Allocation
    circ_t   = fftw_alloc_real ( fftlen );
    cplx_psd = fftw_alloc_complex ( psdlen );

    // FFTW plan
    p = fftw_plan_dft_r2c_1d ( fftlen, circ_t, cplx_psd, FFTW_ESTIMATE);

    // Initialization of input
    std::copy ( tt, (tt + lambda), circ_t );
    std::fill ((circ_t + lambda), (circ_t + fftlen - lambda + 1), 0 );
    std::reverse_copy ( tt + 1, (tt + lambda), (circ_t + fftlen - lambda + 1));

    // Execute FFT plan
    fftw_execute ( p );

    // Compute the PSD values
    for (int i = 0; i < psdlen; ++i) {
        psd[i] = norm * std::abs ( std::complex<double> ( cplx_psd[i][0], cplx_psd[i][1] ));
    }

    // Zero out DC value
    psd[0] = 0;

    // Free allocated memory
    fftw_free ( circ_t );
    fftw_free ( cplx_psd );
    fftw_destroy_plan ( p );
}

double mappraiser::compute_mean ( int samples, double *buf, bool subtract = false ) {
    // Compute the DC level
    double DC = 0;

    for (int i = 0; i < samples; ++i) {
        DC += buf[i];
    }
    DC /= static_cast<double> (samples);

    // Remove it if needed
    if (subtract) {
        for (int i = 0; i < samples; ++i) {
            buf[i] -= DC;
        }
    }
    return DC;
}

double mappraiser::compute_mean_good ( int samples, double *buf, const bool *valid, bool subtract = false ) {
    // Compute the DC level
    double DC   = 0;
    int    good = 0;

    for (int i = 0; i < samples; ++i) {
        if (valid[i]) {
            DC += buf[i];
            ++good;
        }
    }

    DC /= static_cast<double> (good);

    // Remove it if needed
    if (subtract) {
        for (int i = 0; i < samples; ++i) {
            if (valid[i])
                buf[i] -= DC;
        }
    }
    return DC;
}

double mappraiser::compute_variance ( int samples, double mean, const double *buf ) {
    double var = 0;

    for (int i = 0; i < samples; ++i) {
        var += (buf[i] - mean) * (buf[i] - mean);
    }

    var /= static_cast<double> (samples - 1);
    return var;
}

double mappraiser::compute_variance_good ( int samples, double mean, const double *buf, const bool *valid ) {
    double var  = 0;
    int    good = 0;

    for (int i = 0; i < samples; ++i) {
        if (valid[i]) {
            var += (buf[i] - mean) * (buf[i] - mean);
            ++good;
        }
    }

    var /= static_cast<double> (good - 1);
    return var;
}

int mappraiser::find_valid_samples ( Gap *gaps, int samples, int64_t id0, bool *valid ) {
    int gap    = -1;
    int n_good = 0;

    for (int i = 0; i < samples; ++i) {
        // find closest gap before (or at) current sample
        while (gaps->id0gap[gap + 1] <= id0 + i) { ++gap; }

        if (gap == -1)
            // the sample is valid because all gaps are strictly after
            valid[i] = true;
        else {
            // the sample is valid if the gap is too short
            valid[i] = gaps->id0gap[gap] + gaps->lgap[gap] <= id0 + i;
        }

        if (valid[i])
            ++n_good;
    }
    return n_good;
}

void mappraiser::remove_baseline ( double *buf, int samples, double *baseline, const bool *valid, bool rm = true ) {
    const int w = 200; // size of window
    int       c = 0;

    for (int i = 0; i < samples; ++i) {
        double avg   = 0;
        int    count = 0;
        if (valid[i]) {
            // compute the moving average
            for (int j = std::max ( 0, i - w ); j < std::min ( samples, i + w + 1 ); ++j) {
                if (valid[j]) {
                    avg += buf[j];
                    ++count;
                }
            }
            baseline[c] = avg / static_cast<double>(count);
            ++c;
        }
    }

    // remove the baseline
    if (rm) {
        c = 0;
        for (int i = 0; i < samples; ++i) {
            if (valid[i]) {
                buf[i] -= baseline[c];
                ++c;
            }
        }
    }
}

void mappraiser::sim_noise_tod ( int samples,
                                 int lambda,
                                 const double *tt,
                                 double *buf,
                                 u_int64_t realization,
                                 u_int64_t detindx,
                                 u_int64_t obsindx,
                                 u_int64_t telescope,
                                 double var_goal,
                                 bool verbose ) {
    // Logical size of the fft
    // this could be modified to be a power of 2, for example
    int fftlen = samples;

    double    *pdata;
    fftw_plan p;

    // Allocate the input/output buffer
    pdata = fftw_alloc_real ( fftlen );

    // Create a plan for in-place half-complex -> real (HC2R) transform
    p = fftw_plan_r2r_1d ( fftlen, pdata, pdata, FFTW_HC2R, FFTW_ESTIMATE);

    // Generate Re/Im gaussian randoms in a half-complex array
    // (see https://fftw.org/fftw3_doc/The-Halfcomplex_002dformat-DFT.html)
    u_int64_t key1     = realization * 4294967296 + telescope * 65536;
    u_int64_t key2     = obsindx * 4294967296 + detindx;
    u_int64_t counter1 = static_cast<u_int64_t>(lambda) * 4294967296;
    u_int64_t counter2 = 0; // incremented in loop during generation
    mappraiser::rng_dist_normal ( samples, key1, key2, counter1, counter2, pdata );

    // Compute PSD values from the auto-correlation function
    int  psdlen = (fftlen / 2) + 1;
    auto *psd   = static_cast<double *>(std::malloc ( sizeof ( double ) * psdlen ));
    mappraiser::psd_from_tt ( fftlen, lambda, psdlen, tt, psd );

    // Multiply by the PSD
    pdata[0] *= std::sqrt ( psd[0] );

    for (int i = 1; i < (fftlen / 2); ++i) {
        double psdval = std::sqrt ( psd[i] );
        pdata[i] *= psdval;
        pdata[fftlen - i] *= psdval;
    }
    pdata[fftlen / 2] *= std::sqrt ( psd[psdlen - 1] );

    // Execute the FFT plan
    fftw_execute ( p );

    // Backward FFT: 1/N factor not included by FFTW
    for (int i = 0; i < fftlen; ++i) {
        pdata[i] /= fftlen;
    }

    // Copy as many samples as we need into our noise vector
    std::copy ( pdata, (pdata + samples), buf );

    // Zero out DC level
    double DC_subtracted = mappraiser::compute_mean ( samples, buf, true );

    if (verbose)
        std::cout << "  subtracted DC level = " << DC_subtracted << std::endl;

    // Normalize according to desired variance for the TOD
    // In theory this should not be necessary
    if (var_goal >= 0.0) {
        double var     = mappraiser::compute_variance ( samples, 0.0, buf );
        double rescale = std::sqrt ( var_goal / var );

        if (verbose)
            std::cout << "  measured sigma = " << std::sqrt ( var ) << std::endl
                      << "   (goal:" << std::sqrt ( var_goal ) << ")" << std::endl
                      << "  rescale factor = " << rescale << std::endl;

        for (int i = 0; i < samples; ++i)
            buf[i] *= rescale;
    }

    // Free allocated memory
    fftw_free ( pdata );
    fftw_destroy_plan ( p );
    std::free ( psd );
}

void mappraiser::sim_constrained_noise_block ( Tpltz *N_block,
                                               Tpltz *Nm1_block,
                                               double *noise,
                                               Gap *gaps,
                                               u_int64_t realization,
                                               u_int64_t detindx,
                                               u_int64_t obsindx,
                                               u_int64_t telescope,
                                               bool verbose ) {
    // get the number of samples, the global first index and the bandwidth
    const int     samples = N_block->tpltzblocks[0].n;
    const int64_t id0     = N_block->tpltzblocks[0].idv;
    const int     lambda  = N_block->tpltzblocks[0].lambda;

    // copy the original noise vector
    auto *rhs = static_cast<double *>(std::malloc ( sizeof ( double ) * samples ));
    std::copy ( noise, (noise + samples), rhs );

    // locate the valid samples
    auto *valid = static_cast<bool *>(std::malloc ( sizeof ( bool ) * samples ));
    int  n_good = find_valid_samples ( gaps, samples, id0, valid );

    if (verbose)
        std::cout << "  proportion of valid samples = " << 100 * n_good / samples << " %" << std::endl;

    // remove baseline (moving average)
    auto *baseline = static_cast<double *>(std::malloc ( sizeof ( double ) * n_good ));
    mappraiser::remove_baseline ( rhs, samples, baseline, valid, true );

    if (verbose)
        std::cout << "  removed baseline" << std::endl;

/*
    // compute the mean (from valid samples)
    double mean = mappraiser::compute_mean_good ( samples, rhs, valid, false );

    // compute the variance (from valid samples)
    double var = mappraiser::compute_variance_good ( samples, mean, rhs, valid );

    if (verbose) {
        std::cout << "  mean (valid)  = " << mean << std::endl
                  << "  std deviation = " << std::sqrt ( var ) << std::endl;
    }
*/

    // generate random noise realization "xi" with correlations
    auto *xi = static_cast<double *>(std::malloc ( sizeof ( double ) * samples ));
    mappraiser::sim_noise_tod ( samples, lambda, N_block->tpltzblocks[0].T_block, xi, realization, detindx,
                                obsindx, telescope, -1.0, false );

    if (verbose)
        std::cout << "  generated random noise realization" << std::endl;

    // rhs = noise - xi
    for (int i = 0; i < samples; ++i) {
        rhs[i] -= xi[i];
    }

    // invert the system N x = (noise - xi)
    if (verbose)
        std::cout << "  inverting Toeplitz system with PCG..." << std::endl;

    apply_weights ( Nm1_block, N_block, gaps, rhs, ITER, false );

    // compute the constrained realization
    auto *constrained = static_cast<double *>(std::malloc ( sizeof ( double ) * samples ));
    std::copy ( rhs, (rhs + samples), constrained );
    stbmmProd ( N_block, constrained );

    if (verbose)
        std::cout << "  computed full constrained realization" << std::endl;

    int c = 0; // to go through baseline vector

    for (int i = 0; i < samples; ++i) {
        if (valid[i]) {
            constrained[i] += xi[i] + baseline[c];
            ++c;
        } else
            constrained[i] += xi[i]; // + (baseline[c] + baseline[c + 1]) / 2;
    }

    // TODO check quality of constrained realization against original vector?
    // ...

    // copy final result into noise vector
    std::copy ( constrained, (constrained + samples), noise );

    // Free memory
    std::free ( xi );
    std::free ( rhs );
    std::free ( constrained );
    std::free ( baseline );
}

void mappraiser::sim_constrained_noise ( Tpltz *N,
                                         Tpltz *Nm1,
                                         double *noise,
                                         Gap *gaps,
                                         u_int64_t realization,
                                         const u_int64_t *detindxs,
                                         const u_int64_t *obsindxs,
                                         const u_int64_t *telescopes,
                                         bool verbose ) {
    // Loop through toeplitz blocks
    int    t_id = 0;
    double *noise_block;
    Tpltz  N_block, Nm1_block;

    for (int i = 0; i < N->nb_blocks_loc; ++i) {
        if (verbose)
            std::cout << "Processing block n° " << i << std::endl;

        // define single-block Tpltz structures
        set_tpltz_struct ( &N_block, N, &N->tpltzblocks[i] );
        set_tpltz_struct ( &Nm1_block, Nm1, &Nm1->tpltzblocks[i] );

        // pointer to current block in the tod
        noise_block = (noise + t_id);

        // compute a constrained noise realization for the current block
        mappraiser::sim_constrained_noise_block ( &N_block, &Nm1_block, noise_block, gaps, realization, detindxs[i],
                                                  obsindxs[i], telescopes[i], verbose );

        t_id += N->tpltzblocks[i].n;
    }
}

void sim_constrained_noise ( Tpltz *N,
                             Tpltz *Nm1,
                             double *noise,
                             Gap *gaps,
                             u_int64_t realization,
                             const u_int64_t *detindxs,
                             const u_int64_t *obsindxs,
                             const u_int64_t *telescopes,
                             bool verbose ) {
    mappraiser::sim_constrained_noise ( N, Nm1, noise, gaps, realization, detindxs, obsindxs, telescopes, verbose );
}

void gap_filling ( MPI_Comm comm,
                   const int *data_size_proc,
                   int nb_blocks_loc,
                   int *local_blocks_sizes,
                   int nnz,
                   double *tt,
                   double *inv_tt,
                   int lambda,
                   double *noise,
                   int *indices,
                   u_int64_t realization,
                   const u_int64_t *detindxs,
                   const u_int64_t *obsindxs,
                   const u_int64_t *telescopes ) {
    // MPI information

    int rank, size;

    MPI_Comm_rank ( comm, &rank );
    MPI_Comm_size ( comm, &size );
    if (rank == 0) {
        printf ( "--- Gap-filling routine ---\n" );
        fflush ( stdout );
    }

    // Data distribution

    int64_t M             = 0;
    int64_t gif           = 0;
    int     m             = data_size_proc[rank];
    int     nb_blocks_tot = 0;

    for (int i = 0; i < size; i++) {
        M += data_size_proc[i];
    }

    for (int i = 0; i < rank; i++) {
        gif += data_size_proc[i];
    }

    MPI_Allreduce ( &nb_blocks_loc, &nb_blocks_tot, 1, MPI_INT, MPI_SUM, comm );

    if (rank == 0) {
        printf ( "  total tod size  = %ld \n", M );
        printf ( "  local tod size  = %d \n", m );
        printf ( "  total intervals = %d \n", nb_blocks_tot );
        printf ( "  local intervals = %d \n", nb_blocks_loc );
        fflush ( stdout );
    }

    // Pointing matrix initialization

    Mat A;
    Gap G;

    // fake pointing weights
    double *weights = (double *) (calloc ( m, sizeof *weights ));

    A.trash_pix = 0;
    MatInit ( &A, m, nnz, indices, weights, 6, comm );

    // Build pixel-to-time-domain mapping

    G.ngap = build_pixel_to_time_domain_mapping ( &A );

    build_gap_struct ( gif, &G, &A );

    if (rank == 0) {
        printf ( "Pixel-to-time-domain mapping -> detected %d timestream gaps \n", G.ngap );
        fflush ( stdout );
    }

    MatFree ( &A );
    A.indices = NULL;
    A.values  = NULL;

    // Build Toeplitz matrices
    MPI_Barrier ( comm );

    Block *tpltzblocks_N;
    Block *tpltzblocks_Nm1;
    Tpltz N;
    Tpltz Nm1;

    // flags for Toeplitz product strategy
    Flag flag_stgy;
    flag_stgy_init_auto ( &flag_stgy );
    flag_stgy.flag_skip_build_gappy_blocks = 1;

    // Block definition
    tpltzblocks_N   = (Block *) malloc ( nb_blocks_loc * sizeof ( Block ));
    tpltzblocks_Nm1 = (Block *) malloc ( nb_blocks_loc * sizeof ( Block ));
    defineBlocks_avg ( tpltzblocks_N, tt, nb_blocks_loc, local_blocks_sizes,
                       lambda, gif );
    defineBlocks_avg ( tpltzblocks_Nm1, inv_tt, nb_blocks_loc, local_blocks_sizes,
                       lambda, gif );

    // Matrix definition
    defineTpltz_avg ( &N, M, 1, 1, tpltzblocks_N, nb_blocks_loc, nb_blocks_tot,
                      gif, m, flag_stgy, comm );

    defineTpltz_avg ( &Nm1, M, 1, 1, tpltzblocks_Nm1, nb_blocks_loc, nb_blocks_tot,
                      gif, m, flag_stgy, comm );

    // call routine to generate constrained realization
    MPI_Barrier ( comm );
    double start = MPI_Wtime ();

    sim_constrained_noise ( &N, &Nm1, noise, &G, realization, detindxs, obsindxs, telescopes, rank == 0 );

    MPI_Barrier ( comm );
    double end = MPI_Wtime ();

    if (rank == 0) {
        printf ( "Performed gap-filling in %lf seconds\n\n", end - start );
        fflush ( stdout );
    }

    free ( tpltzblocks_N );
    free ( tpltzblocks_Nm1 );
    free ( G.id0gap );
    free ( G.lgap );
}
