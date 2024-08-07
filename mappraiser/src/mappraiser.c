/**
 * @file mappraiser.c
 * @brief High-level map-making routines of the Mappraiser library
 * @authors Simon Biquard
 * @date August 2024
 */

#include <stdlib.h>

#include <mappraiser.h>
#include <mappraiser/make_maps.h>

void MLmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond,
           int Z_2lvl, int pointing_commflag, double tol, int maxiter,
           int enlFac, int ortho_alg, int bs_red, int nside, int gap_stgy,
           bool do_gap_filling, uint64_t realization, void *data_size_proc,
           int nb_blocks_loc, void *local_blocks_sizes, double sample_rate,
           uint64_t *detindxs, uint64_t *obsindxs, uint64_t *telescopes,
           int Nnz, void *pix, void *pixweights, void *signal, double *noise,
           int lambda, double *inv_tt, double *tt) {
    // call the main map-making routine with method = 0
    make_maps(comm, 0, outpath, ref, solver, precond, Z_2lvl, pointing_commflag,
              tol, maxiter, enlFac, ortho_alg, bs_red, nside, gap_stgy,
              do_gap_filling, realization, data_size_proc, nb_blocks_loc,
              local_blocks_sizes, sample_rate, detindxs, obsindxs, telescopes,
              Nnz, pix, pixweights, signal, noise, lambda, inv_tt, tt, 0, 0, 0,
              0, 0, NULL, NULL, NULL, NULL, NULL, NULL, 0);
}

void MTmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond,
           int Z_2lvl, int pointing_commflag, double tol, int maxiter,
           int enlFac, int ortho_alg, int bs_red, int nside, int gap_stgy,
           bool do_gap_filling, uint64_t realization, void *data_size_proc,
           int nb_blocks_loc, void *local_blocks_sizes, double sample_rate,
           uint64_t *detindxs, uint64_t *obsindxs, uint64_t *telescopes,
           int Nnz, void *pix, void *pixweights, void *signal, double *noise,
           int lambda, double *inv_tt, double *tt, int npoly, int nhwp,
           double delta_t, int ground, int n_sss_bins, int **sweeptstamps,
           int *nsweeps, double **az, double *az_min, double *az_max,
           double **hwp_angle, int nces) {
    // call the main map-making routine with method = 1
    make_maps(comm, 1, outpath, ref, solver, precond, Z_2lvl, pointing_commflag,
              tol, maxiter, enlFac, ortho_alg, bs_red, nside, gap_stgy,
              do_gap_filling, realization, data_size_proc, nb_blocks_loc,
              local_blocks_sizes, sample_rate, detindxs, obsindxs, telescopes,
              Nnz, pix, pixweights, signal, noise, lambda, inv_tt, tt, npoly,
              nhwp, delta_t, ground, n_sss_bins, sweeptstamps, nsweeps, az,
              az_min, az_max, hwp_angle, nces);
}
