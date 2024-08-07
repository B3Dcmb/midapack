#ifndef MAPPRAISER_H
#define MAPPRAISER_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

void MLmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond,
           int Z_2lvl, int pointing_commflag, double tol, int maxiter,
           int enlFac, int ortho_alg, int bs_red, int nside, int gap_stgy,
           bool do_gap_filling, uint64_t realization, void *data_size_proc,
           int nb_blocks_loc, void *local_blocks_sizes, double sample_rate,
           uint64_t *detindxs, uint64_t *obsindxs, uint64_t *telescopes,
           int Nnz, void *pix, void *pixweights, void *signal, double *noise,
           int lambda, double *inv_tt, double *tt);

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
           double **hwp_angle, int nces);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_H
