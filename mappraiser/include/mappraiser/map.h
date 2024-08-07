//
// Created by sbiquard on 5/17/23.
//

#ifndef MAPPRAISER_MAP_H
#define MAPPRAISER_MAP_H

#include <mpi.h>

void MLmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond, int Z_2lvl, int pointing_commflag,
           double tol, int maxiter, int enlFac, int ortho_alg, int bs_red, int nside, void *data_size_proc,
           int nb_blocks_loc, void *local_blocks_sizes, int Nnz, void *pix, void *pixweights, void *signal,
           double *noise, int lambda, double *invtt);

void MTmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond, int Z_2lvl, int pointing_commflag, int npoly, int nhwp, double delta_t,
           int ground,int n_sss_bins, double tol, int maxiter, int enlFac, int ortho_alg, int bs_red, int nside, int **sweeptstamps,
           int *nsweeps, double **az, double *az_min, double *az_max, double **hwp_angle, int nces, void *data_size_proc, int nb_blocks_loc,
           void *local_blocks_sizes, int Nnz, void *pix, void *pixweights, void *signal, double *noise, double sampling_freq, double *invtt);

#endif // MAPPRAISER_MAP_H
