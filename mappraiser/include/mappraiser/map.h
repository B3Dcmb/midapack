//
// Created by sbiquard on 5/17/23.
//

#ifndef MAPPRAISER_MAP_H
#define MAPPRAISER_MAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <mpi.h>

void MLmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond, int Z_2lvl, int pointing_commflag,
           double tol, int maxiter, int enlFac, int ortho_alg, int bs_red, int nside, int gap_stgy,
           u_int64_t realization, void *data_size_proc, int nb_blocks_loc, void *local_blocks_sizes, double sample_rate,
           u_int64_t *detindxs, u_int64_t *obsindxs, u_int64_t *telescopes, int Nnz, void *pix, void *pixweights,
           void *signal, double *noise, int lambda, double *inv_tt, double *tt);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_MAP_H
