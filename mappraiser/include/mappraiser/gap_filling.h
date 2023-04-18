#ifndef MAPPRAISER_GAP_FILLING_H
#define MAPPRAISER_GAP_FILLING_H

#include <midapack.h>

#ifdef __cplusplus
namespace mappraiser {
    void psd_from_tt ( int fftlen, int lambda, int psdlen, const double *tt, double *psd, double rate = 200.0 );

    double compute_mean ( int samples, double *buf, bool subtract = false );

    double compute_variance ( int samples, const double &mean, double *buf );

    void sim_noise_tod ( int samples, int lambda, const double *tt, double *buf, u_int64_t realization,
                         u_int64_t detindx, u_int64_t obsindx, u_int64_t telescope, double var_goal = 1.0 );

    void sim_constrained_noise_block ( Tpltz *N_block, Tpltz *Nm1_block, double *noise, Gap *gaps,
                                       u_int64_t realization, u_int64_t detindx, u_int64_t obsindx,
                                       u_int64_t telescope );

    void sim_constrained_noise ( Tpltz *N, Tpltz *Nm1, double *noise, Gap *gaps, u_int64_t realization,
                                 const u_int64_t *detindxs, const u_int64_t *obsindxs, const u_int64_t *telescopes );
}
#endif

// Here the routine destined to be called by C code
#ifdef __cplusplus
extern "C" {
#endif

void sim_constrained_noise ( Tpltz *N, Tpltz *Nm1, double *noise, Gap *gaps, u_int64_t realization,
                             const u_int64_t *detindxs, const u_int64_t *obsindxs, const u_int64_t *telescopes );

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
                   const u_int64_t *telescopes );

#ifdef __cplusplus
}
#endif

#endif //MAPPRAISER_GAP_FILLING_H
