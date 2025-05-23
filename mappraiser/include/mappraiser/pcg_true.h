#ifndef MAPPRAISER_PCG_TRUE_H
#define MAPPRAISER_PCG_TRUE_H

#include <midapack.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <mappraiser/precond.h>
#include <mappraiser/solver_info.h>
#include <mappraiser/weight.h>

// PCG routine for maxL map-making
void PCG_maxL(const Mat *A, const Precond *M, const WeightMatrix *W, double *x,
              const double *data, SolverInfo *si);

// PCG routine for MT map-making
int PCG_GLS_templates(char *outpath, char *ref, const Mat *A, const Precond *M,
                      WeightMatrix *W, TemplateClass *X, double *B,
                      int **sweeptstamps, int npoly, int ground, int nhwp,
                      int *nsweeps, int **az_binned, int n_sss_bins,
                      int *hwp_bins, double ***hwp_mod, double delta_t,
                      int store_hwp, int nces, int *ces_length,
                      int nb_blocks_loc, double *x, double *b, double *noise,
                      double tol, int K, double sampling_freq);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_PCG_TRUE_H
