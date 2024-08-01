#ifndef MAPPRAISER_PCG_TRUE_H
#define MAPPRAISER_PCG_TRUE_H

#include <midapack.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <mappraiser/precond.h>
#include <mappraiser/solver_info.h>
#include <mappraiser/weight.h>

// PCG routine
void PCG_maxL(const Mat *A, const Precond *M, const WeightMatrix *W, double *x,
              const double *data, SolverInfo *si);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_PCG_TRUE_H
