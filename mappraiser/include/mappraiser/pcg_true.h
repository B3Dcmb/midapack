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
void PCG_mm(Mat *A, Precond *M, WeightMatrix *W, double *x, const double *data,
            SolverInfo *si);

// without preconditioning
void CG_mm(Mat *A, const double *pixpond, WeightMatrix *W, double *x,
           const double *data, SolverInfo *si);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_PCG_TRUE_H
