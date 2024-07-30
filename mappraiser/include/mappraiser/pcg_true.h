#ifndef MAPPRAISER_PCG_TRUE_H
#define MAPPRAISER_PCG_TRUE_H

#include <midapack.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <mappraiser/noise_weighting.h>
#include <mappraiser/precond.h>
#include <mappraiser/solver_info.h>

// PCG routine
void PCG_mm(Mat *A, Precond *M, Tpltz *Nm1, Tpltz *N, WeightStgy ws, Gap *G,
            double *x, const double *b, SolverInfo *si);

// without preconditioning
void CG_mm(Mat *A, double *pixpond, Tpltz *Nm1, Tpltz *N, WeightStgy ws, Gap *G,
           double *x, const double *d, SolverInfo *si);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_PCG_TRUE_H
