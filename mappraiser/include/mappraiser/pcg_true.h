#ifndef MAPPRAISER_PCG_TRUE_H
#define MAPPRAISER_PCG_TRUE_H

#include "midapack.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <noise_weighting.h>
#include <precond.h>
#include <solver_info.h>

// PCG routine
void PCG_mm(Mat *A, Precond *M, Tpltz *Nm1, Tpltz *N, WeightStgy ws, Gap *G,
            double *x, const double *b, SolverInfo *si);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_PCG_TRUE_H
