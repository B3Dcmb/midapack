#ifndef MAPPRAISER_PCG_TRUE_H
#define MAPPRAISER_PCG_TRUE_H

#include "midapack.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <noise_weighting.h>
#include <precond.h>

// PCG routine
int PCG_GLS_true(char *outpath, char *ref, Mat *A, Precond *P, Tpltz *Nm1,
                 Tpltz *N, double *x, const double *b, double tol, int K,
                 Gap *G, WeightStgy ws);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_PCG_TRUE_H
