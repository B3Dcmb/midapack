#ifndef MAPPRAISER_ECG_H
#define MAPPRAISER_ECG_H

#include <midapack.h>

#ifdef WITH_ECG
int ECG_GLS(char *outpath, char *ref, Mat *A, Tpltz *Nm1, double *x, double *b, double *noise, double *cond, int *lhits,
            double tol, int maxIter, int enlFac, int ortho_alg, int bs_red, Gap *Gaps, int64_t gif);
#endif

#endif //MAPPRAISER_ECG_H