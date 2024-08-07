#ifndef MAPPRAISER_ECG_H
#define MAPPRAISER_ECG_H

#include <midapack.h>

int ECG_GLS(char *outpath, char *ref, Mat *A, Tpltz *Nm1, Mat *BJ_inv,
            double *pixpond, double *x, double *b, double *noise, double tol,
            int maxIter, int enlFac, int ortho_alg, int bs_red);

#endif // MAPPRAISER_ECG_H
