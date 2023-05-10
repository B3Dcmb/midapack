#ifndef MAPPRAISER_PCG_TRUE_H
#define MAPPRAISER_PCG_TRUE_H

#include <midapack.h>
#include <stdbool.h>

// Pixel share ponderation to deal with overlapping pixels between multiple MPI
// procs
void get_pixshare_pond(Mat *A, double *pixpond);

// PCG routine
int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz *Nm1, double *x,
                 double *b, double *noise, double *cond, int *lhits, double tol,
                 int K, int precond, int Z_2lvl);

#endif // MAPPRAISER_PCG_TRUE_H
