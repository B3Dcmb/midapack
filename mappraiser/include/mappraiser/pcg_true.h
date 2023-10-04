#ifndef MAPPRAISER_PCG_TRUE_H
#define MAPPRAISER_PCG_TRUE_H

#include "midapack.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include <stdint.h>

#include <mapping.h>

// Pixel share ponderation to deal with overlapping pixels between multiple MPI
// procs
void get_pixshare_pond(Mat *A, double *pixpond);

// PCG routine
int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz *Nm1, Tpltz *N,
                 double *x, double *b, double *noise, double *cond, int *lhits,
                 double tol, int K, int precond, int Z_2lvl, GapStrategy gs,
                 bool do_gap_filling, Gap *Gaps, int64_t gif,
                 uint64_t realization, const uint64_t *detindxs,
                 const uint64_t *obsindxs, const uint64_t *telescopes,
                 double sample_rate);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_PCG_TRUE_H
