#ifndef MAPPRAISER_PRECOND_H
#define MAPPRAISER_PRECOND_H

#include <mappraiser/mapping.h>
#include <mappraiser/weight.h>
#include <midapack.h>

typedef enum { BJ = 0, APRIORI = 1, APOSTERIORI = 2 } PrecondType;

typedef struct precond_t {
    PrecondType ptype; // see above
    int n;             // total number of pixels (= n_valid + n_extra)
    int n_valid;       // number of valid pixels
    int n_extra;       // number of extra pixels
    double *pixpond;   // pixel share ponderation
    Mat BJ_inv;        // inverse preconditioner matrix
#if 0
    Mat BJ;          // preconditioner matrix
#endif
    int Zn; // size of deflation space (2-lvl preconditioners)

    /* 2 lvl only (NULL otherwise) */
    double **Z;
    double **AZ;
    double *Em1; // size Zn*Zn
    double *Qg;  // size n
    double *AQg; // size n
    double *Qtx; // size Zn
    double *w;   // size Zn
} Precond;

// Preconditioner constructor (Block Jacobi)
Precond *newPrecondBJ(Mat *A, Tpltz *Nm1, double *cond, int *lhits,
                      GapStrategy gs, Gap *Gaps, int64_t gif,
                      int *local_blocks_sizes);

// 2lvl constructor
// to be called after getting the Block Jacobi from the previous factory
// and instantiating the weighting operator W
void buildPrecond2lvl(Precond *P, Mat *A, WeightMatrix *W, double *x, double *d);

// Product of the preconditioner with a map vector
void applyPrecond(Precond *p, const Mat *A, double *g, double *Cg);

// Free memory of the preconditioner
void PrecondFree(Precond *p);

#endif // MAPPRAISER_PRECOND_H
