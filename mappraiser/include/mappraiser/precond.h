#ifndef MAPPRAISER_PRECOND_H
#define MAPPRAISER_PRECOND_H

#include <mapping.h>
#include <midapack.h>

typedef struct precond_t {
    int precond;     // 0 = BJ, 1 = 2lvl a priori, 2 = 2lvl a posteriori
    int n;           // total number of pixels (= n_valid + n_extra)
    int n_valid;     // number of valid pixels
    int n_extra;     // number of extra pixels
    double *pixpond; // pixel share ponderation
    Mat BJ_inv;      // inverse preconditioner matrix
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

// Block-Jacobi preconditioner
int precondblockjacobilike(Mat *A, Tpltz *Nm1, double *vpixBlock,
                           double *vpixBlock_inv, double *cond, int *lhits);

// Preconditioner constructor
void build_precond(Precond **out_p, double **out_pixpond, Mat *A, Tpltz *Nm1,
                   double **in_out_x, double *b, double *noise, double *cond,
                   int *lhits, double tol, int Zn, int precond, GapStrategy gs,
                   Gap *Gaps, int64_t gif, int *local_blocks_sizes);

void build_BJinv(Mat *A, Tpltz *Nm1, Mat *BJ_inv, double *cond, int *lhits,
                 GapStrategy gs, Gap *Gaps, int64_t gif,
                 int *local_blocks_sizes);

// Product of the preconditioner with a map vector
void apply_precond(Precond *p, const Mat *A, double *g, double *Cg);

// Free memory of the preconditioner
void free_precond(Precond **in_out_p);

// Pixel share ponderation to deal with overlapping pixels between multiple MPI
// procs
void get_pixshare_pond(Mat *A, double *pixpond);

#endif // MAPPRAISER_PRECOND_H
