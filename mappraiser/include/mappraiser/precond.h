#ifndef MAPPRAISER_PRECOND_H
#define MAPPRAISER_PRECOND_H

#include <midapack.h>

typedef struct precond_t {
    int precond;     // 0 = BJ, 1 = 2lvl a priori, 2 = 2lvl a posteriori
    int n;           // total number of pixels (= n_valid + n_extra)
    int n_valid;     // number of valid pixels
    int n_extra;     // number of extra pixels
    double *pixpond; // pixel share ponderation
    Mat BJ_inv;      // inverse preconditioner matrix
    Mat BJ;          // preconditioner matrix
    int Zn;          // size of deflation space (2-lvl preconditioners)

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
int precondblockjacobilike(Mat *A, Tpltz *Nm1, Mat *BJ_inv, Mat *BJ, double *b,
                           double *noise, double *cond, int *lhits, Gap *Gaps,
                           int64_t gif);

// Preconditioner constructor
void build_precond(Precond **out_p, double **out_pixpond, Mat *A, Tpltz *Nm1,
                   double **in_out_x, double *b, double *noise, double *cond,
                   int *lhits, double tol, int Zn, int precond, Gap *Gaps,
                   int64_t gif);

// Product of the preconditioner with a map vector
void apply_precond(Precond *p, const Mat *A, Tpltz *Nm1, double *g, double *Cg);

// Free memory of the preconditioner
void free_precond(Precond **in_out_p);

#endif // MAPPRAISER_PRECOND_H
