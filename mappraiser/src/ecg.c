/**
 * @file  ecg.c
 * @authors Hamza El Bouhargani, Jan Papez
 * @date  July 2019
 * @credit  Adapted from work by Olivier Tissot
 * @description E(Enlarged) C(Conjugate) G(Gradient) solver implementation for
 * MAPPRAISER
 */

/*****************************************************************************/
/*                                  INCLUDE                                  */
/*****************************************************************************/
/* Standard headers */
#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
/* MPI */
#include <mpi.h>
/* MKL */
#include <mkl.h>
/* MAPPRAISER */
#include "mappraiser/ecg.h"
#include "mappraiser/pcg_true.h"
#include "mappraiser/precond.h"
/* preAlps */
#include "overlap_ecg.h"
/*****************************************************************************/

/*****************************************************************************/
/*                           FUNCTIONS PROTOTYPES                            */
/*****************************************************************************/
/* Build right hand side */
int get_rhs(Mat *A, Tpltz *Nm1, double *b, double *noise, double *x,
            double *rhs);
/* Preconditioner operator */
double Opmmpreconditioner(Mat *A, Mat *BJ_inv, double *X, double *Y, int ncol);
/* System matrix operator: A^T * N^-1 * A */
double Opmmmatrix(Mat *A, Tpltz *Nm1, double *X, double *Y, int ncol);
/*****************************************************************************/

/*****************************************************************************/
/*                                 CODE                                      */
/*****************************************************************************/
int ECG_GLS(char *outpath, char *ref, Mat *A, Tpltz *Nm1, double *x, double *b,
            double *noise, double *cond, int *lhits, double tol, int maxIter,
            int enlFac, int ortho_alg, int bs_red) {
    /*================ Get MPI rank & size parameters ================*/
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Force sequential execution on each MPI process
    // OT: I tested and it still works with OpenMP activated
    MKL_Set_Num_Threads(1);
    /*===================== Variables declaration ====================*/
    int i, j, k; // looping indices
    int N;       // Global number of pixels (no overlapping correction)
    int m, n;    // local number of time samples, local number of pixels x nnz
                 // (IQU)
    double *tmp; // temporary pointer
    Mat    *BJ_inv, *BJ; // Block-Jacobi preconditioner

    /*=== Pre-process degenerate pixels & build the preconditioner ===*/
    precondblockjacobilike(A, Nm1, BJ_inv, BJ, b, cond, lhits);
    // Correct the pixels counter after pre-processing
    n = A->lcount - (A->nnz) * (A->trash_pix);
    // Reallocate memory for the well-conditioned map
    tmp = realloc(x, n * sizeof(double));
    if (tmp != NULL) { x = tmp; }
    // Compute global number of pixels (w/o accounting for overlap) --useless
    N += n;
    MPI_Allreduce(MPI_IN_PLACE, &N, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0) { printf("GLOBAL N = %d \n", N); }

    /*=========== Compute the pixels share ponderation =============*/
    double *pixpond;
    pixpond = (double *) malloc(n * sizeof(double));
    get_pixshare_pond(A, pixpond);

    /*======= Build the rhs & compute initial residual norm =======*/
    double *rhs;
    rhs = (double *) malloc(n * sizeof(double));
    get_rhs(A, Nm1, b, noise, x, rhs);

    double normres_init = 0.0;
    for (i = 0; i < n; i++) { normres_init += rhs[i] * rhs[i] * pixpond[i]; }

    MPI_Allreduce(MPI_IN_PLACE, &normres_init, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    normres_init = sqrt(normres_init);
    if (rank == 0) {
        printf("GLOBAL ||r_0|| = %.6f \n", normres_init);
        fflush(stdout);
    }

    /*========================= ECG solve =========================*/
    preAlps_ECG_t ecg;
    // Set parameters
    ecg.comm        = MPI_COMM_WORLD;
    ecg.globPbSize  = N;
    ecg.locPbSize   = n;
    ecg.maxIter     = maxIter;
    ecg.enlFac      = enlFac;
    ecg.tol         = tol;
    ecg.ortho_alg   = (ortho_alg == 0 ? ORTHODIR : ORTHOMIN);
    ecg.bs_red      = (bs_red == 0 ? NO_BS_RED : ADAPT_BS);
    int rci_request = 0;
    int stop        = 0;
    // double* sol = NULL;       //solution = one vector
    // sol = (double*) malloc(n*sizeof(double));
    double *rel_res = NULL; // relative residual
    rel_res         = (double *) malloc((maxIter) * sizeof(double));
    // timings
    double time_ECG_total = 0;
    double time_AV        = 0;
    double time_invMV     = 0;

    time_ECG_total = MPI_Wtime();
    // Allocate memory and initialize variables
    preAlps_oECGInitialize(&ecg, rhs, &rci_request,
                           A->lindices + (A->nnz) * (A->trash_pix));
    ecg.normb  = normres_init;
    rel_res[0] = 1.0;

    // Finish initialization
    time_invMV += Opmmpreconditioner(A, BJ_inv, ecg.R_p, ecg.P_p, ecg.bs);
    // preconditioner ecg.R -> ecg.P
    time_AV += Opmmmatrix(A, Nm1, ecg.P_p, ecg.AP_p, ecg.bs);
    // block operator ecg.P -> ecg.AP

    // Main loop
    while (stop != 1) {
        preAlps_oECGIterate(&ecg, &rci_request, pixpond);
        if (rci_request == 0) {
            time_AV += Opmmmatrix(A, Nm1, ecg.P_p, ecg.AP_p, ecg.bs);
            // block operator ecg.P -> ecg.AP
        } else if (rci_request == 1) {

            preAlps_oECGStoppingCriterion(&ecg, &stop, pixpond);

            rel_res[ecg.iter] = ecg.res / ecg.normb;
            if (stop == 1) break;
            if (ecg.ortho_alg == ORTHOMIN) {
                time_invMV += Opmmpreconditioner(A, BJ_inv, ecg.R_p, ecg.Z_p,
                                                 ecg.enlFac);
                // preconditioner ecg.R -> ecg.Z
            } else if (ecg.ortho_alg == ORTHODIR) {
                time_invMV += Opmmpreconditioner(A, BJ_inv, ecg.AP_p, ecg.Z_p,
                                                 ecg.bs);
                // preconditioner ecg.AP -> ecg.Z
            }
        }
    }
    time_ECG_total = MPI_Wtime() - time_ECG_total;

    if (rank == 0) {
        preAlps_oECGPrint(&ecg, 5);

        printf("*** TIMING, total ECG time = %e s\n", time_ECG_total);
        printf("*** TIMING,     A x V time = %e s\n", time_AV);
        printf("*** TIMING,  invM x V time = %e s\n", time_invMV);

        char filename[256];
        sprintf(filename, "%s/rel_res_%s_enlfac_%d_o=%d.txt", outpath, ref,
                ecg.enlFac, ortho_alg);
        FILE *fp = fopen(filename, "w");
        for (i = 0; i <= ecg.iter; i++) { fprintf(fp, "%.15e\n", rel_res[i]); }
        fclose(fp);
    }
    free(rel_res);

    // Retrieve solution and free memory
    preAlps_oECGFinalize(&ecg, x);
    /*========================= Finalize ==========================*/

    // Free arrays
    if (rhs != NULL) free(rhs);
    // if (sol != NULL) free(sol);

    return 0;
}
/*****************************************************************************/

/*****************************************************************************/
/*                           FUNCTIONS DEFINITION                            */
/*****************************************************************************/
/* Build right hand side */
int get_rhs(Mat *A, Tpltz *Nm1, double *b, double *noise, double *x,
            double *rhs) {                     /* rhs = A^T*Nm1* (b - A*x_0 ) */
    int     i;                                 // some indexes
    int     m, n;
    double *_g;                                // time domain vector

    m = A->m;                                  // number of local time samples
    n = A->lcount - (A->nnz) * (A->trash_pix); // number of local pixels

    _g = (double *) malloc(m * sizeof(double));

    MatVecProd(A, x, _g, 0);
    for (i = 0; i < m; i++)              //
        _g[i] = b[i] + noise[i] - _g[i]; //
    stbmmProd(Nm1, _g);
    TrMatVecProd(A, _g, rhs, 0);

    return 0;
}

/* Preconditioner operator */
double Opmmpreconditioner(Mat *A, Mat *BJ_inv, double *X, double *Y,
                          int ncol) { // Y = M^-1 X
    double timing = MPI_Wtime();

    int     i, j;                              // some indexes
    int     n;
    double *x, *Cg;                            // map domain vector

    n = A->lcount - (A->nnz) * (A->trash_pix); // number of local pixels

    Cg = (double *) malloc(n * sizeof(double));

    for (i = 0; i < ncol; i++) {
        // get column vector x
        x = X + i * n;
        MatVecProd(BJ_inv, x, Cg, 0);
        for (j = 0; j < n; j++) { Y[i * n + j] = Cg[j]; }
    }
    free(Cg);
    return MPI_Wtime() - timing;
}

/* System matrix operator: A^T * N^-1 * A */
double Opmmmatrix(Mat *A, Tpltz *Nm1, double *X, double *Y,
                  int ncol) { /* Y = A^T*Nm1*A * X */
    double timing = MPI_Wtime();

    int     i, j;                              // some indexes
    int     m, n;
    double *_g;                                // time domain vector
    double *x, *g;                             // map domain vector

    m = A->m;                                  // number of local time samples
    n = A->lcount - (A->nnz) * (A->trash_pix); // number of local pixels

    _g = (double *) malloc(m * sizeof(double));
    g  = (double *) malloc(n * sizeof(double));

    for (i = 0; i < ncol; i++) {
        // get column vector x
        x = X + i * n;
        MatVecProd(A, x, _g, 0);
        stbmmProd(Nm1, _g);
        TrMatVecProd(A, _g, g, 0);
        for (j = 0; j < n; j++) { Y[i * n + j] = g[j]; }
    }
    free(_g);
    free(g);
    return MPI_Wtime() - timing;
}
