/**
 * @file pcg_true.c
 * @brief Preconditioned Conjugate Gradient algorithm applied to the map
 * making equation. Can use the block-diagonal Jacobi or Two-level
 * preconditioners.
 * @authors Hamza El Bouhargani (adapted from Frederic Dauvergne), Aygul Jamal
 * @date May 2019
 * @update Mar 2023 by Simon Biquard
 */

#include "mappraiser/pcg_true.h"

#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

double scalar_prod_reduce(MPI_Comm comm, Precond *P, const double *p1,
                          const double *p2, const double *p_sub);

int PCG_GLS_true(char *outpath, char *ref, Mat *A, Precond *P, Tpltz *Nm1,
                 Tpltz *N, double *x, const double *b, double tol, int K,
                 Gap *G, WeightStgy ws) {
    int i, j, k; // some indexes
    int m;       // number of local time samples
    int rank, size;
    double st, t; // timers
    double solve_time = 0.0;
    double res, res0, res_rel;

    double *_g, *Ah, *Nm1Ah;      // time domain vectors
    double *g, *gp, *gt, *Cg, *h; // map domain vectors
    double *AtNm1Ah;              // map domain
    double ro, gamma, coeff;      // scalars
    double g2pix, g2pix_prev, g2pix_polak;

    // if we want to use the true norm to compute the residual
    int TRUE_NORM = 1; // 0: No ; 1: Yes

    FILE *fp;

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);
    m = A->m;

    MPI_Barrier(A->comm);
    if (rank == 0) {
        printf("\n[rank %d] Start main PCG iterations\n", rank);
        fflush(stdout);
    }

    // map domain objects memory allocation
    h = (double *)malloc(P->n * sizeof(double));  // descent direction
    g = (double *)malloc(P->n * sizeof(double));  // residual
    gp = (double *)malloc(P->n * sizeof(double)); // previous residual
    AtNm1Ah = (double *)malloc(P->n * sizeof(double));

    // time domain objects memory allocation
    Ah = (double *)malloc(m * sizeof(double));

    _g = Ah;
    Cg = AtNm1Ah;
    Nm1Ah = Ah;

    st = MPI_Wtime();

    // Compute RHS - initial guess
    MatVecProd(A, x, _g, 0);

    for (i = 0; i < m; i++)
        _g[i] = b[i] - _g[i];

    apply_weights(Nm1, N, G, _g, ws,
                  false); // _g = Nm1 (d-Ax0)  (d = signal + noise)

    TrMatVecProd(A, _g, g, 0); // g = At _g

    apply_precond(P, A, Nm1, g, Cg);

    for (j = 0; j < P->n; j++) // h = Cg
        h[j] = Cg[j];

    // g2pix = (Cg, g)
    g2pix = scalar_prod_reduce(A->comm, P, Cg, g, NULL);

    t = MPI_Wtime();
    solve_time += (t - st);

    if (TRUE_NORM == 1) {
        // res = (g, g)
        res = scalar_prod_reduce(A->comm, P, g, g, NULL);
    } else {
        res = g2pix;
    }

    // double g2pixB = g2pix;
    double tol2rel =
        tol * tol * res; // tol*tol*g2pixB; //*g2pixB; //tol; //*tol*g2;
    res0 = res;
    // Test if already converged
    if (rank == 0) {

        res_rel = sqrt(res) / sqrt(res0);
        printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", 0,
               res, g2pix, res_rel, t - st);
        char filename[FILENAME_MAX];
        sprintf(filename, "%s/pcg_residuals_%s.dat", outpath, ref);
        fp = fopen(filename, "wb");
        if (fp != NULL) {
            fwrite(&res_rel, sizeof(double), 1, fp);
        } else {
            fputs("error when opening pcg_residuals file", stderr);
            exit(EXIT_FAILURE);
        }
    }
    fflush(stdout);

    if (res <= tol) {
        if (rank == 0)
            printf("--> converged (%e < %e)\n", res, tol);
        k = K; // to not enter inside the loop
    }

    st = MPI_Wtime();
    fflush(stdout);

    // PCG Descent Loop *********************************************
    for (k = 1; k < K; k++) {

        // Swap g backup pointers (Ribière-Polak needs g from previous
        // iteration)
        gt = gp;
        gp = g;
        g = gt;

        MatVecProd(A, h, Ah, 0); // Ah = A h

        apply_weights(Nm1, N, G, Nm1Ah, ws,
                      false); // Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)

        TrMatVecProd(A, Nm1Ah, AtNm1Ah, 0); // AtNm1Ah = At Nm1Ah

        coeff = scalar_prod_reduce(A->comm, P, h, AtNm1Ah, NULL);

        ro = g2pix / coeff;

        for (j = 0; j < P->n; j++) // x = x + ro * h
            x[j] = x[j] + ro * h[j];

        for (j = 0; j < P->n; j++)          // g = g + ro * (At Nm1 A) h
            g[j] = gp[j] - ro * AtNm1Ah[j]; // Use Ribière-Polak formula

        apply_precond(P, A, Nm1, g, Cg);

        g2pix_prev = g2pix; // g2p = "res"
        g2pix = scalar_prod_reduce(A->comm, P, Cg, g, NULL);

        g2pix_polak = scalar_prod_reduce(A->comm, P, Cg, g, gp);

        t = MPI_Wtime();
        solve_time += (t - st);

        // Just to check with the true norm:
        if (TRUE_NORM == 1) {
            res = scalar_prod_reduce(A->comm, P, g, g, NULL);
        } else {
            res = g2pix_polak;
        }

        if (rank == 0) { // print iterate info
            res_rel = sqrt(res) / sqrt(res0);
            printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n",
                   k, res, g2pix_polak, res_rel, t - st);
            if (fp != NULL)
                fwrite(&res_rel, sizeof(double), 1, fp);
        }

        fflush(stdout);

        if (res <= tol2rel) {
            if (rank == 0) {
                printf("--> converged (%e < %e) \n", res, tol2rel);
                printf("--> i.e.      (%e < %e) \n", res_rel, tol);
                printf("--> solve time = %lf \n", solve_time);
                fclose(fp);
            }
            break;
        }

        if (g2pix_polak > g2pix_prev) {
            if (rank == 0)
                printf("--> g2pix > g2pixp pb (%e > %e) \n", g2pix, g2pix_prev);
        }

        st = MPI_Wtime();

        // gamma = g2pix / g2pixp;
        gamma = g2pix_polak / g2pix_prev;

        for (j = 0; j < P->n; j++) // h = h * gamma + Cg
            h[j] = h[j] * gamma + Cg[j];

    } // End loop

    if (k == K) { // check unconverged
        if (rank == 0) {
            printf("--> unconverged, max iterate reached (%le > %le)\n", res,
                   tol2rel);
            fclose(fp);
        }
    }

    if (rank == 0)
        printf("--> g2pix = %e\n", g2pix);

    free(h);
    free(g);
    free(gp);
    free(AtNm1Ah);
    free(Ah);

    return 0;
}

double scalar_prod_reduce(MPI_Comm comm, Precond *P, const double *p1,
                          const double *p2, const double *p_sub) {
    double result;

    // contributions from valid pixels
    double local_sum = 0.0;
    for (int i = P->n_extra; i < P->n_valid; i++) {
        if (p_sub == NULL) {
            local_sum += p1[i] * p2[i] * P->pixpond[i - P->n_extra];
        } else {
            local_sum +=
                p1[i] * (p2[i] - p_sub[i]) * P->pixpond[i - P->n_extra];
        }
    }

    // reduce the results for valid pixels
    MPI_Allreduce(&local_sum, &result, 1, MPI_DOUBLE, MPI_SUM, comm);

    // contributions from extra pixels
    for (int i = 0; i < P->n_extra; i++) {
        if (p_sub == NULL) {
            result += p1[i] * p2[i];
        } else {
            result = p1[i] * (p2[i] - p_sub[i]);
        }
    }

    return result;
}
