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
#include "mappraiser/gap_filling.h"
#include "mappraiser/noise_weighting.h"
#include "mappraiser/precond.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz *Nm1, Tpltz *N, double *x, double *b, double *noise,
                 double *cond, int *lhits, double tol, int K, int precond, int Z_2lvl, Gap *Gaps, int64_t gif,
                 int gap_stgy, u_int64_t realization, const u_int64_t *detindxs, const u_int64_t *obsindxs,
                 const u_int64_t *telescopes) {
    int    i, j, k; // some indexes
    int    m, n;    // number of local time samples, number of local pixels
    int    rank, size;
    double localreduce; // reduce buffer
    double st, t;       // timers
    double solve_time = 0.0;
    double res, res0, res_rel;

    double *_g, *Ah, *Nm1Ah;      // time domain vectors
    double *g, *gp, *gt, *Cg, *h; // map domain vectors
    double *AtNm1Ah;              // map domain
    double  ro, gamma, coeff;     // scalars
    double  g2pix, g2pixp, g2pix_polak;

    Precond *p = NULL;
    double  *pixpond;

    // if we want to use the true norm to compute the residual
    int TRUE_NORM = 1; // 0: No ; 1: Yes

    FILE *fp;

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);
    m = A->m;

    st = MPI_Wtime();

    if (Z_2lvl == 0) Z_2lvl = size;

    build_precond(&p, &pixpond, &n, A, Nm1, &x, b, noise, cond, lhits, tol, Z_2lvl, precond, Gaps, gif);

    t = MPI_Wtime();
    if (rank == 0) {
        printf("[rank %d] Preconditioner computation time = %lf\n", rank, t - st);
        printf("[rank %d] trash_pix flag = %d\n", rank, A->trash_pix);
        printf("[rank %d] nbr sky pixels = %d\n", rank, n);
        fflush(stdout);
    }

    // TODO signal = signal + noise?

    WeightStgy strategy;
gap_treatment:
    switch (gap_stgy) {
        case 0: // gap-filling
            if (rank == 0) printf("[proc %d] gap_stgy = %d (gap-filling)\n", rank, gap_stgy);
            sim_constrained_noise(N, Nm1, noise, Gaps, realization, detindxs, obsindxs, telescopes, false);
            strategy = BASIC;
            break;
        case 1: // nested PCG for noise-weighting
            if (rank == 0) printf("[proc %d] gap_stgy = %d (nested PCG)\n", rank, gap_stgy);
            strategy = ITER;
            break;
        case 2: // gap-filling + nested PCG afterwards
            if (rank == 0) printf("[proc %d] gap_stgy = %d (gap-filling + nested PCG)\n", rank, gap_stgy);
            sim_constrained_noise(N, Nm1, noise, Gaps, realization, detindxs, obsindxs, telescopes, false);
            strategy = ITER_IGNORE;
            break;
        default:
            if (rank == 0) printf("[proc %d] invalid gap_stgy (%d), defaulting to 0\n", rank, gap_stgy);
            gap_stgy = 0;
            goto gap_treatment;
    }
    fflush(stdout);

    // map domain objects memory allocation
    h       = (double *) malloc(n * sizeof(double)); // descent direction
    g       = (double *) malloc(n * sizeof(double)); // residual
    gp      = (double *) malloc(n * sizeof(double)); // previous residual
    AtNm1Ah = (double *) malloc(n * sizeof(double));

    // time domain objects memory allocation
    Ah = (double *) malloc(m * sizeof(double));

    _g    = Ah;
    Cg    = AtNm1Ah;
    Nm1Ah = Ah;

    st = MPI_Wtime();

    // Compute RHS - initial guess
    MatVecProd(A, x, _g, 0);

    for (i = 0; i < m; i++) _g[i] = b[i] + noise[i] - _g[i];

    apply_weights(Nm1, N, Gaps, _g, strategy, false); // _g = Nm1 (d-Ax0)  (d = signal + noise)

    TrMatVecProd(A, _g, g, 0); // g = At _g

    apply_precond(p, A, Nm1, g, Cg);

    for (j = 0; j < n; j++) // h = Cg
        h[j] = Cg[j];

    g2pix       = 0.0; // g2pix = "Cg res"
    localreduce = 0.0;
    for (i = 0; i < n; i++) // g2pix = (Cg, g)
        localreduce += Cg[i] * g[i] * pixpond[i];

    MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, A->comm);

    t = MPI_Wtime();
    solve_time += (t - st);

    if (TRUE_NORM == 1) {
        res         = 0.0; // g2 = "res"
        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (g, g)
            localreduce += g[i] * g[i] * pixpond[i];

        MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    } else {
        res = g2pix;
    }

    double g2pixB  = g2pix;
    double tol2rel = tol * tol * res; // tol*tol*g2pixB; //*g2pixB; //tol; //*tol*g2;
    res0           = res;
    // Test if already converged
    if (rank == 0) {

        res_rel = sqrt(res) / sqrt(res0);
        printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", 0, res, g2pix, res_rel, t - st);
        char filename[256];
        sprintf(filename, "%s/pcg_residuals_%s.dat", outpath, ref);
        fp = fopen(filename, "wb");
        if (fp != NULL) fwrite(&res_rel, sizeof(double), 1, fp);
        fflush(stdout);
    }

    if (res <= tol) {
        if (rank == 0) printf("--> converged (%e < %e)\n", res, tol);
        k = K; // to not enter inside the loop
    }

    st = MPI_Wtime();
    fflush(stdout);

    // PCG Descent Loop *********************************************
    for (k = 1; k < K; k++) {

        // Swap g backup pointers (Ribière-Polak needs g from previous iteration)
        gt = gp;
        gp = g;
        g  = gt;

        MatVecProd(A, h, Ah, 0); // Ah = A h

        apply_weights(Nm1, N, Gaps, Nm1Ah, strategy, false); // Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)

        TrMatVecProd(A, Nm1Ah, AtNm1Ah, 0); // AtNm1Ah = At Nm1Ah

        coeff       = 0.0;
        localreduce = 0.0;
        for (i = 0; i < n; i++) localreduce += h[i] * AtNm1Ah[i] * pixpond[i];
        MPI_Allreduce(&localreduce, &coeff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        ro = g2pix / coeff;

        for (j = 0; j < n; j++) // x = x + ro * h
            x[j] = x[j] + ro * h[j];

        for (j = 0; j < n; j++)             // g = g + ro * (At Nm1 A) h
            g[j] = gp[j] - ro * AtNm1Ah[j]; // Use Ribière-Polak formula

        apply_precond(p, A, Nm1, g, Cg);

        g2pixp      = g2pix; // g2p = "res"
        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (Cg, g)
            localreduce += Cg[i] * g[i] * pixpond[i];

        MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, A->comm);

        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (Cg, g)
            localreduce += Cg[i] * (g[i] - gp[i]) * pixpond[i];

        MPI_Allreduce(&localreduce, &g2pix_polak, 1, MPI_DOUBLE, MPI_SUM, A->comm);

        t = MPI_Wtime();
        solve_time += (t - st);

        // Just to check with the true norm:
        if (TRUE_NORM == 1) {
            localreduce = 0.0;
            for (i = 0; i < n; i++) // g2 = (Cg, g)
                localreduce += g[i] * g[i] * pixpond[i];

            MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        } else {
            res = g2pix_polak;
        }

        if (rank == 0) { // print iterate info
            res_rel = sqrt(res) / sqrt(res0);
            printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", k, res, g2pix_polak, res_rel, t - st);
            if (fp != NULL) fwrite(&res_rel, sizeof(double), 1, fp);
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

        if (g2pix_polak > g2pixp) {
            if (rank == 0) printf("--> g2pix > g2pixp pb (%e > %e) \n", g2pix, g2pixp);
        }

        st = MPI_Wtime();

        // gamma = g2pix / g2pixp;
        gamma = g2pix_polak / g2pixp;

        for (j = 0; j < n; j++) // h = h * gamma + Cg
            h[j] = h[j] * gamma + Cg[j];

    } // End loop

    if (k == K) { // check unconverged
        if (rank == 0) {
            printf("--> unconverged, max iterate reached (%le > %le)\n", res, tol2rel);
            fclose(fp);
        }
    }

    if (rank == 0) printf("--> g2pix = %e\n", g2pix);

    free(h);
    free(g);
    free(gp);
    free(AtNm1Ah);
    free(Ah);
    free_precond(&p);

    return 0;
}
