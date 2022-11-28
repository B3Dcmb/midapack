// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// PCG routine applied to the map-making equation
// This can use the block-diagonal jacobi or Two-level preconditionners

/** @file   pcg_true.c
 @author Hamza El Bouhargani
 @date   May 2019
 @credit  Adapted from work by Frederic Dauvergne
 @Last_update June 2020 by Aygul Jamal */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include "mappraiser.h"

int apply_weights(Tpltz *Nm1, Tpltz *Ncov, Gap *Gaps, double *tod, int m);

int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz *Nm1, Tpltz *Ncov, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int K, int precond, int Z_2lvl, Gap *Gaps, int64_t gif)
{
    int i, j, k; // some indexes
    int m, n;    // number of local time samples, number of local pixels
    int rank, size;
    double localreduce; // reduce buffer
    double st, t;       // timers
    double solve_time = 0.0;
    double res, res0, res_rel;

    double *_g, *ACg, *Ah, *Nm1Ah; // time domain vectors
    double *g, *gp, *gt, *Cg, *h;  // map domain vectors
    double *AtNm1Ah;               // map domain
    double ro, gamma, coeff;       // scalars
    double g2pix, g2pixp, g2pix_polak;

    Precond *p = NULL;
    double *pixpond;

    // if we want to use the true norm to compute the residual
    int TRUE_NORM = 1; // 0: No ; 1: Yes

    FILE *fp;

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);
    m = A->m;

    st = MPI_Wtime();

    if (Z_2lvl == 0)
        Z_2lvl = size;

    build_precond(&p, &pixpond, &n, A, Nm1, &x, b, noise, cond, lhits, tol, Z_2lvl, precond, Gaps, gif);

    t = MPI_Wtime();
    if (rank == 0)
    {
        printf("[rank %d] Preconditioner computation time = %lf \n", rank, t - st);
        // printf("[rank %d] trash_pix flag = %d \n", rank, A->trash_pix);
        // printf("[rank %d] nbr sky pixels = %d \n", rank, n);
        fflush(stdout);
    }

    // map domain objects memory allocation
    h = (double *)malloc(n * sizeof(double));  // descent direction
    g = (double *)malloc(n * sizeof(double));  // residual
    gp = (double *)malloc(n * sizeof(double)); // residual of previous iteration
    AtNm1Ah = (double *)malloc(n * sizeof(double));

    // time domain objects memory allocation
    Ah = (double *)malloc(m * sizeof(double));

    _g = Ah;
    Cg = AtNm1Ah;
    Nm1Ah = Ah;

    st = MPI_Wtime();

    // Compute RHS - initial guess
    MatVecProd(A, x, _g, 0);

    for (i = 0; i < m; i++)
        _g[i] = b[i] + noise[i] - _g[i];

    apply_weights(Nm1, Ncov, Gaps, _g, A->m); // _g = Nm1 (d-Ax0)  (d = signal + noise)

    TrMatVecProd(A, _g, g, 0); // g = At _g

    apply_precond(p, A, Nm1, g, Cg);

    for (j = 0; j < n; j++) // h = Cg
        h[j] = Cg[j];

    g2pix = 0.0; // g2pix = "Cg res"
    localreduce = 0.0;
    for (i = 0; i < n; i++) // g2pix = (Cg, g)
        localreduce += Cg[i] * g[i] * pixpond[i];

    MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    t = MPI_Wtime();
    solve_time += (t - st);

    if (TRUE_NORM == 1)
    {
        res = 0.0; // g2 = "res"
        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (g, g)
            localreduce += g[i] * g[i] * pixpond[i];

        MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    else
    {
        res = g2pix;
    }

    double g2pixB = g2pix;
    double tol2rel = tol * tol * res; // tol*tol*g2pixB; //*g2pixB; //tol; //*tol*g2;
    res0 = res;
    // Test if already converged
    if (rank == 0)
    {

        res_rel = sqrt(res) / sqrt(res0);
        printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", 0, res, g2pix, res_rel, t - st);
        char filename[256];
        sprintf(filename, "%s/pcg_residuals_%s.dat", outpath, ref);
        fp = fopen(filename, "wb");
        fwrite(&res_rel, sizeof(double), 1, fp);
        fflush(stdout);
    }

    if (res <= tol)
    {
        if (rank == 0)
            printf("--> converged (%e < %e)\n", res, tol);
        k = K; // to not enter inside the loop
    }

    st = MPI_Wtime();
    fflush(stdout);

    // PCG Descent Loop *********************************************
    for (k = 1; k < K; k++)
    {

        // Swap g backup pointers (Ribière-Polak needs g from previous iteration)
        gt = gp;
        gp = g;
        g = gt;

        MatVecProd(A, h, Ah, 0); // Ah = A h

        apply_weights(Nm1, Ncov, Gaps, Nm1Ah, A->m); // Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)

        TrMatVecProd(A, Nm1Ah, AtNm1Ah, 0); // AtNm1Ah = At Nm1Ah

        coeff = 0.0;
        localreduce = 0.0;
        for (i = 0; i < n; i++)
            localreduce += h[i] * AtNm1Ah[i] * pixpond[i];
        MPI_Allreduce(&localreduce, &coeff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        ro = g2pix / coeff;

        for (j = 0; j < n; j++) // x = x + ro * h
            x[j] = x[j] + ro * h[j];

        for (j = 0; j < n; j++)             // g = g + ro * (At Nm1 A) h
            g[j] = gp[j] - ro * AtNm1Ah[j]; // Use Ribière-Polak formula

        apply_precond(p, A, Nm1, g, Cg);

        g2pixp = g2pix; // g2p = "res"
        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (Cg, g)
            localreduce += Cg[i] * g[i] * pixpond[i];

        MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (Cg, g)
            localreduce += Cg[i] * (g[i] - gp[i]) * pixpond[i];

        MPI_Allreduce(&localreduce, &g2pix_polak, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        t = MPI_Wtime();
        solve_time += (t - st);

        // Just to check with the true norm:
        if (TRUE_NORM == 1)
        {
            localreduce = 0.0;
            for (i = 0; i < n; i++) // g2 = (Cg, g)
                localreduce += g[i] * g[i] * pixpond[i];

            MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        else
        {
            res = g2pix_polak;
        }

        if (rank == 0)
        { // print iterate info
            res_rel = sqrt(res) / sqrt(res0);
            printf("k = %d, res = %e, g2pix = %e, res_rel = %e, time = %lf\n", k, res, g2pix_polak, res_rel, t - st);
            fwrite(&res_rel, sizeof(double), 1, fp);
        }

        fflush(stdout);

        if (res <= tol2rel)
        {
            if (rank == 0)
            {
                printf("--> converged (%e < %e) \n", res, tol2rel);
                printf("--> i.e. \t (%e < %e) \n", res_rel, tol);
                printf("--> solve time = %lf \n", solve_time);
                fclose(fp);
            }
            break;
        }

        if (g2pix_polak > g2pixp)
        {
            if (rank == 0)
                printf("--> g2pix > g2pixp pb (%e > %e) \n", g2pix, g2pixp);
        }

        st = MPI_Wtime();

        // gamma = g2pix / g2pixp;
        gamma = g2pix_polak / g2pixp;

        for (j = 0; j < n; j++) // h = h * gamma + Cg
            h[j] = h[j] * gamma + Cg[j];

    } // End loop

    if (k == K)
    { // check unconverged
        if (rank == 0)
        {
            printf("--> unconverged, max iterate reached (%le > %le)\n", res, tol2rel);
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
    free_precond(&p);

    return 0;
}

/* Weights TOD data according to the adopted noise model*/
int apply_weights(Tpltz *Nm1, Tpltz *Ncov, Gap *Gaps, double *tod, int m)
{
    int t_id; // time sample index in local data
    int i, j;

    if (Gaps->ngap == 0) /* No gaps in the timestream */
    {
        // Use straightforward loop for white noise model
        if (Nm1->tpltzblocks[0].lambda == 1)
        {
            // Here it is assumed that we use a single bandwidth for all TOD intervals, i.e. lambda is the same for all Toeplitz blocks
            t_id = 0;
            for (i = 0; i < Nm1->nb_blocks_loc; i++)
            {
                for (j = 0; j < Nm1->tpltzblocks[i].n; j++)
                {
                    tod[t_id + j] = Nm1->tpltzblocks[i].T_block[0] * tod[t_id + j];
                }
                t_id += Nm1->tpltzblocks[i].n;
            }
        }
        // Use stbmmProd routine for correlated noise model (No det-det correlations for now)
        else
            stbmmProd(Nm1, tod);
    }
    else /* Use PCG + gstbmmProd to apply the noise weights */
    {
        const int kmax = 50;     // maximum iteration count
        const double tol = 1e-6; // convergence criterion
        int k = 0;               // iteration number
        double r0, res, tol2rel; // residuals
        double coef_1, coef_2;   // scalars
        double lval_1, lval_2;   // local values before reduction

        double *x, *p = NULL;      // guess, search direction
        double *_r, *r, *z = NULL; // gradient

        x = malloc((sizeof *x) * m);
        p = malloc((sizeof *p) * m);
        _r = malloc((sizeof *_r) * m);
        r = malloc((sizeof *r) * m);
        z = malloc((sizeof *z) * m);

        if (x == NULL || p == NULL || _r == NULL || r == NULL || z == NULL)
        {
            puts("out of memory: allocation of time-domain vectors for tod weighting failed");
            exit(EXIT_FAILURE);
        }

        // initial guess
        for (i = 0; i < m; ++i)
            _r[i] = tod[i];

        // apply system matrix (_r = Nx0)
        gstbmmProd(Ncov, _r, Gaps);

        // compute initial residual (r = b - Nx0)
        for (i = 0; i < m; ++i)
        {
            r[i] = tod[i] - _r[i];
            z[i] = r[i];
        }

        // apply preconditioner (z0 = M^{-1} * r0)
        gstbmmProd(Nm1, z, Gaps);

        // set initial search direction (p0 = z0)
        for (i = 0; i < m; ++i)
            p[i] = z[i];

        // check convergence
        lval_1 = 0.0;
        for (i = 0; i < n; i++)
            lval_1 += r[i] * r[i];

        res = 0.0;
        MPI_Allreduce(&lval_1, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        r0 = res;
        tol2rel = tol * tol * res;

        if (rank == 0) /* print information on screen */
        {
            printf("[noise weighting PCG] k = %d, r0 = %e\n", k, r0);
        }

        if (res < tol * tol)
        {
            // print info on screen
            if (rank == 0)
                printf("  -> converged (%e < %e)\n", res, tol * tol);

            // do not enter the PCG loop
            k = kmax;
        }

        for (k = 1; k < kmax; ++k)
        {
            fflush(stdout);

            // apply system matrix (_r = Np)
            for (i = 0; i < m; ++i)
                _r[i] = p[i];

            gstbmmProd(Ncov, _r, Gaps);

            // compute alpha = (r,z)/(p,Np)
            lval_1 = 0.0;
            lval_2 = 0.0;
            for (i = 0; i < m; ++i)
            {
                lval_1 += r[i] * z[i];
                lval_2 += p[i] * _r[i];
            }

            coef_1 = 0.0; // (r,z)
            coef_2 = 0.0; // (p,Np)
            MPI_Allreduce(&lval_1, &coef_1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&lval_2, &coef_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // update current vector (x = x + alpha * p)
            // update residual (r = r - alpha * Np)
            for (i = 0; i < m; ++i)
            {
                x[i] = x[i] + (coef_1 / coef_2) * p[i];
                r[i] = r[i] - (coef_1 / coef_2) * _r[i];
                z[i] = r[i];
            }

            // apply preconditioner (z = M^{-1} r)
            gstbmm(Nm1, z, Prod);

            // compute coeff for new search direction
            lval_2 = 0.0;
            for (i = 0; i < m; ++i)
                lval_2 += r[i] * z[i];

            coef_2 = 0.0; // new (r,z)
            MPI_Allreduce(&lval_2, &coef_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // update search direction
            for (i = 0; i < m; ++i)
                p[i] = z[i] + (coef_2 / coef_1) * p[i];

            // check convergence
            lval_1 = 0.0;
            for (i = 0; i < n; i++)
                lval_1 += r[i] * r[i];

            res = 0.0;
            MPI_Allreduce(&lval_1, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            if (rank == 0)
            {
                printf("[k = %d] res = %e, relative = %e\n", k, res, sqrt(res / r0));
            }

            if (res < tol2rel)
            {
                if (rank == 0)
                {
                    printf("  -> converged (%e < %e) \n", res, tol2rel);
                    printf("  -> i.e.      (%e < %e) \n", sqrt(res / r0), tol);
                    // printf("  -> solve time = %lf \n", solve_time);
                }
                break;
            }
        }
        fflush(stdout);
    }

    return 0;
}
