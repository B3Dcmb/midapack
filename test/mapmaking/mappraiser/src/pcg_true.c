// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// PCG routine applied to the mapmaking equation
// This can use the diagonal or the block-diagonal jacobi preconditionners

/** @file   pcg_true.c
 @author Frederic Dauvergne
 @date   November 2012
 @Last_update May 2019 by Hamza El Bouhargani
 @Last_update June 2020 by Aygul Jamal */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <mkl.h>
#include "midapack.h"
#include "mappraiser.h"


int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz Nm1, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int K, int precond)
{
    int    i, j, k;     // some indexes
    int    m, n;        // number of local time samples, number of local pixels
    int    rank, size;
    double localreduce; // reduce buffer
    double st, t;       // timers
    double solve_time = 0.0;
    double res, res0, res_rel;
    
    double *_g, *ACg, *Ah, *Nm1Ah; // time domain vectors
    double *g, *gp, *gt, *Cg, *h;  // map domain vectors
    double *AtNm1Ah;               // map domain
    double ro, gamma, coeff;       // scalars
    double g2pix, g2pixp, g2pix_polak;

    struct Precond *p = NULL;
    double *pixpond;

    // if we want to use the true norm to compute the residual
    int TRUE_NORM = 1;  //0: No ; 1: Yes

    FILE *fp;

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);
    m = A->m;                    

    st = MPI_Wtime();
    
    //void build_precond(struct Precond **out_p, double **out_pixpond, int *out_n, Mat *A, Tpltz *Nm1, double **in_out_x, double *b, const double *noise, double *cond, int *lhits, double tol, int Zn, int precond)
    build_precond(&p, &pixpond, &n, A, &Nm1, &x, b, noise, cond, lhits, tol, size /* Zn */, precond);

    // map domain
    h = (double *) malloc(n * sizeof(double)); // descent direction
    g = (double *) malloc(n * sizeof(double));
    gp = (double *) malloc(n * sizeof(double));
    AtNm1Ah = (double *) malloc(n * sizeof(double));
    
    // time domain
    Ah = (double *) malloc(m * sizeof(double));

    _g = Ah;
    Cg = AtNm1Ah;
    Nm1Ah = Ah;

    //printf("n=%d, m=%d, A.nnz=%d \n", n, m, A->nnz);

    t = MPI_Wtime();
    if (rank == 0) {
        printf("[rank %d] compute precond time = %lf \n", rank, t - st);
        fflush(stdout);
    }

    st = MPI_Wtime();

    MatVecProd(A, x, _g, 0);

    for (i = 0; i < m; i++)
        _g[i] = b[i] + noise[i] - _g[i];

    stbmmProd(Nm1, _g); // _g = Nm1 (Ax-b)

    TrMatVecProd(A, _g, g, 0); // g = At _g
    
    //void apply_precond(struct precond *p, Mat *A, Tpltz *Nm1, double *g, double *Cg)
    apply_precond(p, A, &Nm1, g, Cg);

    for (j = 0; j < n; j++) // h = -Cg
        h[j] = Cg[j];
    
    g2pix = 0.0; // g2 = "res"
    localreduce = 0.0;
    for (i = 0; i < n; i++) // g2 = (Cg, g)
        localreduce += Cg[i] * g[i] * pixpond[i];

    MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    t = MPI_Wtime();
    solve_time += (t - st);

    // Just to check with the true norm:
    if (TRUE_NORM == 1) {
        res = 0.0; // g2 = "res"
        localreduce = 0.0;
        for (i = 0; i < n; i++) // g2 = (Cg, g)
            localreduce += g[i] * g[i] * pixpond[i];

        MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    else {        
        res = g2pix;
    }

    double g2pixB = g2pix;
    double tol2rel = tol * tol * res; //tol*tol*g2pixB; //*g2pixB; //tol; //*tol*g2;
    res0 = res;
    // Test if already converged
    if (rank == 0) {
        printf("res = %e, res0 = %e, g2pix = %e, g2pixB = %e\n", res, res0, g2pix, g2pixB);
        res_rel = sqrt(res) / sqrt(res0);
        printf("k=%d res_g2pix=%e res_g2pix_rel=%e res_rel=%e time=%lf\n", 0, g2pix, sqrt(g2pix) / sqrt(g2pixB), sqrt(res) / sqrt(res0), t - st);
        char filename[256];
        sprintf(filename,"%s/pcg_residuals_%s.dat", outpath, ref);
        fp = fopen(filename, "wb");
        fwrite(&res_rel, sizeof(double), 1, fp);
    }

    if (res <= tol) {
        if (rank == 0)
            printf("--> converged (%e < %e)\n", res, tol);
        k = K; // to not enter inside the loop
    }

    st = MPI_Wtime();
    fflush(stdout);
    
    
    // PCG Descent Loop *********************************************
    for (k = 1; k < K ; k++){
        
        // Swap g backup pointers (Ribière-Polak needs g from previous iteration)
        gt = gp;
        gp = g;
        g = gt;
        
        MatVecProd(A, h, Ah, 0);            // Ah = A h
        stbmmProd(Nm1, Nm1Ah);              // Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)
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

	//void apply_precond(struct precond *p, Mat *A, Tpltz *Nm1, double *g, double *Cg)
	apply_precond(p, A, &Nm1, g, Cg);

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
        if (TRUE_NORM == 1) {
            localreduce = 0.0;
            for (i = 0; i < n; i++) // g2 = (Cg, g)
                localreduce += g[i] * g[i] * pixpond[i];
            
            MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
        else {
            res = g2pix_polak;
        }
        
        if (rank == 0){ //print iterate info
            printf("res = %e, res0 = %e, g2pix = %e, g2pixB = %e\n", res, res0, g2pix_polak, g2pixB);
            res_rel = sqrt(res) / sqrt(res0);
            printf("k=%d res_g2pix=%e res_g2pix_rel=%e res_rel=%e time=%lf \n", k, g2pix_polak, sqrt(g2pix_polak) / sqrt(g2pixB), sqrt(res) / sqrt(res0), t - st);
            fwrite(&res_rel, sizeof(double), 1, fp);
        }
        fflush(stdout);
                
        if (res <= tol2rel) {
            if (rank == 0) {
                printf("--> converged (%e < %e) \n", res, tol2rel);
                printf("--> i.e. \t (%e < %e) \n", sqrt(res / res0), tol);
                printf("--> solve time = %lf \n", solve_time);
                fclose(fp);
            }
            break;
        }
  
        if (g2pix_polak > g2pixp) {
            if (rank == 0)
                printf("--> g2pix > g2pixp pb (%e > %e) \n", g2pix, g2pixp);
        }

        st = MPI_Wtime();
        
        //gamma = g2pix / g2pixp;
        gamma = g2pix_polak / g2pixp;

        for (j = 0; j < n; j++) // h = h * gamma + Cg
            h[j] = h[j] * gamma + Cg[j];

    } // End loop

    if (k == K) { // check unconverged
        if (rank == 0) {
            printf("--> unconverged, max iterate reached (%le > %le)\n", g2pix, tol2rel);
            fclose(fp);
        }
    }

    if (rank == 0)
        printf("--> res_g2pix=%e\n", g2pix);

    free(h);
    free(Ah);
    free(g);
    free(AtNm1Ah);
    free_precond(&p);

    return 0;
}
