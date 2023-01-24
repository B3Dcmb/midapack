#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include "mappraiser.h"

#define PCG_TOEPLITZ_KMAX 100
#define PCG_TOEPLITZ_TOL 1e-6

void reset_tod_gaps(double *tod, Tpltz *N, Gap *Gaps);
void set_tpltz_struct(Tpltz *single_block_struct, Tpltz *full_struct, Block *block);
void PCG_single_block(Tpltz *N_block, Tpltz *Nm1_block, Gap *Gaps, double *tod_block, int verbose);
void PCG_global(Tpltz *N, Tpltz *Nm1, Gap *Gaps, double *tod, int m, int rank, int verbose);

/**
 * @brief Weight TOD with the adopted noise model.
 * An PCG is used unless the noise model has lambda=1 (white noise).
 *
 * @param Nm1 Approximate Toeplitz inverse noise covariance
 * @param N Toeplitz noise covariance
 * @param Gaps timestream gaps
 * @param tod data vector
 */
void apply_weights(Tpltz *Nm1, Tpltz *N, Gap *Gaps, double *tod)
{
    int t_id = 0; // time sample index in local data
    int i, j;

    int rank, size;
    MPI_Comm_rank(N->comm, &rank);
    MPI_Comm_size(N->comm, &size);

    if (Nm1->tpltzblocks[0].lambda == 1) /* Use straightforward loop for white noise model */
    {
        // Here it is assumed that we use a single bandwidth for all TOD intervals
        // i.e. lambda is the same for all Toeplitz blocks
        for (i = 0; i < Nm1->nb_blocks_loc; i++)
        {
            for (j = 0; j < Nm1->tpltzblocks[i].n; j++)
            {
                tod[t_id + j] = Nm1->tpltzblocks[i].T_block[0] * tod[t_id + j];
            }
            t_id += Nm1->tpltzblocks[i].n;
        }
    }
    // else if (tot_ngap == 0) /* No gaps in the timestream */
    // {
    //     // Use stbmmProd routine for correlated noise model
    //     // (no det-det correlations for now)
    //     stbmmProd(Nm1, tod);
    // }
    // else
    // {
    //     PCG_global(N, Nm1, Gaps, tod, N->local_V_size, rank, (rank == 0));
    // }
    else /* Iteratively solve Nx=b to apply the noise weights */
    {
        // variables for a single toeplitz block
        Block block, block_m1;
        Tpltz N_block, Nm1_block;
        double *tod_block = NULL;

        for (i = 0; i < N->nb_blocks_loc; ++i)
        {
            // pick the single Toeplitz block
            block = N->tpltzblocks[i];
            block_m1 = Nm1->tpltzblocks[i];

            // if (rank == 0)
            // {
            //     printf("process toeplitz block %d of %d\n", i+1, N->nb_blocks_loc);
            //     printf("  block.n = %d\n", block.n);
            //     printf("  block.idv = %ld\n", block.idv);
            //     printf("  offset = %d\n", t_id);
            //     fflush(stdout);
            // }

            // define Tpltz structures for the single block
            set_tpltz_struct(&N_block, N, &block);
            set_tpltz_struct(&Nm1_block, Nm1, &block_m1);

            // pointer to current block in the tod
            tod_block = (tod + t_id);

            // apply the weights with a PCG
            PCG_single_block(&N_block, &Nm1_block, Gaps, tod_block, (rank == 0));

            // update our index of local samples
            t_id += block.n;
        }
    }
}

void set_tpltz_struct(Tpltz *single_block_struct, Tpltz *full_struct, Block *block)
{
    single_block_struct->nrow = full_struct->nrow; // does not matter anyway
    single_block_struct->m_cw = full_struct->m_cw; // 1
    single_block_struct->m_rw = full_struct->m_rw; // 1
    single_block_struct->tpltzblocks = block;
    single_block_struct->nb_blocks_loc = 1;
    single_block_struct->nb_blocks_tot = 1;
    single_block_struct->idp = block->idv;
    single_block_struct->local_V_size = block->n;
    single_block_struct->flag_stgy = full_struct->flag_stgy;

    // TODO: can not set to NULL, maybe create a communicator with only 1 process?
    single_block_struct->comm = full_struct->comm;
}

/**
 * @brief Iteratively solve Nx=b in order to compute an inverse noise-weighted timestream.
 * Here it is assumed that we work with a single toeplitz block.
 *
 * @param N_block [in] Tpltz structure of the single block
 * @param Nm1_block [in] Tpltz structure of the single inverse block
 * @param Gaps [in] structure of the timestream gaps
 * @param tod_block [in] pointer to the RHS vector [out] solution of the system
 * @param m_block [in] size of the block
 * @param verbose
 */
void PCG_single_block(Tpltz *N_block, Tpltz *Nm1_block, Gap *Gaps, double *tod_block, int verbose)
{
    int k = 0;                               // iteration number
    double r0, res, tol2rel;                 // residuals
    double coef_1, coef_2;                   // scalars
    double st, t;                            // timing variables
    bool converged = false;                  // convergence state
    int ng = Gaps->ngap;                     // number of gaps
    int m_block = N_block->tpltzblocks[0].n; // size of the data
    int i;                                   // loop index

    // TODO: implement stopping criterion (divergence)
    // TODO: store solver results in a SolverInfo structure passed as argument

    st = MPI_Wtime();

    // reset gaps of the rhs
    if (ng > 0)
    {
        reset_tod_gaps(tod_block, N_block, Gaps);
    }

    double *_r = malloc((sizeof *_r) * m_block);
    double *r = malloc((sizeof *r) * m_block); // residual

    if (_r == NULL || r == NULL)
    {
        puts("out of memory: allocation of time-domain vectors for tod weighting failed");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m_block; ++i)
    {
        // set initial guess
        // _r[i] = 0.;
        _r[i] = tod_block[i];

        // need rhs to compute first residual
        r[i] = tod_block[i];
        tod_block[i] = _r[i];
    }

    // apply system matrix (_r = Nx0)
    if (ng > 0)
    {
        gstbmmProd(N_block, _r, Gaps);
    }
    else
    {
        stbmmProd(N_block, _r);
    }

    // compute initial residual (r = b - Nx0)
    // check for instant convergence
    res = 0.0; // (r,r)
    for (i = 0; i < m_block; ++i)
    {
        r[i] = r[i] - _r[i];
        res += r[i] * r[i];
    }

    r0 = res;
    tol2rel = PCG_TOEPLITZ_TOL * PCG_TOEPLITZ_TOL * res;

    t = MPI_Wtime();

    if (verbose > 0) /* print information on screen */
    {
        printf("  apply_weights: k = %d, r0 = %e, time=%lf\n", k, r0, t - st);
        fflush(stdout);
    }

    if (res < PCG_TOEPLITZ_TOL * PCG_TOEPLITZ_TOL)
    {
        // print info on screen
        if (verbose > 0)
        {
            printf("  -> converged (%e < %e)\n", res, PCG_TOEPLITZ_TOL * PCG_TOEPLITZ_TOL);
            fflush(stdout);
        }

        // return immediately
        free(_r);
        free(r);
        return;
    }

    double *p = malloc((sizeof *p) * m_block);   // search direction
    double *z = malloc((sizeof *z) * m_block);   // M^{-1} * r
    double *zp = malloc((sizeof *zp) * m_block); // previous z
    double *zt = NULL;                           // backup pointer for zp

    if (p == NULL || z == NULL || zp == NULL)
    {
        puts("out of memory: allocation of time-domain vectors for tod weighting failed");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // apply preconditioner (z0 = M^{-1} * r0)
    for (i = 0; i < m_block; ++i)
    {
        z[i] = r[i];
    }

    if (ng > 0)
    {
        gstbmmProd(Nm1_block, z, Gaps);
    }
    else
    {
        stbmmProd(Nm1_block, z);
    }

    // set initial search direction (p0 = z0)
    for (i = 0; i < m_block; ++i)
    {
        p[i] = z[i];
    }

    for (k = 1; k < PCG_TOEPLITZ_KMAX + 1; ++k)
    {
        // apply system matrix (_r = Np)
        for (i = 0; i < m_block; ++i)
        {
            _r[i] = p[i];
        }

        if (ng > 0)
        {
            gstbmmProd(N_block, _r, Gaps);
        }
        else
        {
            stbmmProd(N_block, _r);
        }

        // compute alpha = (r,z)/(p,Np)
        coef_1 = 0.0; // (r,z)
        coef_2 = 0.0; // (p,Np)
        for (i = 0; i < m_block; ++i)
        {
            coef_1 += r[i] * z[i];
            coef_2 += p[i] * _r[i];
        }
        
        // swap pointers to store previous z before updating
        zt = zp;
        zp = z;
        z = zt;

        // update current vector (x = x + alpha * p)
        // update residual (r = r - alpha * Np)
        for (i = 0; i < m_block; ++i)
        {
            tod_block[i] = tod_block[i] + (coef_1 / coef_2) * p[i];
            r[i] = r[i] - (coef_1 / coef_2) * _r[i];
            z[i] = r[i];
        }

        // apply preconditioner (z = M^{-1} r)
        if (ng > 0)
        {
            gstbmmProd(Nm1_block, z, Gaps);
        }
        else
        {
            stbmmProd(Nm1_block, z);
        }

        // compute coeff for new search direction
        coef_2 = 0.0; // new (r,z)
        for (i = 0; i < m_block; ++i)
        {
            // Fletcher-Reeves
            // coef_2 += r[i] * z[i];

            // Polak-Ribière
            coef_2 += r[i] * (z[i] - zp[i]);
        }

        // update search direction
        res = 0.0; // (r,r)
        for (i = 0; i < m_block; ++i)
        {
            p[i] = z[i] + (coef_2 / coef_1) * p[i];
            res += r[i] * r[i];
        }

        if (verbose > 0)
        {
            printf("  -> k = %d, res = %e, relative = %e\n", k, res, sqrt(res / r0));
            fflush(stdout);
        }

        if (res < tol2rel)
        {
            t = MPI_Wtime();
            converged = true;
            break;
        }
    }

    t = MPI_Wtime();

    if (verbose > 0)
    {
        if (converged)
        {
            printf("  -> converged (%e < %e)\n", res, tol2rel);
            printf("  -> i.e.      (%e < %e)\n", sqrt(res / r0), PCG_TOEPLITZ_TOL);
            printf("  -> n_iter = %d, solve time = %lf s\n", k, t - st);
        }
        else
        {
            printf("  -> no convergence (%e > %e)\n", res, tol2rel);
            printf("  -> i.e.           (%e > %e)\n", sqrt(res / r0), PCG_TOEPLITZ_TOL);
            printf("  -> n_iter = %d, time = %lf s\n", k - 1, t - st);
        }
        fflush(stdout);
    }

    // free memory
    free(p);
    free(_r);
    free(r);
    free(z);
    free(zp);
}

void reset_tod_gaps(double *tod, Tpltz *N, Gap *Gaps)
{
    reset_gaps(&tod, N->idp, N->local_V_size, N->m_cw, N->nrow, N->m_rw, Gaps->id0gap, Gaps->lgap, Gaps->ngap);
}

void PCG_global(Tpltz *N, Tpltz *Nm1, Gap *Gaps, double *tod, int m, int rank, int verbose)
{
    int k = 0;               // iteration number
    double r0, res, tol2rel; // residuals
    double coef_1, coef_2;   // scalars
    double lval_1, lval_2;   // local variables before reduce
    double st, t;            // timing variables
    bool converged = false;  // convergence state
    int ng = Gaps->ngap;
    int i;

    double *p = NULL, *_r = NULL, *r = NULL, *z = NULL; // time-domain objects

    st = MPI_Wtime();

    // reset gaps of the rhs
    if (ng > 0)
        reset_tod_gaps(tod, N, Gaps);

    _r = malloc((sizeof *_r) * m);
    r = malloc((sizeof *r) * m);

    if (_r == NULL || r == NULL)
    {
        puts("out of memory: allocation of time-domain vectors for tod weighting failed");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m; ++i)
    {
        // set initial guess
        // _r[i] = 0.;
        _r[i] = tod[i];

        // need rhs to compute first residual
        r[i] = tod[i];
        tod[i] = _r[i];
    }

    // apply system matrix (_r = Nx0)
    if (ng > 0)
        gstbmmProd(N, _r, Gaps);
    else
        stbmmProd(N, _r);

    // compute initial residual (r = b - Nx0)
    // check for instant convergence
    lval_1 = 0.0;
    for (i = 0; i < m; ++i)
    {
        r[i] = r[i] - _r[i];
        lval_1 += r[i] * r[i];
    }

    res = 0.0;
    MPI_Allreduce(&lval_1, &res, 1, MPI_DOUBLE, MPI_SUM, N->comm);
    r0 = res;
    tol2rel = PCG_TOEPLITZ_TOL * PCG_TOEPLITZ_TOL * res;

    t = MPI_Wtime();

    if ((rank == 0) && (verbose > 0))
    {
        printf("  apply_weights: k = %d, r0 = %e, time=%lf\n", k, r0, t - st);
        fflush(stdout);
    }

    if (res < PCG_TOEPLITZ_TOL * PCG_TOEPLITZ_TOL)
    {
        // print info on screen
        if ((rank == 0) && (verbose > 0))
        {
            printf("  -> converged (%e < %e)\n", res, PCG_TOEPLITZ_TOL * PCG_TOEPLITZ_TOL);
            fflush(stdout);
        }

        // return immediately
        free(_r);
        free(r);
        return;
    }

    p = malloc((sizeof *p) * m);
    z = malloc((sizeof *z) * m);
    if (p == NULL || z == NULL)
    {
        puts("out of memory: allocation of time-domain vectors for tod weighting failed");
        fflush(stdout);
        exit(EXIT_FAILURE);
    }

    // apply preconditioner (z0 = M^{-1} * r0)
    for (i = 0; i < m; ++i)
        z[i] = r[i];

    if (ng > 0)
        gstbmmProd(Nm1, z, Gaps);
    else
        stbmmProd(Nm1, z);

    // set initial search direction (p0 = z0)
    for (i = 0; i < m; ++i)
        p[i] = z[i];

    for (k = 1; k < PCG_TOEPLITZ_KMAX + 1; ++k)
    {
        // apply system matrix (_r = Np)
        for (i = 0; i < m; ++i)
            _r[i] = p[i];

        if (ng > 0)
            gstbmmProd(N, _r, Gaps);
        else
            stbmmProd(N, _r);

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
        MPI_Allreduce(&lval_1, &coef_1, 1, MPI_DOUBLE, MPI_SUM, N->comm);
        MPI_Allreduce(&lval_2, &coef_2, 1, MPI_DOUBLE, MPI_SUM, N->comm);

        // update current vector (x = x + alpha * p)
        // update residual (r = r - alpha * Np)
        for (i = 0; i < m; ++i)
        {
            tod[i] = tod[i] + (coef_1 / coef_2) * p[i];
            r[i] = r[i] - (coef_1 / coef_2) * _r[i];
            z[i] = r[i];
        }

        // apply preconditioner (z = M^{-1} r)
        if (ng > 0)
            gstbmmProd(Nm1, z, Gaps);
        else
            stbmmProd(Nm1, z);

        // compute coeff for new search direction
        lval_2 = 0.0;
        for (i = 0; i < m; ++i)
            lval_2 += r[i] * z[i];

        coef_2 = 0.0; // new (r,z)
        MPI_Allreduce(&lval_2, &coef_2, 1, MPI_DOUBLE, MPI_SUM, N->comm);

        // update search direction
        lval_1 = 0.0;
        for (i = 0; i < m; ++i)
        {
            p[i] = z[i] + (coef_2 / coef_1) * p[i];
            lval_1 += r[i] * r[i];
        }

        res = 0.0;
        MPI_Allreduce(&lval_1, &res, 1, MPI_DOUBLE, MPI_SUM, N->comm);

        if ((rank == 0) && (verbose > 0))
        {
            printf("  -> k = %d, res = %e, relative = %e\n", k, res, sqrt(res / r0));
            fflush(stdout);
        }

        if (res < tol2rel)
        {
            t = MPI_Wtime();
            converged = true;
            if (rank == 0)
            {
                if (verbose > 0)
                {
                    printf("  -> converged (%e < %e)\n", res, tol2rel);
                    printf("  -> i.e.      (%e < %e)\n", sqrt(res / r0), PCG_TOEPLITZ_TOL);
                    printf("  -> n_iter = %d, solve time = %lf s\n", k, t - st);
                }
                else
                {
                    printf("  -> applied noise weights in %d iterations (%lf s)\n", k - 1, t - st);
                }
            }
            break;
        }
    }

    if (!converged && rank == 0)
    {
        t = MPI_Wtime();
        printf("  -> did not converge in %d iterations (%lf s)\n", k - 1, t - st);
    }

    fflush(stdout);

    // free memory
    free(p);
    free(_r);
    free(r);
    free(z);
}
