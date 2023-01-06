#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include "mappraiser.h"

void reset_tod_gaps(double *tod, Tpltz *N, Gap *Gaps)
{
    reset_gaps(&tod, N->idp, N->local_V_size, N->m_cw, N->nrow, N->m_rw, Gaps->id0gap, Gaps->lgap, Gaps->ngap);
}

/* Weight TOD data according to the adopted noise model */
void apply_weights(Tpltz *Nm1, Tpltz *N, Gap *Gaps, double *tod, int m, int tot_ngap, MPI_Comm comm, int verbose)
{
    int t_id; // time sample index in local data
    int i, j;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (Nm1->tpltzblocks[0].lambda == 1) /* Use straightforward loop for white noise model */
    {
        // Here it is assumed that we use a single bandwidth for all TOD intervals
        // i.e. lambda is the same for all Toeplitz blocks
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
    else if (tot_ngap == 0) /* No gaps in the timestream */
    {
        // Use stbmmProd routine for correlated noise model
        // (no det-det correlations for now)
        stbmmProd(Nm1, tod);
    }
    else /* Iteratively solve Nx=b to apply the noise weights */
    {
        const int kmax = 50;    // maximum iteration count
        const double tol = 1e-6; // convergence criterion
        int k = 0;               // iteration number
        double r0, res, tol2rel; // residuals
        double coef_1, coef_2;   // scalars
        double lval_1, lval_2;   // local values before reduction
        double st, t;            // timing variables
        bool converged = false;

        double *p = NULL, *_r = NULL, *r = NULL, *z = NULL; // time-domain objects

        st = MPI_Wtime();

        // the rhs should have gaps full of 0
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
            _r[i] = 0.;    // set initial guess to zero
            r[i] = tod[i]; // need rhs to compute first residual
            tod[i] = _r[i];
        }

        // apply system matrix (_r = Nx0)
        if (Gaps->ngap > 0)
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
        MPI_Allreduce(&lval_1, &res, 1, MPI_DOUBLE, MPI_SUM, comm);
        r0 = res;
        tol2rel = tol * tol * res;

        t = MPI_Wtime();

        if ((rank == 0) && (verbose > 0)) /* print information on screen */
        {
            printf("  apply_weights: k = %d, r0 = %e, time=%lf\n", k, r0, t - st);
            fflush(stdout);
        }

        if (res < tol * tol)
        {
            // print info on screen
            if ((rank == 0) && (verbose > 0))
            {
                printf("  -> converged (%e < %e)\n", res, tol * tol);
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
        
        if (Gaps->ngap > 0)
            gstbmmProd(Nm1, z, Gaps);
        else
            stbmmProd(Nm1, z);

        // set initial search direction (p0 = z0)
        for (i = 0; i < m; ++i)
            p[i] = z[i];

        for (k = 1; k < kmax; ++k)
        {
            // apply system matrix (_r = Np)
            for (i = 0; i < m; ++i)
                _r[i] = p[i];

            if (Gaps->ngap > 0)
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
            MPI_Allreduce(&lval_1, &coef_1, 1, MPI_DOUBLE, MPI_SUM, comm);
            MPI_Allreduce(&lval_2, &coef_2, 1, MPI_DOUBLE, MPI_SUM, comm);

            // update current vector (x = x + alpha * p)
            // update residual (r = r - alpha * Np)
            for (i = 0; i < m; ++i)
            {
                tod[i] = tod[i] + (coef_1 / coef_2) * p[i];
                r[i] = r[i] - (coef_1 / coef_2) * _r[i];
                z[i] = r[i];
            }

            // apply preconditioner (z = M^{-1} r)
            if (Gaps->ngap > 0)
                gstbmmProd(Nm1, z, Gaps);
            else
                stbmmProd(Nm1, z);

            // compute coeff for new search direction
            lval_2 = 0.0;
            for (i = 0; i < m; ++i)
                lval_2 += r[i] * z[i];

            coef_2 = 0.0; // new (r,z)
            MPI_Allreduce(&lval_2, &coef_2, 1, MPI_DOUBLE, MPI_SUM, comm);

            // update search direction
            lval_1 = 0.0;
            for (i = 0; i < m; ++i)
            {
                p[i] = z[i] + (coef_2 / coef_1) * p[i];
                lval_1 += r[i] * r[i];
            }

            res = 0.0;
            MPI_Allreduce(&lval_1, &res, 1, MPI_DOUBLE, MPI_SUM, comm);

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
                        printf("  -> i.e.      (%e < %e)\n", sqrt(res / r0), tol);
                        printf("  -> n_iter = %d, solve time = %lf s\n", k, t - st);
                    }
                    else
                    {
                        printf("  -> applied noise weights in %d iterations (%lf s)\n", k, t - st);
                    }
                }
                break;
            }
        }

        if (!converged && rank == 0)
        {
            t = MPI_Wtime();
            printf("  -> did not converge in %d iterations (%lf s)\n", k, t - st);
        }

        fflush(stdout);

        // free memory
        free(p);
        free(_r);
        free(r);
        free(z);
    }
}
