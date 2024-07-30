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
#include "memutils.h"

#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

double scalar_prod_reduce(MPI_Comm comm, int n, const double *pixpond,
                          const double *p1, const double *p2,
                          const double *p_sub) {
    double result;

    // contributions from valid pixels
    double local_sum = 0.0;
    for (int i = 0; i < n; i++) {
        if (p_sub == NULL) {
            local_sum += p1[i] * p2[i] * pixpond[i];
        } else {
            local_sum += p1[i] * (p2[i] - p_sub[i]) * pixpond[i];
        }
    }

    // reduce the results for ALL pixels
    MPI_Allreduce(&local_sum, &result, 1, MPI_DOUBLE, MPI_SUM, comm);

    return result;
}

void build_rhs(Mat *A, Tpltz *Nm1, Tpltz *N, Gap *G, WeightStgy ws,
               const double *d, double *rhs) {
    // Build the right-hand side of the equation
    // rhs = A^t * N^{-1} * d

    // number of local samples
    int m = A->m;

    // time domain vector
    double *_t = SAFEMALLOC(sizeof *_t * m);

    // multiply by N^{-1}
    for (int i = 0; i < m; i++)
        _t[i] = d[i];
    apply_weights(Nm1, N, G, _t, ws, false);

    // multiply by A^t
    TrMatVecProd(A, _t, rhs, 0);

    FREE(_t);
}

void opmm(Mat *A, Tpltz *Nm1, Tpltz *N, Gap *G, WeightStgy ws, double *x,
          double *y) {
    // Apply system matrix, i.e. compute
    // y = A^t * N^{-1} * A * x

    // time domain vector
    double *_t = SAFEMALLOC(sizeof *_t * A->m);

    MatVecProd(A, x, _t, 0);
    apply_weights(Nm1, N, G, _t, ws, false);
    TrMatVecProd(A, _t, y, 0);

    FREE(_t);
}

#if 0
double maxabs(double *x, int n) {
    double r;
    double fmax = fabs(x[0]);
    for (int j = 1; j < n; j++) {
        r = fabs(x[j]);
        if (r > fmax)
            fmax = r;
    }
    return fmax;
}

void print_residual_max(double *x, int n_extra, int n_valid) {
    printf("  extra: max = %lf\n", maxabs(x, n_extra));
    printf("  valid: max = %lf\n", maxabs(x + n_extra, n_valid));
}
#endif

/**
 * @brief Solve the map making equation with a preconditioned conjugate gradient
 * algorithm.
 *
 * @param A pointing matrix
 * @param M preconditioner
 * @param Nm1 inverse noise covariance
 * @param N noise covariance
 * @param G timestream gaps
 * @param x [in] starting map [out] estimated solution
 * @param d data vector
 * @param si SolverInfo structure
 */
void PCG_mm(Mat *A, Precond *M, Tpltz *Nm1, Tpltz *N, WeightStgy ws, Gap *G,
            double *x, const double *d, SolverInfo *si) {
    // MPI information
    int rank, size;
    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);

    // initialize the SolverInfo struct
    solverinfo_init(si);

    // store starting time
    MPI_Barrier(A->comm);
    si->start_time = MPI_Wtime();

    int k = 0;                      // iteration number
    int n = get_actual_map_size(A); // map size

    double res;            // norm of residual
    double coef_1, coef_2; // scalars
    double wtime;          // timing variable

    bool stop = false; // stop iteration or continue

    double *rhs = NULL; // right-hand side
    double *r = NULL;   // residual
    double *z = NULL;   // M^{-1} * r
    double *p = NULL;   // search direction
    double *Sp = NULL;  // [system matrix] * p
    double *zp = NULL;  // previous z
    double *zt = NULL;  // backup pointer for zp

    rhs = SAFEMALLOC(sizeof *rhs * n);
    r = SAFEMALLOC(sizeof *r * n);

    // right-hand side
    build_rhs(A, Nm1, N, G, ws, d, rhs);

    // compute initial residual
    opmm(A, Nm1, N, G, ws, x, r);
    for (int i = 0; i < n; i++) {
        r[i] = rhs[i] - r[i]; // r = b - A * x0
    }
    res = scalar_prod_reduce(A->comm, M->n, M->pixpond, r, r, NULL);

    // update SolverInfo
    MPI_Barrier(A->comm);
    wtime = MPI_Wtime();
    solverinfo_update(si, &stop, k, res, wtime);

    if (stop) {
        // one stop condition already met, no iteration
        solverinfo_finalize(si);
        FREE(rhs);
        FREE(r);
        return;
    }

    // iterate to solve the system

    // allocate buffers needed for the iteration
    p = SAFEMALLOC(sizeof *p * n);
    Sp = SAFEMALLOC(sizeof *Sp * n);
    z = SAFEMALLOC(sizeof *z * n);
    zp = SAFEMALLOC(sizeof *zp * n);

    // apply preconditioner (z0 = M^{-1} * r0)
    apply_precond(M, A, r, z);

    // set initial search direction (p0 = z0)
    for (int i = 0; i < n; i++) {
        p[i] = z[i];
    }

    double alpha, beta;

    // iteration loop
    while (!stop) {
        // we are doing one more iteration step
        k++;

        // apply system matrix
        opmm(A, Nm1, N, G, ws, p, Sp);

        // compute (r,z) and (p,_p)
        coef_1 = scalar_prod_reduce(A->comm, M->n, M->pixpond, r, z, NULL);
        coef_2 = scalar_prod_reduce(A->comm, M->n, M->pixpond, p, Sp, NULL);

        alpha = coef_1 / coef_2;

        // swap pointers to store previous z before updating
        zt = zp;
        zp = z;
        z = zt;

        // update current vector (x = x + alpha * p)
        // update residual (r = r - alpha * _p)
        for (int i = 0; i < n; i++) {
            x[i] = x[i] + alpha * p[i];
            r[i] = r[i] - alpha * Sp[i];
        }

        if (si->use_exact_residual) {
            // compute explicit residual
#if 0
            printf("\nimplicit r_%d\n", k);
            print_residual_max(r, M->n_extra, M->n_valid);
#endif
            opmm(A, Nm1, N, G, ws, x, r);
            for (int i = 0; i < n; i++) {
                r[i] = rhs[i] - r[i]; // r = b - A * x_k
            }
#if 0
            printf("explicit r_%d\n", k);
            print_residual_max(r, M->n_extra, M->n_valid);
#endif
        }

        // apply preconditioner (z = M^{-1} * r)
        apply_precond(M, A, r, z);

        // compute updated (r,z)
        // use Polak-Ribière formula (r,z-zp)
        coef_2 = scalar_prod_reduce(A->comm, M->n, M->pixpond, r, z, zp);

        beta = coef_2 / coef_1;

        // take beta = max(beta, 0)
        if (beta < 0) {
            beta = 0;
        }

        // update search direction
        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }

        // compute residual
        res = scalar_prod_reduce(A->comm, M->n, M->pixpond, r, r, NULL);

        // update SolverInfo structure
        // and check stop conditions
        MPI_Barrier(A->comm);
        wtime = MPI_Wtime();
        solverinfo_update(si, &stop, k, res, wtime);
    }

    // free memory after the iteration
    FREE(rhs);
    FREE(r);
    FREE(p);
    FREE(Sp);
    FREE(z);
    FREE(zp);

    solverinfo_finalize(si);
}

/**
 * @brief Solve the map making equation with a simple conjugate gradient
 * algorithm.
 *
 * @param A pointing matrix
 * @param pixpond pixel share ponderation
 * @param Nm1 inverse noise covariance
 * @param N noise covariance
 * @param G timestream gaps
 * @param x [in] starting map [out] estimated solution
 * @param d data vector
 * @param si SolverInfo structure
 */
void CG_mm(Mat *A, double *pixpond, Tpltz *Nm1, Tpltz *N, WeightStgy ws, Gap *G,
           double *x, const double *d, SolverInfo *si) {
    // MPI information
    int rank, size;
    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);

    // initialize the SolverInfo struct
    solverinfo_init(si);

    // store starting time
    MPI_Barrier(A->comm);
    si->start_time = MPI_Wtime();

    int k = 0;                      // iteration number
    int n = get_actual_map_size(A); // map size

    double res;            // norm of residual
    double coef_1, coef_2; // scalars
    double wtime;          // timing variable

    bool stop = false; // stop iteration or continue

    double *rhs = NULL; // right-hand side
    double *r = NULL;   // residual
    double *p = NULL;   // search direction
    double *Sp = NULL;  // [system matrix] * p

    rhs = SAFEMALLOC(sizeof *rhs * n);
    r = SAFEMALLOC(sizeof *r * n);

    // right-hand side
    build_rhs(A, Nm1, N, G, ws, d, r);

    // compute initial residual
    opmm(A, Nm1, N, G, ws, x, r);
    for (int i = 0; i < n; i++) {
        r[i] = rhs[i] - r[i]; // r = b - A * x0
    }
    res = scalar_prod_reduce(A->comm, n, pixpond, r, r, NULL);

    // update SolverInfo
    MPI_Barrier(A->comm);
    wtime = MPI_Wtime();
    solverinfo_update(si, &stop, k, res, wtime);

    if (stop) {
        // one stop condition already met, no iteration
        solverinfo_finalize(si);
        FREE(rhs);
        FREE(r);
        return;
    }

    // iterate to solve the system

    // allocate buffers needed for the iteration
    p = SAFEMALLOC(sizeof *p * n);
    Sp = SAFEMALLOC(sizeof *Sp * n);

    // set initial search direction (p0 = r0)
    for (int i = 0; i < n; i++) {
        p[i] = r[i];
    }

    // iteration loop
    while (!stop) {
        // we are doing one more iteration step
        k++;

        // apply system matrix
        opmm(A, Nm1, N, G, ws, p, Sp);

        // compute (r,r) and (p,_p)
        coef_1 = scalar_prod_reduce(A->comm, n, pixpond, r, r, NULL);
        coef_2 = scalar_prod_reduce(A->comm, n, pixpond, p, Sp, NULL);

        // update current vector (x = x + alpha * p)
        // update residual (r = r - alpha * _p)
        for (int i = 0; i < n; i++) {
            x[i] = x[i] + (coef_1 / coef_2) * p[i];
            r[i] = r[i] - (coef_1 / coef_2) * Sp[i];
        }

        if (si->use_exact_residual) {
            // compute explicit residual
            opmm(A, Nm1, N, G, ws, x, r);
            for (int i = 0; i < n; i++) {
                r[i] = rhs[i] - r[i]; // r = b - A * x_k
            }
        }

        // compute updated (r,r)
        coef_2 = scalar_prod_reduce(A->comm, n, pixpond, r, r, NULL);

        // update search direction
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + (coef_2 / coef_1) * p[i];
        }

        // compute residual
        res = scalar_prod_reduce(A->comm, n, pixpond, r, r, NULL);

        // update SolverInfo structure
        // and check stop conditions
        MPI_Barrier(A->comm);
        wtime = MPI_Wtime();
        solverinfo_update(si, &stop, k, res, wtime);
    }

    // free memory after the iteration
    FREE(rhs);
    FREE(r);
    FREE(p);
    FREE(Sp);

    solverinfo_finalize(si);
}
