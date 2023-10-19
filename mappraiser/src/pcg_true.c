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

#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

double scalar_prod_reduce(MPI_Comm comm, int n, int n_extra,
                          const double *pixpond, const double *p1,
                          const double *p2, const double *p_sub) {
    double result;

    // contributions from valid pixels
    double local_sum = 0.0;
    for (int i = n_extra; i < n; i++) {
        if (p_sub == NULL) {
            local_sum += p1[i] * p2[i] * pixpond[i - n_extra];
        } else {
            local_sum += p1[i] * (p2[i] - p_sub[i]) * pixpond[i - n_extra];
        }
    }

    // reduce the results for valid pixels
    MPI_Allreduce(&local_sum, &result, 1, MPI_DOUBLE, MPI_SUM, comm);

    // contributions from extra pixels
    for (int i = 0; i < n_extra; i++) {
        if (p_sub == NULL) {
            result += p1[i] * p2[i];
        } else {
            result = p1[i] * (p2[i] - p_sub[i]);
        }
    }

    return result;
}

int build_rhs(Mat *A, Tpltz *Nm1, Tpltz *N, Gap *G, WeightStgy ws,
              const double *d, double *x0, double *rhs) {
    // Build the right-hand side of the equation
    // rhs = A^t * N^{-1} * (b - A * x_0 )

    int m = A->m; // number of local samples
    double *_t;   // time domain vector

    _t = malloc((sizeof *_t) * m);
    if (_t == NULL) {
        fprintf(stderr, "malloc of _t failed in get_rhs");
        return 1;
    }

    MatVecProd(A, x0, _t, 0);
    for (int i = 0; i < m; i++)
        _t[i] = d[i] - _t[i];
    apply_weights(Nm1, N, G, _t, ws, false);
    TrMatVecProd(A, _t, rhs, 0);

    free(_t);
    return 0;
}

int opmm(Mat *A, Tpltz *Nm1, Tpltz *N, Gap *G, WeightStgy ws, double *x,
         double *y) {
    // Apply system matrix, i.e. compute
    // y = A^t * N^{-1} * A * x

    double *_t = NULL; // time domain vector

    _t = malloc((sizeof *_t) * A->m);
    if (_t == NULL) {
        fprintf(stderr, "malloc of _t failed in opmm");
        return 1;
    }

    MatVecProd(A, x, _t, 0);
    apply_weights(Nm1, N, G, _t, ws, false);
    TrMatVecProd(A, _t, y, 0);

    free(_t);
    return 0;
}

void saveArrayToFile(const void *array, size_t elementSize, size_t length,
                     const char *filename) {
    FILE *file = fopen(filename, "wb");

    if (file == NULL) {
        perror("Failed to open the file");
        return;
    }

    size_t elements_written = fwrite(array, elementSize, length, file);

    if (elements_written != length) {
        fprintf(stderr, "Error writing to file\n");
    } else {
        printf("Array saved to %s\n", filename);
    }

    fclose(file);
}

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

    int info;
    int k = 0;                      // iteration number
    int n = get_actual_map_size(A); // map size

    double res;            // norm of residual
    double coef_1, coef_2; // scalars
    double wtime;          // timing variable

    bool stop = false; // stop iteration or continue

    double *r = NULL; // residual
    double *z = NULL; // M^{-1} * r
    double *p = NULL; // search direction
    double *_p = NULL;
    double *zp = NULL; // previous z
    double *zt = NULL; // backup pointer for zp

    // allocate first buffer
    r = malloc((sizeof *r) * n);
    if (r == NULL) {
        si->has_failed = true;
        fprintf(stderr, "[proc %d] malloc of r failed", rank);
        solverinfo_finalize(si);
        return;
    }

    // build rhs
    info = build_rhs(A, Nm1, N, G, ws, d, x, r);
    if (info != 0) {
        fprintf(stderr, "build_rhs routine returned a non-zero code (%d)",
                info);
        exit(EXIT_FAILURE);
    }

    // compute initial residual
    res = scalar_prod_reduce(A->comm, M->n, M->n_extra, M->pixpond, r, r, NULL);

    // update SolverInfo
    MPI_Barrier(A->comm);
    wtime = MPI_Wtime();
    solverinfo_update(si, &stop, k, res, wtime);

    if (stop) {
        // one stop condition already met, no iteration
        solverinfo_finalize(si);
        free(r);
        return;
    }

    // iterate to solve the system

    // allocate buffers needed for the iteration
    p = malloc((sizeof *p) * n);
    _p = malloc((sizeof *_p) * n);
    z = malloc((sizeof *z) * n);
    zp = malloc((sizeof *zp) * n);

    if (p == NULL || _p == NULL || z == NULL || zp == NULL) {
        // free memory if possible
        if (p)
            free(p);
        if (_p)
            free(_p);
        if (z)
            free(z);
        if (zp)
            free(zp);

        si->has_failed = true;
        fprintf(stderr, "[proc %d] malloc of p, _p, z or zp failed", rank);
        solverinfo_finalize(si);
        return;
    }

    // apply preconditioner (z0 = M^{-1} * r0)
    apply_precond(M, A, r, z);

    // set initial search direction (p0 = z0)
    for (int i = 0; i < n; i++) {
        p[i] = z[i];
    }

    // iteration loop
    while (!stop) {
        // we are doing one more iteration step
        k++;

#if 0
        if (rank == 0) {
            char filename[FILENAME_MAX];
            sprintf(filename, "/home/sbiquard/test/data/cond/x_%d", k);
            saveArrayToFile(x, sizeof *x, n, filename);
        }
#endif

        // apply system matrix
        info = opmm(A, Nm1, N, G, ws, p, _p);
        if (info != 0) {
            fprintf(stderr, "opmm routine returned a non-zero code (%d)", info);
            exit(EXIT_FAILURE);
        }

        // compute (r,z) and (p,_p)
        coef_1 = scalar_prod_reduce(A->comm, M->n, M->n_extra, M->pixpond, r, z,
                                    NULL);
        coef_2 = scalar_prod_reduce(A->comm, M->n, M->n_extra, M->pixpond, p,
                                    _p, NULL);

        // swap pointers to store previous z before updating
        zt = zp;
        zp = z;
        z = zt;

        // update current vector (x = x + alpha * p)
        // update residual (r = r - alpha * _p)
        for (int i = 0; i < n; i++) {
            x[i] = x[i] + (coef_1 / coef_2) * p[i];
            r[i] = r[i] - (coef_1 / coef_2) * _p[i];
        }

        // apply preconditioner (z = M^{-1} * r)
        apply_precond(M, A, r, z);

        // compute updated (r,z)
        // use Polak-Ribière formula (r,z-zp)
        coef_2 =
            scalar_prod_reduce(A->comm, M->n, M->n_extra, M->pixpond, r, z, zp);

        // update search direction
        for (int i = 0; i < n; i++) {
            p[i] = z[i] + (coef_2 / coef_1) * p[i];
        }

        // compute residual
        res = scalar_prod_reduce(A->comm, M->n, M->n_extra, M->pixpond, r, r,
                                 NULL);

        // update SolverInfo structure
        // and check stop conditions
        MPI_Barrier(A->comm);
        wtime = MPI_Wtime();
        solverinfo_update(si, &stop, k, res, wtime);
    }

    // free memory after the iteration
    free(r);
    free(p);
    free(_p);
    free(z);
    free(zp);

    solverinfo_finalize(si);
}
