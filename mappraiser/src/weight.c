/**
 * @file weight.c
 * @author Simon Biquard
 * @brief Inverse noise weighting of the map-making procedure with an iterative
 * approach
 * @date Nov 2022
 * @last_update Jul 2024
 */

#ifndef NDEBUG
#include <assert.h>
#endif

#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>

#include <mappraiser/mapping.h>
#include <mappraiser/solver_info.h>
#include <mappraiser/weight.h>
#include <midapack/memutils.h>

WeightMatrix createWeightMatrix(Tpltz *Nm1, Tpltz *N, Gap *G, WeightStgy stgy) {
    // assume everything already allocated
    WeightMatrix W;
    W.Nm1 = Nm1;
    W.N = N;
    W.G = G;
    W.stgy = stgy;
    return W;
}

__attribute__((unused)) void reset_tod_gaps(double *tod, Tpltz *N, Gap *Gaps);

void gappy_tpltz_mult(Tpltz *tmat_block, double *tod, Gap *gaps);

void PCG_single_block(Tpltz *N_block, Tpltz *Nm1_block, Gap *Gaps,
                      double *tod_block, const double *x0, SolverInfo *si,
                      bool ignore_gaps);

// void PCG_global ( Tpltz *N, Tpltz *Nm1, Gap *Gaps, double *tod, int m, int
// rank, int verbose );

/**
 * @brief Weight TOD with the adopted noise model.
 * An PCG is used unless the noise model has lambda=1 (white noise).
 *
 * @param W WeightMatrix structure containing the noise model and gap info
 * @param tod data vector
 * @return number of PCG iterations if there is a single data block, 0 otherwise
 */
int applyWeightMatrix(WeightMatrix *W, double *tod) {
    int t_id = 0; // time sample index in local data

    int rank, size;
    MPI_Comm_rank(W->N->comm, &rank);
    MPI_Comm_size(W->N->comm, &size);

    if (W->stgy == BASIC) {
        // assume no gaps (this is the case if we use a gap-filling procedure,
        // for example)

        if (W->Nm1->tpltzblocks[0].lambda == 1) {
            // Use straightforward loop for white noise model

            // Here it is assumed that we use a single bandwidth for all TOD
            // intervals i.e. lambda is the same for all Toeplitz blocks
            for (int i = 0; i < W->Nm1->nb_blocks_loc; i++) {
                for (int j = 0; j < W->Nm1->tpltzblocks[i].n; j++) {
                    tod[t_id + j] =
                        W->Nm1->tpltzblocks[i].T_block[0] * tod[t_id + j];
                }
                t_id += W->Nm1->tpltzblocks[i].n;
            }
        } else {
            // Use stbmmProd routine for correlated noise model
            // (no det-det correlations for now)
            stbmmProd(W->Nm1, tod);
        }
        return 0;
    }

    // strategy is not BASIC
    // --> iteratively solve Nx=b to apply the noise weights

    bool ignore_gaps = false;
    if (W->stgy == ITER) {
        // don't ignore the gaps
        ignore_gaps = false;
    } else if (W->stgy == ITER_IGNORE) {
        // gap-filling procedure was performed
        // we will ignore the gaps but still iterate
        ignore_gaps = true;
    }

    int n_blocks = W->N->nb_blocks_loc;

    // create a SolverInfo structure for the solver parameters and output
    SolverInfo si;
    solverinfo_set_defaults(&si);

    if (n_blocks == 1) {
        PCG_single_block(W->N, W->Nm1, W->G, tod, NULL, &si, ignore_gaps);
        return si.n_iter;
    }

    // variables for a single toeplitz block
    Block block, block_m1;
    Tpltz N_block, Nm1_block;
    double *tod_block = NULL;

    for (int i = 0; i < W->N->nb_blocks_loc; ++i) {
        // pick the single Toeplitz block
        block = W->N->tpltzblocks[i];
        block_m1 = W->Nm1->tpltzblocks[i];

        // if (rank == 0)
        // {
        //     printf("process toeplitz block %d/%d\n", i,
        //     (N->nb_blocks_loc) - 1); printf("  block.n = %d\n",
        //     block.n); printf("  block.idv = %ld\n", block.idv);
        //     printf("  offset = %d\n", t_id);
        //     fflush(stdout);
        // }

        // define Tpltz structures for the single block
        set_tpltz_struct(&N_block, W->N, &block);
        set_tpltz_struct(&Nm1_block, W->Nm1, &block_m1);

        // pointer to current block in the tod
        tod_block = (tod + t_id);

        // apply the weights with a PCG
        PCG_single_block(&N_block, &Nm1_block, W->G, tod_block, NULL, &si,
                         ignore_gaps);

        // do something with solver output
        //
        //

        if (si.store_hist)
            solverinfo_free(&si);

        // update our index of local samples
        t_id += block.n;
    }
    return 0;
}

/// @brief Initialize a dedicated Tpltz structure for a single data block
/// @param single_block_struct [out] Tpltz structure to initialize
/// @param full_struct [in] Tpltz structure of the whole data
/// @param block [in] the Block structure representing the single data block
void set_tpltz_struct(Tpltz *single_block_struct, Tpltz *full_struct,
                      Block *block) {
    single_block_struct->nrow = full_struct->nrow; // does not matter anyway
    single_block_struct->m_cw = full_struct->m_cw; // 1
    single_block_struct->m_rw = full_struct->m_rw; // 1
    single_block_struct->tpltzblocks = block;
    single_block_struct->nb_blocks_loc = 1;
    single_block_struct->nb_blocks_tot = 1;
    single_block_struct->idp = block->idv;
    single_block_struct->local_V_size = block->n;
    single_block_struct->flag_stgy = full_struct->flag_stgy;
    single_block_struct->comm = full_struct->comm;
}

__attribute__((unused)) void reset_tod_gaps(double *tod, Tpltz *N, Gap *Gaps) {
#ifdef DEBUG
    double start = MPI_Wtime();
#endif
    reset_gaps(&tod, N->idp, N->local_V_size, N->m_cw, N->nrow, N->m_rw,
               Gaps->id0gap, Gaps->lgap, Gaps->ngap);
#ifdef DEBUG
    double duration = MPI_Wtime() - start;
    printf(" reset_tod_gaps (size=%d) in %lf s\n", N->local_V_size, duration);
    fflush(stdout);
#endif
}

void gappy_tpltz_mult(Tpltz *tmat_block, double *tod, Gap *gaps) {
#ifndef NDEBUG
    assert(tmat_block->nb_blocks_loc == 1);
#endif
    reset_relevant_gaps(tod, tmat_block, gaps);
    stbmmProd(tmat_block, tod);
    reset_relevant_gaps(tod, tmat_block, gaps);
}

/**
 * @brief Iteratively solve Nx=b in order to compute an inverse noise-weighted
 * timestream. Here it is assumed that we work with a single toeplitz block.
 *
 * @param N_block [in] Tpltz structure of the single block
 * @param Nm1_block [in] Tpltz structure of the single inverse block
 * @param Gaps [in] structure of the timestream gaps
 * @param tod_block [in] pointer to the RHS vector [out] solution of the system
 * @param x0 [in] starting vector; if NULL, the RHS is used as initial guess
 * @param si SolverInfo structure [in] contains solver parameters [out] contains
 * iteration info
 */
void PCG_single_block(Tpltz *N_block, Tpltz *Nm1_block, Gap *Gaps,
                      double *tod_block, const double *x0, SolverInfo *si,
                      bool ignore_gaps) {
#ifndef NDEBUG
    assert(N_block->nb_blocks_loc == 1 && Nm1_block->nb_blocks_loc == 1);
#endif
    // initialize the SolverInfo struct
    solverinfo_init(si);

    // if (si->print)
    //     solverinfo_print(si);

    int k = 0;                         // iteration number
    int ng = Gaps->ngap;               // number of gaps
    int m = N_block->tpltzblocks[0].n; // size of the data

    double res;            // norm of residual
    double coef_1, coef_2; // scalars
    double wtime;          // timing variable

    bool stop = false; // stop iteration or continue

    if (ng == 0)
        ignore_gaps = true;

    double *_r = NULL;
    double *r = NULL;  // residual
    double *p = NULL;  // search direction
    double *z = NULL;  // M^{-1} * r
    double *zp = NULL; // previous z
    double *zt = NULL; // backup pointer for zp

    // store starting time
    si->start_time = MPI_Wtime();

    // reset gaps of the rhs
    // if (!ignore_gaps) reset_tod_gaps(tod_block, N_block, Gaps);
    if (!ignore_gaps)
        reset_relevant_gaps(tod_block, N_block, Gaps);

    _r = SAFEMALLOC(sizeof *_r * m);
    r = SAFEMALLOC(sizeof *r * m);

    for (int i = 0; i < m; ++i) {
        if (x0 != NULL) {
            // use starting vector if provided
            _r[i] = x0[i];
        } else {
            // else, use the RHS as initial guess
            // one might also simply initialize at zero
            _r[i] = tod_block[i];
        }

        // store the RHS, needed to compute the initial residual
        r[i] = tod_block[i];

        // set TOD to starting vector
        tod_block[i] = _r[i];
    }

    // apply system matrix (_r = Nx0)
    if (ignore_gaps) {
        stbmmProd(N_block, _r);
    } else {
        // gstbmmProd(N_block, _r, Gaps);
        gappy_tpltz_mult(N_block, _r, Gaps);
    }

    // compute initial residual (r = b - Nx0)
    res = 0.0; // (r,r)
    for (int i = 0; i < m; ++i) {
        r[i] = r[i] - _r[i];
        res += r[i] * r[i];
    }

    // update SolverInfo structure
    wtime = MPI_Wtime();
    solverinfo_update(si, &stop, k, res, wtime);

    if (stop) {
        // one stop condition already met, no iteration
        // free the allocated memory
        FREE(_r);
        FREE(r);
        solverinfo_finalize(si);
        return;
    }

    // iterate to solve the system

    // allocate buffers needed for the iteration
    p = SAFEMALLOC(sizeof *p * m);
    z = SAFEMALLOC(sizeof *z * m);
    zp = SAFEMALLOC(sizeof *zp * m);

    // apply preconditioner (z0 = M^{-1} * r0)
    for (int i = 0; i < m; ++i)
        z[i] = r[i];

    if (ignore_gaps) {
        stbmmProd(Nm1_block, z);
    } else {
        // gstbmmProd(Nm1_block, z, Gaps);
        gappy_tpltz_mult(Nm1_block, z, Gaps);
    }

    // set initial search direction (p0 = z0)
    for (int i = 0; i < m; ++i)
        p[i] = z[i];

    // iteration loop
    while (!stop) {
        // we are doing one more iteration step
        ++k;

        // apply system matrix (_r = Np)
        for (int i = 0; i < m; ++i)
            _r[i] = p[i];

        if (ignore_gaps)
            stbmmProd(N_block, _r);
        else {
            // gstbmmProd(N_block, _r, Gaps);
            gappy_tpltz_mult(N_block, _r, Gaps);
        }

        // compute alpha = (r,z)/(p,Np)
        coef_1 = 0.0; // (r,z)
        coef_2 = 0.0; // (p,Np)
        for (int i = 0; i < m; ++i) {
            coef_1 += r[i] * z[i];
            coef_2 += p[i] * _r[i];
        }

        // swap pointers to store previous z before updating
        zt = zp;
        zp = z;
        z = zt;

        // update current vector (x = x + alpha * p)
        // update residual (r = r - alpha * Np)
        for (int i = 0; i < m; ++i) {
            tod_block[i] = tod_block[i] + (coef_1 / coef_2) * p[i];
            r[i] = r[i] - (coef_1 / coef_2) * _r[i];
            z[i] = r[i];
        }

        // apply preconditioner (z = M^{-1} r)
        if (ignore_gaps) {
            stbmmProd(Nm1_block, z);
        } else {
            // gstbmmProd(Nm1_block, z, Gaps);
            gappy_tpltz_mult(Nm1_block, z, Gaps);
        }

        // compute coeff for new search direction
        coef_2 = 0.0; // new (r,z)
        for (int i = 0; i < m; ++i) {
            // Fletcher-Reeves
            // coef_2 += r[i] * z[i];

            // Polak-Ribière
            coef_2 += r[i] * (z[i] - zp[i]);
        }

        // update search direction
        res = 0.0; // (r,r)
        for (int i = 0; i < m; ++i) {
            p[i] = z[i] + (coef_2 / coef_1) * p[i];
            res += r[i] * r[i];
        }

        // update SolverInfo structure
        // and check stop conditions
        wtime = MPI_Wtime();
        solverinfo_update(si, &stop, k, res, wtime);
    }

    // free memory after the iteration
    FREE(p);
    FREE(_r);
    FREE(r);
    FREE(z);
    FREE(zp);

    solverinfo_finalize(si);
}
