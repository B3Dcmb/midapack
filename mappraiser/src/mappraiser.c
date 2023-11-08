/**
 * @file mappraiser.c
 * @brief Process pointing, signal and noise data arrays to produce maps in FITS
 * format
 * @authors Hamza El Bouhargani
 * @date May 2019
 * @update June 2020 by Aygul Jamal
 */

#include "mappraiser/create_toeplitz.h"
#include "mappraiser/gap_filling.h"
#include "mappraiser/iofiles.h"
#include "mappraiser/map.h"
#include "mappraiser/mapping.h"
#include "mappraiser/noise_weighting.h"
#include "mappraiser/pcg_true.h"

#ifdef WITH_ECG
#include "mappraiser/ecg.h"
#endif

#include <fitsio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

WeightStgy handle_gaps(Gap *Gaps, Mat *A, Tpltz *Nm1, Tpltz *N, GapStrategy gs,
                       double *b, const double *noise, bool do_gap_filling,
                       uint64_t realization, const uint64_t *detindxs,
                       const uint64_t *obsindxs, const uint64_t *telescopes,
                       double sample_rate);

void x2map_pol(double *mapI, double *mapQ, double *mapU, double *Cond,
               int *hits, const double *x, const int *lstid, const double *cond,
               const int *lhits, int xsize);

void MLmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond,
           int Z_2lvl, int pointing_commflag, double tol, int maxiter,
           int enlFac, int ortho_alg, int bs_red, int nside, int gap_stgy,
           bool do_gap_filling, uint64_t realization, void *data_size_proc,
           int nb_blocks_loc, void *local_blocks_sizes, double sample_rate,
           uint64_t *detindxs, uint64_t *obsindxs, uint64_t *telescopes,
           int Nnz, void *pix, void *pixweights, void *signal, double *noise,
           int lambda, double *inv_tt, double *tt) {
    int64_t M;             // Global number of rows
    int m, Nb_t_Intervals; // local number of rows of the pointing matrix A, nbr
                           // of stationary intervals
    int64_t gif;           // global indice for the first local line
    int i;
    Mat A;                   // pointing matrix structure
    int nbr_valid_pixels;    // nbr of valid pixel indices
    int nbr_extra_pixels;    // nbr of extra pixel indices
    Gap Gaps;                // timestream gaps structure
    double *x, *cond = NULL; // pixel domain vectors
    int *lhits = NULL;
    double st, t; // timer, start time
    int rank, size;
    MPI_Status status;

    // mkl_set_num_threads(1); // Circumvent an MKL bug

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (rank == 0) {
        printf("\n############# MAPPRAISER : MidAPack PaRAllel Iterative Sky "
               "EstimatoR vDev, May 2019 "
               "################\n");
        printf("Last compiled on %s at %s\n", __DATE__, __TIME__);
        printf("[MPI info] rank = %d, size = %d\n", rank, size);
        puts("##### Initialization ####################");
        fflush(stdout);
    }

    // total length of the time domain signal
    M = 0;
    for (i = 0; i < size; i++) {
        M += ((int *)data_size_proc)[i];
    }

    // compute distribution indexes over the processes
    m = ((int *)data_size_proc)[rank];
    gif = 0;
    for (i = 0; i < rank; i++) {
        gif += ((int *)data_size_proc)[i];
    }

    // Print information on data distribution
    int Nb_t_Intervals_loc = nb_blocks_loc;
    MPI_Allreduce(&nb_blocks_loc, &Nb_t_Intervals, 1, MPI_INT, MPI_SUM, comm);
    if (rank == 0) {
        printf("[Data] global M = %ld (%d intervals)\n", M, Nb_t_Intervals);
        printf("[Data] local  m = %d (%d intervals)\n", m, Nb_t_Intervals_loc);
        fflush(stdout);
    }

    GapStrategy gs = gap_stgy;

    // Set flag to ignore extra pixels when not marginalizing
    A.flag_ignore_extra = !(gs == MARG_LOCAL_SCAN || gs == MARG_PROC);

    // Create extra pixels according to the chosen strategy
    create_extra_pix(pix, pixweights, Nnz, nb_blocks_loc, local_blocks_sizes,
                     gs);

    if (rank == 0) {
        printf("[Gaps] strategy: ");
        print_gap_stgy(gs);
        fflush(stdout);
    }

    // ____________________________________________________________
    // Pointing matrix initialization + mapping

    MPI_Barrier(comm);
    st = MPI_Wtime();

    MatInit(&A, m, Nnz, pix, pixweights, pointing_commflag, comm);
    Gaps.ngap = build_pixel_to_time_domain_mapping(&A);

    MPI_Barrier(comm);
    t = MPI_Wtime();

    nbr_extra_pixels = A.trash_pix * A.nnz;
    nbr_valid_pixels = A.lcount - nbr_extra_pixels;

    if (rank == 0) {
        printf("Initialized pointing matrix in %lf s\n", t - st);
        printf("[proc %d] sky pixels = %d", rank, A.lcount / A.nnz);
        printf(" (%d valid + %d extra)\n", nbr_valid_pixels / A.nnz,
               nbr_extra_pixels / A.nnz);
        printf("[proc %d] local timestream gaps = %d\n", rank, Gaps.ngap);
        fflush(stdout);
    }

    // ____________________________________________________________
    // Map objects memory allocation

    // Size of map that will be estimated by the solver
    int solver_map_size = get_actual_map_size(&A);

    x = calloc(solver_map_size, sizeof *x);
    cond = malloc((sizeof *cond) * solver_map_size / A.nnz);
    lhits = malloc((sizeof *lhits) * solver_map_size / A.nnz);

    if (x == NULL || cond == NULL || lhits == NULL) {
        fprintf(stderr, "[rank %d] memory allocation of map objects failed",
                rank);
        exit(1);
    }

    // ____________________________________________________________
    // Create piecewise Toeplitz matrix

    // specifics parameters:
    int nb_blocks_tot = Nb_t_Intervals;
    int lambda_block_avg = lambda;

    // flags for Toeplitz product strategy
    Flag flag_stgy;
    flag_stgy_init_auto(&flag_stgy);

    // skip build gappy blocks
    flag_stgy.flag_skip_build_gappy_blocks = 1;

    // to print something on screen
    // flag_stgy.flag_verbose = 1;

    // define Toeplitz blocks list and structure for Nm1
    Block *tpltzblocks;
    Tpltz Nm1;

    // dependants parameters:
    int64_t nrow = M;
    int mcol = 1;

    int64_t id0 = gif;
    int local_V_size = m;

    // Block definition
    tpltzblocks = (Block *)malloc(nb_blocks_loc * sizeof(Block));
    defineBlocks_avg(tpltzblocks, inv_tt, nb_blocks_loc, local_blocks_sizes,
                     lambda_block_avg, id0);
    defineTpltz_avg(&Nm1, nrow, 1, mcol, tpltzblocks, nb_blocks_loc,
                    nb_blocks_tot, id0, local_V_size, flag_stgy, comm);

    // define the noise covariance matrix
    Tpltz N;
    Block *tpltzblocks_N = (Block *)malloc(nb_blocks_loc * sizeof(Block));
    defineBlocks_avg(tpltzblocks_N, tt, nb_blocks_loc, local_blocks_sizes,
                     lambda_block_avg, id0);
    defineTpltz_avg(&N, nrow, 1, mcol, tpltzblocks_N, nb_blocks_loc,
                    nb_blocks_tot, id0, local_V_size, flag_stgy, comm);

    // print Toeplitz parameters for information
    if (rank == 0) {
        printf("Noise model: Banded block Toeplitz");
        printf(" (half bandwidth = %d)\n", lambda_block_avg);
        fflush(stdout);
    }

    // ____________________________________________________________
    // Compute the system preconditioner

    MPI_Barrier(comm);
    if (rank == 0) {
        puts("##### Preconditioner ####################");
        fflush(stdout);
    }

    Precond *P = NULL;
    double *pixpond;

    st = MPI_Wtime();

    if (Z_2lvl == 0)
        Z_2lvl = size;

    bool use_precond;
    if (precond < 0) {
        // we won't be using preconditioning in the solver
        use_precond = false;
        precond = 0;
    } else {
        use_precond = true;
    }

    build_precond(&P, &pixpond, &A, &Nm1, &x, signal, noise, &cond, &lhits, tol,
                  Z_2lvl, precond, gs, &Gaps, gif, local_blocks_sizes);

    MPI_Barrier(A.comm);
    t = MPI_Wtime();

    if (rank == 0) {
        printf("Precond built for %d sky pixels (%d valid + %d extra)\n",
               P->n / A.nnz, P->n_valid / A.nnz, P->n_extra / A.nnz);
        printf("Total time = %lf s\n", t - st);
        fflush(stdout);
    }

    // Guard against using ECG in some cases
    if (P->n_extra > 0 && solver == 1) {
        if (rank == 0) {
            fprintf(stderr,
                    "ECG solver does not support solving for extra pixels. "
                    "Choose another gap strategy, or use solver=0 (PCG).\n");
        }
        exit(EXIT_FAILURE);
    }

    //____________________________________________________________
    // Gap treatment

    MPI_Barrier(comm);
    if (rank == 0) {
        puts("##### Gap treatment ####################");
        fflush(stdout);
    }

    WeightStgy ws =
        handle_gaps(&Gaps, &A, &Nm1, &N, gs, signal, noise, do_gap_filling,
                    realization, detindxs, obsindxs, telescopes, sample_rate);

    // ____________________________________________________________
    // Solve the system

    MPI_Barrier(comm);
    if (rank == 0) {
        puts("##### Main solver ####################");
        fflush(stdout);
    }

    if (solver == 0) {
        // set up SolverInfo structure
        SolverInfo si;
        solverinfo_set_defaults(&si);
        si.store_hist = true;
        si.print = rank == 0;
        si.rel_res_reduct = tol;
        si.max_steps = maxiter;

        // solve the equation
        if (use_precond) {
#ifdef DEBUG
            if (rank == 0) {
                puts("call PCG routine (with preconditioning)");
                fflush(stdout);
            }
#endif

            PCG_mm(&A, P, &Nm1, &N, ws, &Gaps, x, signal, &si);
        } else {
#ifdef DEBUG
            if (rank == 0) {
                puts("call CG routine (no preconditioning)");
                fflush(stdout);
            }
#endif
            CG_mm(&A, P->pixpond, &Nm1, &N, ws, &Gaps, x, signal, &si);
        }

        solverinfo_free(&si);

    } else if (solver == 1) {
#ifdef WITH_ECG
        ECG_GLS(outpath, ref, &A, &Nm1, &(P->BJ_inv), P->pixpond, x, signal,
                noise, tol, maxiter, enlFac, ortho_alg, bs_red);
#else
        if (rank == 0)
            fprintf(stderr,
                    "The choice of solver is 1 (=ECG), but the ECG source "
                    "file has not been compiled.\n");
        exit(EXIT_FAILURE);
#endif

    } else {
        if (rank == 0) {
            char msg[] =
                "Incorrect solver parameter. Reminder: solver=0 -> PCG, "
                "solver=1 -> ECG.";
            fputs(msg, stderr);
        }
        exit(EXIT_FAILURE);
    }

    // free tpltz blocks
    free(tpltzblocks);
    free(tpltzblocks_N);

    // free Gap structure
    free(Gaps.id0gap);
    free(Gaps.lgap);

    // free preconditioner structure
    free_precond(&P);

    // ____________________________________________________________
    // Write output to fits files

    MPI_Barrier(comm);
    if (rank == 0) {
        puts("##### Write products ####################");
        fflush(stdout);
    }

    st = MPI_Wtime();

    // throw away estimated extra pixels if there are any

    int map_size = get_valid_map_size(&A);
    int extra = get_actual_map_size(&A) - map_size;

    if (extra > 0) {
#ifdef DEBUG
        double *extra_map = (double *)malloc(extra * sizeof(double));
        memcpy(extra_map, x, extra * sizeof(double));
        if (rank == 0) {
            printf("extra map with %d pixels (T only)\n {", extra / Nnz);
            for (int j = 0; j < extra; j += A.nnz) {
                printf(" %e", extra_map[j]);
            }
            puts(" }");
        }
        fflush(stdout);
#endif
        // valid map
        memmove(x, (x + extra), map_size * sizeof(double));
        memmove(lhits, lhits + extra / Nnz, (map_size / Nnz) * sizeof(int));
        memmove(cond, cond + extra / Nnz, (map_size / Nnz) * sizeof(double));
        double *tmp_x = realloc(x, (sizeof tmp_x) * map_size);
        int *tmp_hits = realloc(lhits, (sizeof tmp_hits) * map_size / Nnz);
        double *tmp_cond = realloc(cond, (sizeof tmp_cond) * map_size / Nnz);
        if (tmp_x == NULL || tmp_hits == NULL || tmp_cond == NULL) {
            fprintf(stderr, "[proc %d] realloc of x, lhits or cond failed",
                    rank);
            exit(EXIT_FAILURE);
        }
        x = tmp_x;
        lhits = tmp_hits;
        cond = tmp_cond;
    }

    // get maps from all processes and combine them

    int *lstid = malloc((sizeof lstid) * map_size);
    if (lstid == NULL) {
        fprintf(stderr, "[proc %d] memory allocation of lstid failed", rank);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < map_size; i++) {
        lstid[i] = A.lindices[i + Nnz * A.trash_pix];
    }

    if (rank != 0) {
        MPI_Send(&map_size, 1, MPI_INT, 0, 0, comm);
        MPI_Send(lstid, map_size, MPI_INT, 0, 1, comm);
        MPI_Send(x, map_size, MPI_DOUBLE, 0, 2, comm);
        MPI_Send(cond, map_size / Nnz, MPI_DOUBLE, 0, 3, comm);
        MPI_Send(lhits, map_size / Nnz, MPI_INT, 0, 4, comm);
    }

    if (rank == 0) {
        int npix = 12 * nside * nside;
        int oldsize;

        double *mapI;
        mapI = (double *)calloc(npix, sizeof(double));
        double *mapQ;
        mapQ = (double *)calloc(npix, sizeof(double));
        double *mapU;
        mapU = (double *)calloc(npix, sizeof(double));
        int *hits;
        hits = (int *)calloc(npix, sizeof(int));
        double *Cond;
        Cond = (double *)calloc(npix, sizeof(double));

        for (i = 0; i < size; i++) {
            if (i != 0) {
                oldsize = map_size;
                MPI_Recv(&map_size, 1, MPI_INT, i, 0, comm, &status);
                if (oldsize != map_size) {
                    int *tmp1, *tmp4;
                    double *tmp2, *tmp3;
                    tmp1 = (int *)realloc(lstid, map_size * sizeof(int));
                    tmp2 = (double *)realloc(x, map_size * sizeof(double));
                    tmp3 = (double *)realloc(cond, map_size * sizeof(double));
                    tmp4 = (int *)realloc(lhits, map_size * sizeof(int));
                    if (tmp1 == NULL || tmp2 == NULL || tmp3 == NULL ||
                        tmp4 == NULL) {
                        fprintf(
                            stderr,
                            "realloc failed while receiving data from proc %d",
                            i);
                        exit(EXIT_FAILURE);
                    } else {
                        lstid = tmp1;
                        x = tmp2;
                        cond = tmp3;
                        lhits = tmp4;
                    }
                }
                MPI_Recv(lstid, map_size, MPI_INT, i, 1, comm, &status);
                MPI_Recv(x, map_size, MPI_DOUBLE, i, 2, comm, &status);
                MPI_Recv(cond, map_size / Nnz, MPI_DOUBLE, i, 3, comm, &status);
                MPI_Recv(lhits, map_size / Nnz, MPI_INT, i, 4, comm, &status);
            }
            x2map_pol(mapI, mapQ, mapU, Cond, hits, x, lstid, cond, lhits,
                      map_size);
        }
        puts("Checking output directory... old files will be overwritten");
        char Imap_name[FILENAME_MAX];
        char Qmap_name[FILENAME_MAX];
        char Umap_name[FILENAME_MAX];
        char Condmap_name[FILENAME_MAX];
        char Hitsmap_name[FILENAME_MAX];
        char nest = 1;
        char *cordsys = "C";
        int ret, w = 1;

        sprintf(Imap_name, "%s/mapI_%s.fits", outpath, ref);
        sprintf(Qmap_name, "%s/mapQ_%s.fits", outpath, ref);
        sprintf(Umap_name, "%s/mapU_%s.fits", outpath, ref);
        sprintf(Condmap_name, "%s/Cond_%s.fits", outpath, ref);
        sprintf(Hitsmap_name, "%s/Hits_%s.fits", outpath, ref);

        if (access(Imap_name, F_OK) != -1) {
            ret = remove(Imap_name);
            if (ret != 0) {
                printf("Error: unable to delete the file %s\n", Imap_name);
                w = 0;
            }
        }

        if (access(Qmap_name, F_OK) != -1) {
            ret = remove(Qmap_name);
            if (ret != 0) {
                printf("Error: unable to delete the file %s\n", Qmap_name);
                w = 0;
            }
        }

        if (access(Umap_name, F_OK) != -1) {
            ret = remove(Umap_name);
            if (ret != 0) {
                printf("Error: unable to delete the file %s\n", Umap_name);
                w = 0;
            }
        }

        if (access(Condmap_name, F_OK) != -1) {
            ret = remove(Condmap_name);
            if (ret != 0) {
                printf("Error: unable to delete the file %s\n", Condmap_name);
                w = 0;
            }
        }

        if (access(Hitsmap_name, F_OK) != -1) {
            ret = remove(Hitsmap_name);
            if (ret != 0) {
                printf("Error: unable to delete the file %s\n", Hitsmap_name);
                w = 0;
            }
        }

        if (w == 1) {
            printf("Writing HEALPix maps FITS files to %s...\n", outpath);
            write_map(mapI, TDOUBLE, nside, Imap_name, nest, cordsys);
            write_map(mapQ, TDOUBLE, nside, Qmap_name, nest, cordsys);
            write_map(mapU, TDOUBLE, nside, Umap_name, nest, cordsys);
            write_map(Cond, TDOUBLE, nside, Condmap_name, nest, cordsys);
            write_map(hits, TINT, nside, Hitsmap_name, nest, cordsys);
        } else {
            fprintf(stderr, "IO Error: Could not overwrite old files, map "
                            "results will not be stored ;(\n");
        }

        free(mapI);
        free(mapQ);
        free(mapU);
        free(Cond);
        free(hits);
    }

    t = MPI_Wtime();
    if (rank == 0) {
        printf("Total time = %lf s\n", t - st);
        fflush(stdout);
    }

    // free memory
    free(x);
    free(cond);
    free(lhits);

    MatFree(&A);
    A.indices = NULL;
    A.values = NULL;
    free(lstid);

    // MPI_Finalize();
}

WeightStgy handle_gaps(Gap *Gaps, Mat *A, Tpltz *Nm1, Tpltz *N, GapStrategy gs,
                       double *b, const double *noise, bool do_gap_filling,
                       uint64_t realization, const uint64_t *detindxs,
                       const uint64_t *obsindxs, const uint64_t *telescopes,
                       double sample_rate) {
    int my_rank;
    MPI_Comm_rank(A->comm, &my_rank);

    compute_gaps_per_block(Gaps, Nm1->nb_blocks_loc, Nm1->tpltzblocks);
    copy_gap_info(Nm1->nb_blocks_loc, Nm1->tpltzblocks, N->tpltzblocks);

#if 0
    if (my_rank == 0) {
        puts("gap informations");
        for (int i = 0; i < Nm1->nb_blocks_loc; ++i) {
            printf("block %d: first %d last %d\n", i,
                   Nm1->tpltzblocks[i].first_gap, Nm1->tpltzblocks[i].last_gap);
        }
        fflush(stdout);
    }

    MPI_Barrier(A->comm);
#endif

    WeightStgy ws;

    switch (gs) {

    case COND:
        // set noise weighting strategy
        ws = BASIC;

        if (my_rank == 0) {
            puts("[Gaps/conditioning] weighting strategy = BASIC");
        }

        // set signal in all gaps to zero
        reset_relevant_gaps(b, Nm1, Gaps);

        // this is not needed any more
        // condition_extra_pix_zero(A);

        // recombine signal and noise
        for (int i = 0; i < A->m; ++i) {
            b[i] += noise[i];
        }

        if (do_gap_filling) {
            perform_gap_filling(A->comm, N, Nm1, b, Gaps, realization, detindxs,
                                obsindxs, telescopes, sample_rate, true);
        } else {
            // perfect noise reconstruction
            if (my_rank == 0) {
                puts("[Gaps/conditioning] perfect noise reconstruction");
            }
        }

        break;

    case MARG_LOCAL_SCAN:
        // set noise weighting strategy
        ws = BASIC;

        if (my_rank == 0) {
            puts("[Gaps/marginalization] weighting strategy = BASIC");
        }

        // set signal in all gaps to zero
        reset_relevant_gaps(b, Nm1, Gaps);

        // recombine signal and noise
        for (int i = 0; i < A->m; ++i) {
            b[i] += noise[i];
        }

        if (do_gap_filling) {
            perform_gap_filling(A->comm, N, Nm1, b, Gaps, realization, detindxs,
                                obsindxs, telescopes, sample_rate, true);
        } else {
            // perfect noise reconstruction
            if (my_rank == 0) {
                puts("[Gaps/marginalization] perfect noise reconstruction");
            }
        }

        break;

    case NESTED_PCG:
        // set noise weighting strategy
        ws = ITER;

        if (my_rank == 0) {
            puts("[Gaps/nested] weighting strategy = ITER");
        }

        // recombine signal and noise
        for (int i = 0; i < A->m; ++i) {
            b[i] += noise[i];
        }

        break;

    case NESTED_PCG_NO_GAPS:
        // set noise weighting strategy
        ws = ITER_IGNORE;

        if (my_rank == 0) {
            puts("[Gaps/nested-ignore] weighting strategy = ITER_IGNORE");
        }

        // set signal in all gaps to zero
        reset_relevant_gaps(b, Nm1, Gaps);

        // recombine signal and noise
        for (int i = 0; i < A->m; ++i) {
            b[i] += noise[i];
        }

        if (do_gap_filling) {
            perform_gap_filling(A->comm, N, Nm1, b, Gaps, realization, detindxs,
                                obsindxs, telescopes, sample_rate, true);
        } else {
            // perfect noise reconstruction
            if (my_rank == 0) {
                puts("[Gaps/nested-ignore] perfect noise reconstruction");
            }
        }

        break;

    case MARG_PROC:
        // set noise weighting strategy
        ws = BASIC;

        if (my_rank == 0) {
            puts("[Gaps/marginalization] weighting strategy = BASIC");
        }

        // set signal in all gaps to zero
        reset_relevant_gaps(b, Nm1, Gaps);

        // recombine signal and noise
        for (int i = 0; i < A->m; ++i) {
            b[i] += noise[i];
        }

        if (do_gap_filling) {
            perform_gap_filling(A->comm, N, Nm1, b, Gaps, realization, detindxs,
                                obsindxs, telescopes, sample_rate, true);
        } else {
            // perfect noise reconstruction
            if (my_rank == 0) {
                puts("[Gaps/marginalization] perfect noise reconstruction");
            }
        }

        break;
    }
    fflush(stdout);

    return ws;
}

// FIXME handle nnz != 3
void x2map_pol(double *mapI, double *mapQ, double *mapU, double *Cond,
               int *hits, const double *x, const int *lstid, const double *cond,
               const int *lhits, int xsize) {
    for (int i = 0; i < xsize; i++) {
        if (i % 3 == 0) {
            mapI[(int)(lstid[i] / 3)] = x[i];
            hits[(int)(lstid[i] / 3)] = lhits[(int)(i / 3)];
            Cond[(int)(lstid[i] / 3)] = cond[(int)(i / 3)];
        } else if (i % 3 == 1)
            mapQ[(int)(lstid[i] / 3)] = x[i];
        else
            mapU[(int)(lstid[i] / 3)] = x[i];
    }
}
