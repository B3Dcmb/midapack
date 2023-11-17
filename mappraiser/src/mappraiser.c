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
           int Z_2lvl, double tol, int maxiter, int enlFac, int ortho_alg,
           int bs_red, int nside, int gap_stgy, bool do_gap_filling,
           uint64_t realization, void *data_size_proc, int nb_blocks_loc,
           void *local_blocks_sizes, double sample_rate, uint64_t *detindxs,
           uint64_t *obsindxs, uint64_t *telescopes, int nnz, void *pix,
           void *pixweights, uint8_t *flags, void *signal, double *noise,
           int lambda, double *inv_tt, double *tt) {
    // ____________________________________________________________
    // MPI info

    int rank, size;
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

    // total length of the time domain signal (global number of rows)
    int64_t M = 0;
    for (int i = 0; i < size; i++) {
        M += ((int *)data_size_proc)[i];
    }

    // number of local samples
    int m = ((int *)data_size_proc)[rank];

    // global index of first local row
    int64_t gif = 0;
    for (int i = 0; i < rank; i++) {
        gif += ((int *)data_size_proc)[i];
    }

    // Print information on data distribution
    int Nb_t_Intervals_loc = nb_blocks_loc;
    int Nb_t_Intervals;
    MPI_Allreduce(&nb_blocks_loc, &Nb_t_Intervals, 1, MPI_INT, MPI_SUM, comm);
    if (rank == 0) {
        printf("[Data] global M = %ld (%d intervals)\n", M, Nb_t_Intervals);
        printf("[Data] local  m = %d (%d intervals)\n", m, Nb_t_Intervals_loc);
        fflush(stdout);
    }

    GapStrategy gs = gap_stgy;

    // Create extra pixels according to the chosen strategy
    create_extra_pix(pix, pixweights, nnz, nb_blocks_loc, local_blocks_sizes,
                     gs);

    if (rank == 0) {
        printf("[Gaps] strategy: ");
        print_gap_stgy(gs);
        fflush(stdout);
    }

    // ____________________________________________________________
    // Pointing matrix initialization + mapping

    MPI_Barrier(comm);
    double st = MPI_Wtime();

    Mat A;

    // Set flag to ignore extra pixels when not marginalizing
    A.ignore_extra = !(gs == MARG_LOCAL_SCAN || gs == MARG_PROC);

    MatInit(&A, m, nnz, pix, pixweights, flags, comm);

    Gap Gaps;
    Gaps.ngap = build_pixel_to_time_domain_mapping(&A);

    MPI_Barrier(comm);
    double t = MPI_Wtime();

    if (rank == 0) {
        printf("Initialized pointing matrix in %lf s\n", t - st);
        printf("[proc %d] sky pixels = %d", rank, A.lcount);
        printf(" (%d valid + %d extra)\n", A.lcount - A.trash_pix, A.trash_pix);
        printf("[proc %d] local timestream gaps = %d\n", rank, Gaps.ngap);
        fflush(stdout);
    }

    // ____________________________________________________________
    // Map objects memory allocation

    // Size of map that will be estimated by the solver
    int npix_solve = get_actual_map_size(&A) / nnz;

    double *x = calloc(nnz * npix_solve, sizeof *x);
    double *cond = malloc((sizeof *cond) * npix_solve);
    int *lhits = malloc((sizeof *lhits) * npix_solve);

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
               P->n / nnz, P->n_valid / nnz, P->n_extra / nnz);
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
        si.use_exact_residual = true;

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

    int npix_iqu = get_valid_map_size(&A);
    int extra_iqu = get_actual_map_size(&A) - npix_iqu;
    int npix = npix_iqu / nnz;
    int extra = extra_iqu / nnz;

    if (extra > 0) {
#ifdef DEBUG
        double *extra_map = (double *)malloc(extra_iqu * sizeof(double));
        memcpy(extra_map, x, extra_iqu * sizeof(double));
        if (rank == 0) {
            printf("extra map with %d pixels (T only)\n {", extra);
            for (int j = 0; j < extra; j++) {
                printf(" %e", extra_map[nnz * j]);
            }
            puts(" }");
        }
        fflush(stdout);
        free(extra_map);
#endif
        // valid map
        memmove(x, x + extra_iqu, npix_iqu * sizeof(double));
        memmove(lhits, lhits + extra, npix * sizeof(int));
        memmove(cond, cond + extra, npix * sizeof(double));
        double *tmp_x = realloc(x, (sizeof *tmp_x) * npix_iqu);
        int *tmp_hits = realloc(lhits, (sizeof *tmp_hits) * npix);
        double *tmp_cond = realloc(cond, (sizeof *tmp_cond) * npix);
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

    int *lstid = malloc((sizeof *lstid) * npix);
    if (lstid == NULL) {
        fprintf(stderr, "[proc %d] memory allocation of lstid failed", rank);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < npix; i++) {
        lstid[i] = A.lindices[i + A.trash_pix];
    }

    if (rank != 0) {
        MPI_Send(&npix, 1, MPI_INT, 0, 0, comm);
        MPI_Send(lstid, npix, MPI_INT, 0, 1, comm);
        MPI_Send(x, npix_iqu, MPI_DOUBLE, 0, 2, comm);
        MPI_Send(cond, npix, MPI_DOUBLE, 0, 3, comm);
        MPI_Send(lhits, npix, MPI_INT, 0, 4, comm);
    }

    if (rank == 0) {
        int npix_fullsky = 12 * nside * nside;

        double *mapI;
        mapI = (double *)calloc(npix_fullsky, sizeof(double));
        double *mapQ;
        mapQ = (double *)calloc(npix_fullsky, sizeof(double));
        double *mapU;
        mapU = (double *)calloc(npix_fullsky, sizeof(double));
        int *hits;
        hits = (int *)calloc(npix_fullsky, sizeof(int));
        double *Cond;
        Cond = (double *)calloc(npix_fullsky, sizeof(double));

        for (int i = 0; i < size; i++) {
            if (i != 0) {
                int npix_p = npix;
                MPI_Recv(&npix, 1, MPI_INT, i, 0, comm, MPI_STATUS_IGNORE);
                npix_iqu = nnz * npix;
                if (npix_p != npix) {
                    int *tmp1, *tmp4;
                    double *tmp2, *tmp3;
                    tmp1 = (int *)realloc(lstid, npix * sizeof(int));
                    tmp2 = (double *)realloc(x, npix_iqu * sizeof(double));
                    tmp3 = (double *)realloc(cond, npix * sizeof(double));
                    tmp4 = (int *)realloc(lhits, npix * sizeof(int));
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
                MPI_Recv(lstid, npix, MPI_INT, i, 1, comm, MPI_STATUS_IGNORE);
                MPI_Recv(x, npix_iqu, MPI_DOUBLE, i, 2, comm,
                         MPI_STATUS_IGNORE);
                MPI_Recv(cond, npix, MPI_DOUBLE, i, 3, comm, MPI_STATUS_IGNORE);
                MPI_Recv(lhits, npix, MPI_INT, i, 4, comm, MPI_STATUS_IGNORE);
            }
            x2map_pol(mapI, mapQ, mapU, Cond, hits, x, lstid, cond, lhits,
                      npix);
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
    int p; // sky pixel index
    for (int i = 0; i < xsize; i++) {
        p = lstid[i];
        mapI[p] = x[3 * i];
        mapQ[p] = x[3 * i + 1];
        mapU[p] = x[3 * i + 2];
        hits[p] = lhits[i];
        Cond[p] = cond[i];
    }
}
