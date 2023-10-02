/**
 * @file mappraiser.c
 * @brief Process pointing, signal and noise data arrays to produce maps in FITS
 * format
 * @authors Hamza El Bouhargani
 * @date May 2019
 * @update June 2020 by Aygul Jamal
 */

#include "mappraiser/create_toeplitz.h"
#include "mappraiser/iofiles.h"
#include "mappraiser/map.h"
#include "mappraiser/mapping.h"
#include "mappraiser/pcg_true.h"

#ifdef WITH_ECG
#include "mappraiser/ecg.h"
#endif

#include <fitsio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void x2map_pol(double *mapI, double *mapQ, double *mapU, double *Cond,
               int *hits, const double *x, const int *lstid, const double *cond,
               const int *lhits, int xsize);

void MLmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond,
           int Z_2lvl, int pointing_commflag, double tol, int maxiter,
           int enlFac, int ortho_alg, int bs_red, int nside, int gap_stgy,
           uint64_t realization, void *data_size_proc, int nb_blocks_loc,
           void *local_blocks_sizes, double sample_rate, uint64_t *detindxs,
           uint64_t *obsindxs, uint64_t *telescopes, int Nnz, void *pix,
           void *pixweights, void *signal, double *noise, int lambda,
           double *inv_tt, double *tt) {
    int64_t M;             // Global number of rows
    int m, Nb_t_Intervals; // local number of rows of the pointing matrix A, nbr
                           // of stationary intervals
    int64_t gif;           // global indice for the first local line
    int i, j, k;
    Mat A;                   // pointing matrix structure
    int nbr_valid_pixels;    // nbr of valid pixel indices
    int nbr_extra_pixels;    // nbr of extra pixel indices
    int ngap;                // nbr of timestream gaps
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
        printf("rank = %d, size = %d\n", rank, size);
        fflush(stdout);
    }

    // total length of the time domain signal
    M = 0;
    for (i = 0; i < size; i++) {
        M += ((int *)data_size_proc)[i];
        fflush(stdout);
    }
    if (rank == 0) {
        printf("[rank %d] M=%ld\n", rank, M);
        fflush(stdout);
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
    int nb_proc_shared_one_interval = 1; // max(1, size/Nb_t_Intervals );
    if (rank == 0) {
        printf("[rank %d] size=%d \t m=%d \t Nb_t_Intervals=%d \n", rank, size,
               m, Nb_t_Intervals);
        printf("[rank %d] Nb_t_Intervals_loc=%d \n", rank, Nb_t_Intervals_loc);
        fflush(stdout);
    }

    // Hardcode this for the moment
    ExtraPixStgy stg = MARG_LOCAL_SCAN;

    // Create extra pixels for marginalization
    create_extra_pix(pix, Nnz, nb_blocks_loc, local_blocks_sizes, stg);

    // ____________________________________________________________
    // Pointing matrix initialization

    st = MPI_Wtime();

    A.trash_pix = 0;
    MatInit(&A, m, Nnz, pix, pixweights, pointing_commflag, comm);

#if 0
    printf("rank %d: %d extra pixels\n", rank, A.trash_pix);
    fflush(stdout);
#endif

    nbr_extra_pixels = A.trash_pix * A.nnz;
    nbr_valid_pixels = A.lcount - nbr_extra_pixels;

    t = MPI_Wtime();
    if (rank == 0) {
        printf("[rank %d] Initializing pointing matrix time = %lf\n", rank,
               t - st);
        printf("Treatment of gaps: %d ", stg);
        puts(stg == COND ? "(conditioning)" : "(marginalization)");
        printf("  -> nbr of sky pixels = %d\n", A.lcount);
        printf("  -> valid pixels = %d\n", nbr_valid_pixels);
        printf("  -> extra pixels = %d\n", nbr_extra_pixels);
        fflush(stdout);
    }

    // ____________________________________________________________
    // Build pixel-to-time-domain mapping

    st = MPI_Wtime();

    ngap = build_pixel_to_time_domain_mapping(&A);

    t = MPI_Wtime();
    if (rank == 0) {
        printf("[rank %d] Pixel-to-time-domain mapping time = %lf s\n", rank,
               t - st);
        printf("  -> # of timestream gaps [local]  = %d\n", ngap);
        fflush(stdout);
    }

    // ____________________________________________________________
    // Map objects memory allocation

    // Decide size of map that will be solved
    // We have computed the amount of valid/extra pixels beforehand
    int solver_map_size;

    switch (stg) {
    case COND:
        solver_map_size = nbr_valid_pixels;
        break;
    case MARG_LOCAL_SCAN:
        solver_map_size = nbr_valid_pixels + nbr_extra_pixels;
        break;
    }

    x = (double *)malloc(solver_map_size * sizeof(double));
    cond = (double *)malloc((int)(solver_map_size / A.nnz) * sizeof(double));
    lhits = (int *)malloc((int)(solver_map_size / A.nnz) * sizeof(int));
    if (x == NULL || cond == NULL || lhits == NULL) {
        printf("memory allocation failed");
        exit(1);
    }

    for (j = 0; j < solver_map_size; j++) {
        x[j] = 0.;
        if (j % A.nnz == 0) {
            lhits[(int)(j / A.nnz)] = 0;
            cond[(int)(j / A.nnz)] = 0.;
        }
    }

#if 0
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
        printf("[rank %d] Noise model: Banded block Toeplitz, half bandwidth = "
               "%d\n",
               rank, lambda_block_avg);
    }

    // ____________________________________________________________
    // Solve the system

    MPI_Barrier(comm);
    if (rank == 0) {
        printf("##### Start PCG ####################\n");
    }
    fflush(stdout);

    st = MPI_Wtime();
    // Conjugate Gradient
    if (solver == 0) {
        PCG_GLS_true(outpath, ref, &A, &Nm1, &N, x, signal, noise, cond, lhits,
                     tol, maxiter, precond, Z_2lvl, &Gaps, gif, gap_stgy,
                     realization, detindxs, obsindxs, telescopes, sample_rate);
    } else if (solver == 1) {
#ifdef WITH_ECG
        ECG_GLS(outpath, ref, &A, &Nm1, x, signal, noise, cond, lhits, tol,
                maxiter, enlFac, ortho_alg, bs_red, &Gaps, gif);
#else
        if (rank == 0)
            fprintf(stderr,
                    "The choice of solver is 1 (=ECG), but the ECG source "
                    "file has not been compiled.\n");
        exit(EXIT_FAILURE);
#endif
    } else {
        printf("Incorrect solver parameter.\n");
        printf("Reminder: solver = 0 -> PCG, solver = 1 -> ECG\n");
        exit(EXIT_FAILURE);
    }

    MPI_Barrier(comm);
    t = MPI_Wtime();
    if (rank == 0) {
        printf("##### End PCG ####################\n");
        printf("[rank %d] Total PCG time=%lf \n", rank, t - st);
    }
    fflush(stdout);

    // free tpltz blocks
    free(tpltzblocks);
    free(tpltzblocks_N);

    // free Gap structure
    free(Gaps.id0gap);
    free(Gaps.lgap);
#endif

    // ____________________________________________________________
    // Write output to fits files

    st = MPI_Wtime();

    // size of the estimated map without extra pixels
    int map_size = A.lcount - (A.nnz) * (A.trash_pix);

    switch (stg) {
    case COND:
        /* map only contains valid pixels, nothing more to do */
        break;

    case MARG_LOCAL_SCAN:
        solver_map_size = A.lcount;

        // map contains estimates for extra pixels which we don't want to keep
        int extra = A.nnz * A.trash_pix;
        double *extra_map = (double *)malloc(extra * sizeof(double));
        memcpy(extra_map, x, extra * sizeof(double));

        // valid map
        memmove(x, (x + extra), map_size * sizeof(double));
        double *tmp = realloc(x, map_size * sizeof(double));
        if (tmp == NULL) {
            fprintf(stderr, "Map reallocation failed");
            exit(EXIT_FAILURE);
        }
        x = tmp;

        // TODO what to do with the extra map?
#if 0
        if (rank == 0) {
            printf("extra map with %d pixels\nI component = {", extra);
            for (j = 0; j < extra; j += A.nnz) {
                printf("%lf ", extra_map[j]);
            }
            puts("}");
        }
        fflush(stdout);
#endif
        break;
    }

    int *lstid = (int *)malloc(map_size * sizeof(int));
    for (i = 0; i < map_size; i++) {
        lstid[i] = A.lindices[i + (A.nnz) * (A.trash_pix)];
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
        printf("Checking output directory ... old files will be overwritten\n");
        char Imap_name[256];
        char Qmap_name[256];
        char Umap_name[256];
        char Condmap_name[256];
        char Hitsmap_name[256];
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
            printf("IO Error: Could not overwrite old files, map results will "
                   "not be stored ;(\n");
        }
    }

    t = MPI_Wtime();
    if (rank == 0) {
        printf("[rank %d] Write output files time=%lf \n", rank, t - st);
        fflush(stdout);
    }
    st = MPI_Wtime();

    // free map domain objects
    free(x);
    free(cond);
    free(lhits);

    MPI_Barrier(comm);
    t = MPI_Wtime();
    if (rank == 0) {
        printf("[rank %d] Free memory time=%lf \n", rank, t - st);
        fflush(stdout);
    }
    // MPI_Finalize();
}

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
