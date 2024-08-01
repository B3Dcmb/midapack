/**
 * @file mappraiser.c
 * @brief Process pointing, signal and noise data arrays to produce maps in FITS
 * format
 * @authors Hamza El Bouhargani
 * @date May 2019
 * @update June 2020 by Aygul Jamal
 */

#include "precond.h"
#include <fitsio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <mappraiser/create_toeplitz.h>
#include <mappraiser/gap_filling.h>
#include <mappraiser/iofiles.h>
#include <mappraiser/map.h>
#include <mappraiser/mapping.h>
#include <mappraiser/pcg_true.h>
#include <mappraiser/weight.h>
#include <midapack/memutils.h>

#ifdef WITH_ECG
#include <mappraiser/ecg.h>
#endif

WeightStgy handle_gaps(Gap *Gaps, Mat *A, Tpltz *Nm1, Tpltz *N, GapStrategy gs,
                       double *b, const double *noise, bool do_gap_filling,
                       uint64_t realization, const uint64_t *detindxs,
                       const uint64_t *obsindxs, const uint64_t *telescopes,
                       double sample_rate);

void x2map_pol(double *mapI, double *mapQ, double *mapU, double *Cond,
               int *hits, const double *x, const int *lstid, const double *cond,
               const int *lhits, int xsize, int nnz);

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
    double st = MPI_Wtime();

    MatInit(&A, m, Nnz, pix, pixweights, pointing_commflag, comm);
    Gaps.ngap = build_pixel_to_time_domain_mapping(&A);

    MPI_Barrier(comm);
    double elapsed = MPI_Wtime() - st;

    nbr_extra_pixels = A.trash_pix * A.nnz;
    nbr_valid_pixels = A.lcount - nbr_extra_pixels;

    if (rank == 0) {
        printf("Initialized pointing matrix in %lf s\n", elapsed);
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

    cond = SAFEMALLOC(sizeof *cond * solver_map_size / A.nnz);
    lhits = SAFEMALLOC(sizeof *lhits * solver_map_size / A.nnz);

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
    tpltzblocks = SAFEMALLOC(sizeof *tpltzblocks * nb_blocks_loc);
    defineBlocks_avg(tpltzblocks, inv_tt, nb_blocks_loc, local_blocks_sizes,
                     lambda_block_avg, id0);
    defineTpltz_avg(&Nm1, nrow, 1, mcol, tpltzblocks, nb_blocks_loc,
                    nb_blocks_tot, id0, local_V_size, flag_stgy, comm);

    // define the noise covariance matrix
    Tpltz N;
    Block *tpltzblocks_N = SAFEMALLOC(sizeof *tpltzblocks_N * nb_blocks_loc);
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

    st = MPI_Wtime();

    // first build the BJ preconditioner
    Precond *P =
        newPrecondBJ(&A, &Nm1, cond, lhits, gs, &Gaps, gif, local_blocks_sizes);

    // Allocate memory for the map with the right number of pixels
    x = SAFECALLOC(P->n, sizeof *x);

    MPI_Barrier(comm);
    elapsed = MPI_Wtime() - st;

    if (rank == 0) {
        printf("Block Jacobi preconditioner built for %d sky pixels (%d valid "
               "+ %d extra)\n",
               P->n / A.nnz, P->n_valid / A.nnz, P->n_extra / A.nnz);
        printf("Total time = %lf s\n", elapsed);
        fflush(stdout);
    }

    // Guard against using ECG with extra pixels to estimate
    if (P->n_extra > 0 && solver == 1) {
        if (rank == 0) {
            fprintf(stderr,
                    "ECG solver does not support solving for extra pixels. "
                    "Choose another gap strategy, or use solver=0 (PCG).\n");
        }
        exit(EXIT_FAILURE);
    }

    // Gap treatment can happen now

    MPI_Barrier(comm);
    if (rank == 0) {
        puts("##### Gap treatment ####################");
        fflush(stdout);
    }

    WeightStgy ws =
        handle_gaps(&Gaps, &A, &Nm1, &N, gs, signal, noise, do_gap_filling,
                    realization, detindxs, obsindxs, telescopes, sample_rate);

    // final weighting operator
    WeightMatrix W = createWeightMatrix(&Nm1, &N, &Gaps, ws);

    // ____________________________________________________________
    // Now build the 2lvl part of the preconditioner if needed

    if (precond != BJ) {
        MPI_Barrier(comm);
        if (rank == 0) {
            puts("##### 2lvl preconditioner ####################");
            fflush(stdout);
        }

        st = MPI_Wtime();

        P->ptype = precond;
        P->Zn = Z_2lvl == 0 ? size : Z_2lvl;
        buildPrecond2lvl(P, &A, &W, x, signal);

        MPI_Barrier(comm);
        elapsed = MPI_Wtime() - st;

        if (rank == 0) {
            printf("2lvl preconditioner construction took %lf s\n", elapsed);
            fflush(stdout);
        }
    }

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
        PCG_maxL(&A, P, &W, x, signal, &si);

        // Write PCG residuals to disk
        if (rank == 0) {
            char fname[FILENAME_MAX];
            sprintf(fname, "%s/residuals_%s.dat", outpath, ref);
            int info = solverinfo_write(&si, fname);
            if (info != 0) {
                fputs("Problem writing residuals to file", stderr);
            }
        }

        // Free SolverInfo structure
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
    FREE(tpltzblocks);
    FREE(tpltzblocks_N);

    // free Gap structure
    FREE(Gaps.id0gap);
    FREE(Gaps.lgap);

    // free memory allocated for preconditioner
    PrecondFree(P);

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
        double *extra_map = SAFEMALLOC(sizeof *extra_map * extra);
        memcpy(extra_map, x, extra * sizeof(double));
        if (rank == 0) {
            printf("extra map with %d pixels (T only)\n {", extra / Nnz);
            for (int j = 0; j < extra; j += Nnz) {
                printf(" %e", extra_map[j]);
            }
            puts(" }");
        }
        fflush(stdout);
        FREE(extra_map);
#endif
        // valid map
        memmove(x, x + extra, sizeof *x * map_size);
        memmove(lhits, lhits + extra / Nnz, sizeof *lhits * map_size / Nnz);
        memmove(cond, cond + extra / Nnz, sizeof *cond * map_size / Nnz);
        x = SAFEREALLOC(x, sizeof *x * map_size);
        lhits = SAFEREALLOC(lhits, sizeof *lhits * map_size / Nnz);
        cond = SAFEREALLOC(cond, sizeof *cond * map_size / Nnz);
    }

    // get maps from all processes and combine them

    int *lstid = SAFEMALLOC(sizeof *lstid * map_size);
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

        double *mapI = NULL;
        if (Nnz == 3) {
            mapI = SAFECALLOC(npix, sizeof *mapI);
        }
        double *mapQ = SAFECALLOC(npix, sizeof *mapQ);
        double *mapU = SAFECALLOC(npix, sizeof *mapU);
        int *hits = SAFECALLOC(npix, sizeof *hits);
        double *Cond = SAFECALLOC(npix, sizeof *Cond);

        for (i = 0; i < size; i++) {
            if (i != 0) {
                oldsize = map_size;
                MPI_Recv(&map_size, 1, MPI_INT, i, 0, comm, &status);
                if (oldsize != map_size) {
                    lstid = SAFEREALLOC(lstid, sizeof *lstid * map_size);
                    x = SAFEREALLOC(x, sizeof *x * map_size);
                    cond = SAFEREALLOC(cond, sizeof *cond * map_size);
                    lhits = SAFEREALLOC(lhits, sizeof *lhits * map_size);
                }
                MPI_Recv(lstid, map_size, MPI_INT, i, 1, comm, &status);
                MPI_Recv(x, map_size, MPI_DOUBLE, i, 2, comm, &status);
                MPI_Recv(cond, map_size / Nnz, MPI_DOUBLE, i, 3, comm, &status);
                MPI_Recv(lhits, map_size / Nnz, MPI_INT, i, 4, comm, &status);
            }
            x2map_pol(mapI, mapQ, mapU, Cond, hits, x, lstid, cond, lhits,
                      map_size, Nnz);
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

        if (Nnz == 3) {
            sprintf(Imap_name, "%s/mapI_%s.fits", outpath, ref);
            if (access(Imap_name, F_OK) != -1) {
                ret = remove(Imap_name);
                if (ret != 0) {
                    printf("Error: unable to delete the file %s\n", Imap_name);
                    w = 0;
                }
            }
        }

        sprintf(Qmap_name, "%s/mapQ_%s.fits", outpath, ref);
        sprintf(Umap_name, "%s/mapU_%s.fits", outpath, ref);
        sprintf(Condmap_name, "%s/Cond_%s.fits", outpath, ref);
        sprintf(Hitsmap_name, "%s/Hits_%s.fits", outpath, ref);

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
            if (Nnz == 3) {
                write_map(mapI, TDOUBLE, nside, Imap_name, nest, cordsys);
            }
            write_map(mapQ, TDOUBLE, nside, Qmap_name, nest, cordsys);
            write_map(mapU, TDOUBLE, nside, Umap_name, nest, cordsys);
            write_map(Cond, TDOUBLE, nside, Condmap_name, nest, cordsys);
            write_map(hits, TINT, nside, Hitsmap_name, nest, cordsys);
        } else {
            fprintf(stderr, "IO Error: Could not overwrite old files, map "
                            "results will not be stored ;(\n");
        }

        FREE(mapI);
        FREE(mapQ);
        FREE(mapU);
        FREE(Cond);
        FREE(hits);
    }

    elapsed = MPI_Wtime() - st;
    if (rank == 0) {
        printf("Total time = %lf s\n", elapsed);
        fflush(stdout);
    }

    // free memory
    FREE(x);
    FREE(cond);
    FREE(lhits);

    MatFree(&A);
    A.indices = NULL;
    A.values = NULL;
    FREE(lstid);

    // MPI_Finalize();
}

void MTmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond, int Z_2lvl, int pointing_commflag, int npoly, int nhwp, double delta_t,
           int ground,int n_sss_bins, double tol, int maxiter, int enlFac, int ortho_alg, int bs_red, int nside, int **sweeptstamps,
           int *nsweeps, double **az, double *az_min, double *az_max, double **hwp_angle, int nces, void *data_size_proc, int nb_blocks_loc,
           void *local_blocks_sizes, int Nnz, void *pix, void *pixweights, void *signal, double *noise, double sampling_freq, double *invtt)
{
    int64_t M;                          // Global number of rows
    int     m, Nb_t_Intervals;          // local number of rows of the pointing matrix A, nbr
                                        // of stationary intervals
    int64_t    gif;                     // global indice for the first local line
    int        i, j, k, l;
    Mat        A;                       // pointing matrix structure
    int       *id_last_pix, *ll = NULL; // pixel-to-time-domain mapping
    int        nbr_valid_pixels;        // nbr of valid pixel indices
    double    *x, *cond = NULL;         // pixel domain vectors
    int       *lhits = NULL;
    double     st, t;                   // timer, start time
    int        rank, size;
    MPI_Status status;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    #ifdef HAVE_MKL
    mkl_set_num_threads(1); // required for some obscure reason when working with mkl routines
    #endif

    if (rank == 0) {
        printf("\n############# MAPPRAISER : MidAPack PaRAllel Iterative Sky EstimatoR v2.1, July 2024 "
               "################\n");
        printf("Last compiled on %s at %s\n", __DATE__, __TIME__);
        printf("rank = %d, size = %d\n", rank, size);
        fflush(stdout);
    }

    // total length of the time domain signal
    M = 0;
    for (i = 0; i < size; i++) {
        M += ((int *) data_size_proc)[i];
        fflush(stdout);
    }
    if (rank == 0) {
        printf("[rank %d] M=%ld\n", rank, M);
        fflush(stdout);
    }

    // compute distribution indexes over the processes
    m   = ((int *) data_size_proc)[rank];
    gif = 0;
    for (i = 0; i < rank; i++) { gif += ((int *) data_size_proc)[i]; }

    // Print information on data distribution
    int Nb_t_Intervals_loc = nb_blocks_loc;
    MPI_Allreduce(&nb_blocks_loc, &Nb_t_Intervals, 1, MPI_INT, MPI_SUM, comm);
    int nb_proc_shared_one_interval = 1; // max(1, size/Nb_t_Intervals );
    if (rank == 0) {
        printf("[rank %d] size=%d \t m=%d \t Nb_t_Intervals=%d \n", rank, size, m, Nb_t_Intervals);
        printf("[rank %d] Nb_t_Intervals_loc=%d \n", rank, Nb_t_Intervals_loc);
        fflush(stdout);
    }

    // Pointing matrix initialization

    st          = MPI_Wtime();
    A.trash_pix = 0;
    MatInit(&A, m, Nnz, pix, pixweights, pointing_commflag, comm);

    // reduce the size of map-domain objects in case the trash_pix flag is
    // raised
    nbr_valid_pixels = A.lcount;
    if (A.trash_pix) { nbr_valid_pixels -= A.nnz; }

    t = MPI_Wtime();
    if (rank == 0) {
        printf("[rank %d] Initializing pointing matrix time = %lf\n", rank, t - st);
        printf("  -> nbr of sky pixels = %d\n", A.lcount);
        printf("  -> valid sky pixels  = %d\n", nbr_valid_pixels);
        printf("  -> trash_pix flag    = %d\n", A.trash_pix);
        fflush(stdout);
    }

    // Build pixel-to-time-domain mapping

    st = MPI_Wtime();

    // index of last sample pointing to each pixel
    id_last_pix = (int *) malloc(nbr_valid_pixels / (A.nnz) * sizeof(int));

    // linked list of time samples indexes
    ll = (int *) malloc(m * sizeof(int));

    if (id_last_pix == NULL || ll == NULL) {
        printf("memory allocation failed");
        exit(1);
    }

    // initialize the mapping arrays to -1
    for (i = 0; i < m; i++) { ll[i] = -1; }
    for (j = 0; j < nbr_valid_pixels / (A.nnz); j++) { id_last_pix[j] = -1; }

    // build the linked list chain of time samples corresponding to each pixel
    int vpix_i = 0;
    for (i = 0; i < m; i++) {
        // only do something for valid pixels
        if (!A.trash_pix || (A.trash_pix && (A.indices[i * A.nnz] != 0))) {
            vpix_i = A.indices[i * A.nnz] / A.nnz - A.trash_pix;
            if (id_last_pix[vpix_i] == -1) {
                id_last_pix[vpix_i] = i;
            } else {
                ll[i]               = id_last_pix[vpix_i];
                id_last_pix[vpix_i] = i;
            }
        }
    }

    A.id_last_pix = id_last_pix;
    A.ll          = ll;

    t = MPI_Wtime();
    if (rank==0) {
      printf("[rank %d] Total pixel-to-time domain mapping time=%lf \n", rank, t-st);
      fflush(stdout);
    }

    // Map objects memory allocation
    // MatInit gives A.lcount which is used to compute nbr_valid_pixels

    x     = (double *) malloc(nbr_valid_pixels * sizeof(double));
    cond  = (double *) malloc((int) (nbr_valid_pixels / A.nnz) * sizeof(double));
    lhits = (int *) malloc((int) (nbr_valid_pixels / A.nnz) * sizeof(int));
    if (x == NULL || cond == NULL || lhits == NULL) {
        printf("memory allocation failed");
        exit(1);
    }

    for (j = 0; j < nbr_valid_pixels; j++) {
        x[j] = 0.;
        if (j % A.nnz == 0) {
            lhits[(int) (j / A.nnz)] = 0;
            cond[(int) (j / A.nnz)]  = 0.;
        }
    }

    // Create piecewise Toeplitz matrix
    // specifics parameters:
    int nb_blocks_tot    = Nb_t_Intervals;
    int lambda_block_avg = 1;

    // flags for Toeplitz product strategy
    Flag flag_stgy;
    flag_stgy_init_auto(&flag_stgy);

    // to print something on screen
    flag_stgy.flag_verbose = 1;

    // define Toeplitz blocks list and structure for Nm1
    Block *tpltzblocks;
    Tpltz  Nm1;

    // dependants parameters:
    int64_t nrow = M;
    int     mcol = 1;

    int64_t id0          = gif;
    int     local_V_size = m;

    // Block definition
    tpltzblocks = (Block *) malloc(nb_blocks_loc * sizeof(Block));
    defineBlocks_avg(tpltzblocks, invtt, nb_blocks_loc, local_blocks_sizes, lambda_block_avg, id0);
    defineTpltz_avg(&Nm1, nrow, 1, mcol, tpltzblocks, nb_blocks_loc, nb_blocks_tot, id0, local_V_size, flag_stgy, comm);

    // print Toeplitz parameters for information
    if (rank == 0) {
        printf("[rank %d] Noise model: Banded block Toeplitz, half bandwidth = "
               "%d\n",
               rank, lambda_block_avg);
        fflush(stdout);
    }

    // Templates classes initialization
    st=MPI_Wtime();

    int store_hwp = 0;
    int n_class = 0;
    double sigma2;
    hwpss_w hwpss_wghts;
    int ndet = nb_blocks_loc / nces;
    int **detnsweeps = (int **) malloc(nces * sizeof(int*));
    int *ces_length = (int *) malloc(nces * sizeof(int));
    int *hwp_bins = (int *) malloc(nces * sizeof(int));

    // Polynomial templates metadata
    for(i=0;i<nces;i++){
        detnsweeps[i] = (int *) calloc(ndet, sizeof(int));
        for(j=0;j<ndet;j++)
          detnsweeps[i][j] = nsweeps[i];
          ces_length[i] = sweeptstamps[i][nsweeps[i]];
    }

    // Binned boresight azimuth
    int **az_binned;
    az_binned = bin_az(az, az_min, az_max, ces_length, ground, n_sss_bins, nces);

    // HWP harmonics
    double ***hwp_mod = (double ***) malloc(nces * sizeof(double **));
    if(nhwp){
        for(i=0;i<nces;i++){
            hwp_mod[i] = (double **) malloc(2 * sizeof(double *));
            hwp_mod[i][0] = (double *) calloc(ces_length[i], sizeof(double));//hwp_cos[i];
            hwp_mod[i][1] = (double *) calloc(ces_length[i], sizeof(double));//hwp_sin[i];
            for(j=0;j<ces_length[i];j++){
                //hwp_angle_bis = (double)(2*M_PI*hwp_f*j)/sampling_freq;
                hwp_mod[i][0][j] = cos(hwp_angle[i][j]);
                hwp_mod[i][1][j] = sin(hwp_angle[i][j]);
            }
        }
    }

    // Set number of template classes
    n_class = npoly + ground + 2*nhwp;

    // Allocate memory to the templates classes instances
    TemplateClass *X = (TemplateClass *) malloc(n_class * nb_blocks_loc * sizeof(TemplateClass));

    // Initialize templates classes list
    Tlist_init(X, ndet, nces, (int *)local_blocks_sizes, detnsweeps, ces_length,
    sweeptstamps, n_sss_bins, az_binned, sampling_freq, npoly, ground, nhwp, delta_t,
    store_hwp, hwp_mod);

    // Allocate memory for the list of kernel blocks and inv block container
    int global_size_kernel = 0;
    int id_kernelblock = 0;
    for(i=0;i<nces;i++){
      hwp_bins[i] = (int)ceil(ces_length[i]/(delta_t*sampling_freq));
      global_size_kernel += ndet * (npoly*nsweeps[i] + ground*n_sss_bins + 2*nhwp*hwp_bins[i]) * (npoly*nsweeps[i] + ground*n_sss_bins + 2*nhwp*hwp_bins[i]);
    }
    double *B = (double *) calloc(global_size_kernel, sizeof(double));
    double *Binv = (double *) calloc((npoly*nsweeps[0] + ground*n_sss_bins + 2*nhwp*hwp_bins[0])*(npoly*nsweeps[0] + ground*n_sss_bins + 2*nhwp*hwp_bins[0]), sizeof(double));

    // Build the list of inverse kernel blocks
    for(i=0;i<nces;i++){
      if(store_hwp == 0)
        build_hwpss_w(&hwpss_wghts, hwp_mod[i], ces_length[i], nhwp, i);

      if(i!=0){
        Binv = (double *) realloc(Binv, (npoly*nsweeps[i] + ground*n_sss_bins + 2*nhwp*hwp_bins[i])*(npoly*nsweeps[i] + ground*n_sss_bins + 2*nhwp*hwp_bins[i])*sizeof(double));
        // init to zero
        for(k=0;k<(npoly*nsweeps[i] + ground*n_sss_bins + 2*nhwp*hwp_bins[i])*(npoly*nsweeps[i] + ground*n_sss_bins + 2*nhwp*hwp_bins[i]);k++)
          Binv[k] = 0;
      }

      // Processing detector blocks
      for(j=0;j<ndet;j++){
        BuildKernel(X+(i*ndet+j)*n_class, n_class, B+id_kernelblock, Nm1.tpltzblocks[i*ndet+j].T_block[0], sweeptstamps[i], az_binned[i], hwpss_wghts, delta_t, sampling_freq);

        if(rank==0){
          printf("[rank %d] Effective rank of local kernel block %d = %d\n",rank, i*ndet+j, InvKernel(B+id_kernelblock, npoly*nsweeps[i] + ground*n_sss_bins + 2*nhwp*hwp_bins[i], Binv));
          fflush(stdout);
        }
        else
          InvKernel(B+id_kernelblock, npoly*nsweeps[i] + ground*n_sss_bins +  2*nhwp*hwp_bins[i], Binv);
        for(k=0;k<(npoly*nsweeps[i] + ground*n_sss_bins +  2*nhwp*hwp_bins[i])*(npoly*nsweeps[i] + ground*n_sss_bins +  2*nhwp*hwp_bins[i]);k++){
          B[id_kernelblock+k] = Binv[k];
          Binv[k] = 0;
        }
        id_kernelblock += (npoly*nsweeps[i] + ground*n_sss_bins +  2*nhwp*hwp_bins[i])*(npoly*nsweeps[i] + ground*n_sss_bins +  2*nhwp*hwp_bins[i]);
      }
      if(store_hwp == 0)
        free_hwpss_w(&hwpss_wghts, nhwp);
    }

    t=MPI_Wtime();
    if (rank==0) {
      printf("[rank %d] Total time building Templates classes and inverse kernel blocks=%lf \n", rank, t-st);
      fflush(stdout);
    }

    MPI_Barrier(comm);
    if (rank == 0) { printf("##### Start PCG ####################\n"); }
    fflush(stdout);

    st=MPI_Wtime();
    // Conjugate Gradient
    if(solver==0)
        PCG_GLS_templates(outpath, ref, &A, &Nm1, X, B, sweeptstamps, npoly, ground, nhwp, nsweeps, az_binned, n_sss_bins, hwp_bins, hwp_mod, delta_t, store_hwp, nces, ces_length, nb_blocks_loc, x, signal, noise, cond, lhits, tol, maxiter, sampling_freq, precond, Z_2lvl);
    else{
      printf("ECG unavailable at this stage please choose the PCG solver: solver=0\n");
      exit(1);
    }

    MPI_Barrier(comm);
    t = MPI_Wtime();
    if (rank == 0) {
        printf("##### End PCG ####################\n");
        printf("[rank %d] Total PCG time=%lf \n", rank, t - st);
    }
    fflush(stdout);

    // write output to fits files:
    st          = MPI_Wtime();
    int mapsize = A.lcount - (A.nnz) * (A.trash_pix);
    int map_id  = rank;

    int *lstid;
    lstid = (int *) calloc(mapsize, sizeof(int));
    for (i = 0; i < mapsize; i++) { lstid[i] = A.lindices[i + (A.nnz) * (A.trash_pix)]; }

    if (rank != 0) {
        MPI_Send(&mapsize, 1, MPI_INT, 0, 0, comm);
        MPI_Send(lstid, mapsize, MPI_INT, 0, 1, comm);
        MPI_Send(x, mapsize, MPI_DOUBLE, 0, 2, comm);
        MPI_Send(cond, mapsize / Nnz, MPI_DOUBLE, 0, 3, comm);
        MPI_Send(lhits, mapsize / Nnz, MPI_INT, 0, 4, comm);
    }

    if (rank == 0) {
        int npix = 12 * nside * nside;
        int oldsize;

        double *mapI;
        mapI = (double *) calloc(npix, sizeof(double));
        double *mapQ;
        mapQ = (double *) calloc(npix, sizeof(double));
        double *mapU;
        mapU = (double *) calloc(npix, sizeof(double));
        int *hits;
        hits = (int *) calloc(npix, sizeof(int));
        double *Cond;
        Cond = (double *) calloc(npix, sizeof(double));

        for (i = 0; i < size; i++) {
            if (i != 0) {
                oldsize = mapsize;
                MPI_Recv(&mapsize, 1, MPI_INT, i, 0, comm, &status);
                if (oldsize != mapsize) {
                    int    *tmp1, *tmp4;
                    double *tmp2, *tmp3;
                    tmp1 = (int *) realloc(lstid, mapsize * sizeof(int));
                    tmp2 = (double *) realloc(x, mapsize * sizeof(double));
                    tmp3 = (double *) realloc(cond, mapsize * sizeof(double));
                    tmp4 = (int *) realloc(lhits, mapsize * sizeof(int));
                    if (tmp1 == NULL || tmp2 == NULL || tmp3 == NULL || tmp4 == NULL) {
                        fprintf(stderr, "realloc failed while receiving data from proc %d", i);
                        exit(EXIT_FAILURE);
                    } else {
                        lstid = tmp1;
                        x     = tmp2;
                        cond  = tmp3;
                        lhits = tmp4;
                    }
                }
                MPI_Recv(lstid, mapsize, MPI_INT, i, 1, comm, &status);
                MPI_Recv(x, mapsize, MPI_DOUBLE, i, 2, comm, &status);
                MPI_Recv(cond, mapsize / Nnz, MPI_DOUBLE, i, 3, comm, &status);
                MPI_Recv(lhits, mapsize / Nnz, MPI_INT, i, 4, comm, &status);
            }
            x2map_pol(mapI, mapQ, mapU, Cond, hits, x, lstid, cond, lhits, mapsize);
        }
        printf("Checking output directory ... old files will be overwritten\n");
        char  Imap_name[256];
        char  Qmap_name[256];
        char  Umap_name[256];
        char  Condmap_name[256];
        char  Hitsmap_name[256];
        char  nest    = 1;
        char *cordsys = "C";
        int   ret, w = 1;

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

    MatFree(&A);
    A.indices = NULL;
    A.values  = NULL; // free memory
    free(x);
    free(cond);
    free(lhits);
    free(tpltzblocks);
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
               const int *lhits, int xsize, int nnz) {
    for (int i = 0; i < xsize; i++) {
        if (nnz == 3) {
            // I, Q and U maps
            if (i % nnz == 0) {
                mapI[(int)(lstid[i] / nnz)] = x[i];
                hits[(int)(lstid[i] / nnz)] = lhits[(int)(i / nnz)];
                Cond[(int)(lstid[i] / nnz)] = cond[(int)(i / nnz)];
            } else if (i % nnz == 1) {
                mapQ[(int)(lstid[i] / nnz)] = x[i];
            } else {
                mapU[(int)(lstid[i] / nnz)] = x[i];
            }
        } else {
            // only Q and U maps are estimated
            if (i % nnz == 0) {
                mapQ[(int)(lstid[i] / nnz)] = x[i];
                hits[(int)(lstid[i] / nnz)] = lhits[(int)(i / nnz)];
                Cond[(int)(lstid[i] / nnz)] = cond[(int)(i / nnz)];
            } else {
                mapU[(int)(lstid[i] / nnz)] = x[i];
            }
        }
    }
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

    // When not doing gap-filling, set signal in the gaps to zero
    const bool reset_signal_in_gaps = !do_gap_filling;

    switch (gs) {

    case COND:
        // set noise weighting strategy
        ws = BASIC;

        if (my_rank == 0) {
            puts("[Gaps/conditioning] weighting strategy = BASIC");
        }

        // set signal in all gaps to zero
        if (reset_signal_in_gaps) {
            reset_relevant_gaps(b, Nm1, Gaps);
        }

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
        if (reset_signal_in_gaps) {
            reset_relevant_gaps(b, Nm1, Gaps);
        }

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
        if (reset_signal_in_gaps) {
            reset_relevant_gaps(b, Nm1, Gaps);
        }

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
        if (reset_signal_in_gaps) {
            reset_relevant_gaps(b, Nm1, Gaps);
        }

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
