/**
 * @file mappraiser.c
 * @brief Process pointing, signal and noise data arrays to produce maps in
 * FITS format
 * @authors Hamza El Bouhargani
 * @date May 2019
 * @update June 2020 by Aygul Jamal
 */

#include "fitsio.h"
#include "mappraiser/create_toeplitz.h"
#include "mappraiser/iofiles.h"
#include "mappraiser/pcg_true.h"
#include <errno.h>
#include <math.h>
#include <midapack.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

int x2map_pol(double *mapI, double *mapQ, double *mapU, double *Cond, int *hits,
              int npix, double *x, int *lstid, double *cond, int *lhits,
              int xsize);

void MLmap(MPI_Comm comm, char *outpath, char *ref, int solver, int precond,
           int Z_2lvl, int pointing_commflag, double tol, int maxiter,
           int enlFac, int ortho_alg, int bs_red, int nside,
           void *data_size_proc, int nb_blocks_loc, void *local_blocks_sizes,
           int Nnz, void *pix, void *pixweights, void *signal, double *noise,
           int lambda, double *invtt) {
    int64_t M;             // Global number of rows
    int m, Nb_t_Intervals; // local number of rows of the pointing matrix A, nbr
                           // of stationary intervals
    int64_t    gif;        // global indice for the first local line
    int        i, j, k;
    Mat        A;          // pointing matrix structure
    int       *id_last_pix, *ll = NULL; // pixel-to-time-domain mapping
    int        nbr_valid_pixels;        // nbr of valid pixel indices
    double    *x, *cond = NULL;         // pixel domain vectors
    int       *lhits = NULL;
    double     st, t;                   // timer, start time
    int        rank, size;
    MPI_Status status;

    // mkl_set_num_threads(1); // Circumvent an MKL bug

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (rank == 0) {
        printf("\n############# MAPPRAISER : MidAPack PaRAllel Iterative Sky "
               "EstimatoR vDev, May 2019 ################\n");
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
        printf("[rank %d] size=%d \t m=%d \t Nb_t_Intervals=%d \n", rank, size,
               m, Nb_t_Intervals);
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
        printf("[rank %d] Initializing pointing matrix time = %lf\n", rank,
               t - st);
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
    if (rank == 0) {
        printf("[rank %d] Total pixel-to-time-domain mapping time=%lf \n", rank,
               t - st);
    }

    // Map objects memory allocation
    // MatInit gives A.lcount which is used to compute nbr_valid_pixels

    x     = (double *) malloc(nbr_valid_pixels * sizeof(double));
    cond  = (double *) malloc((int) (nbr_valid_pixels / 3) * sizeof(double));
    lhits = (int *) malloc((int) (nbr_valid_pixels / 3) * sizeof(int));
    if (x == NULL || cond == NULL || lhits == NULL) {
        printf("memory allocation failed");
        exit(1);
    }

    for (j = 0; j < nbr_valid_pixels; j++) {
        x[j] = 0.;
        if (j % 3 == 0) {
            lhits[(int) (j / 3)] = 0;
            cond[(int) (j / 3)]  = 0.;
        }
    }

    // Create piecewise Toeplitz matrix
    // specifics parameters:
    int nb_blocks_tot    = Nb_t_Intervals;
    int lambda_block_avg = lambda;

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
    defineBlocks_avg(tpltzblocks, invtt, nb_blocks_loc, local_blocks_sizes,
                     lambda_block_avg, id0);
    defineTpltz_avg(&Nm1, nrow, 1, mcol, tpltzblocks, nb_blocks_loc,
                    nb_blocks_tot, id0, local_V_size, flag_stgy, comm);

    // print Toeplitz parameters for information
    if (rank == 0) {
        printf("[rank %d] Noise model: Banded block Toeplitz, half bandwidth = "
               "%d\n",
               rank, lambda_block_avg);
    }

    MPI_Barrier(comm);
    if (rank == 0) { printf("##### Start PCG ####################\n"); }
    fflush(stdout);

    st = MPI_Wtime();
    // Conjugate Gradient
    if (solver == 0) {
        PCG_GLS_true(outpath, ref, &A, &Nm1, x, signal, noise, cond, lhits, tol,
                     maxiter, precond, Z_2lvl);
    } else if (solver == 1) {
#ifdef WITH_ECG
        ECG_GLS(outpath, ref, &A, &Nm1, x, signal, noise, cond, lhits, tol,
                maxiter, enlFac, ortho_alg, bs_red);
#else
        if (rank == 0)
            printf("To use solver=1 (ECG), use configure option '--with "
                   "ecg'.\n");
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

    // write output to fits files:
    st          = MPI_Wtime();
    int mapsize = A.lcount - (A.nnz) * (A.trash_pix);
    int map_id  = rank;

    int *lstid;
    lstid = (int *) calloc(mapsize, sizeof(int));
    for (i = 0; i < mapsize; i++) {
        lstid[i] = A.lindices[i + (A.nnz) * (A.trash_pix)];
    }

    if (rank != 0) {
        MPI_Send(&mapsize, 1, MPI_INT, 0, 0, comm);
        MPI_Send(lstid, mapsize, MPI_INT, 0, 1, comm);
        MPI_Send(x, mapsize, MPI_DOUBLE, 0, 2, comm);
        MPI_Send(cond, mapsize / Nnz, MPI_DOUBLE, 0, 3, comm);
        MPI_Send(lhits, mapsize / Nnz, MPI_INT, 0, 4, comm);
    }

    if (rank == 0) {
        int npix = 12 * pow(nside, 2);
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
                    lstid = (int *) realloc(lstid, mapsize * sizeof(int));
                    x     = (double *) realloc(x, mapsize * sizeof(double));
                    cond  = (double *) realloc(cond, mapsize * sizeof(double));
                    lhits = (int *) realloc(lhits, mapsize * sizeof(int));
                }
                MPI_Recv(lstid, mapsize, MPI_INT, i, 1, comm, &status);
                MPI_Recv(x, mapsize, MPI_DOUBLE, i, 2, comm, &status);
                MPI_Recv(cond, mapsize / Nnz, MPI_DOUBLE, i, 3, comm, &status);
                MPI_Recv(lhits, mapsize / Nnz, MPI_INT, i, 4, comm, &status);
            }
            x2map_pol(mapI, mapQ, mapU, Cond, hits, npix, x, lstid, cond, lhits,
                      mapsize);
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

int x2map_pol(double *mapI, double *mapQ, double *mapU, double *Cond, int *hits,
              int npix, double *x, int *lstid, double *cond, int *lhits,
              int xsize) {

    int i;

    for (i = 0; i < xsize; i++) {
        if (i % 3 == 0) {
            mapI[(int) (lstid[i] / 3)] = x[i];
            hits[(int) (lstid[i] / 3)] = lhits[(int) (i / 3)];
            Cond[(int) (lstid[i] / 3)] = cond[(int) (i / 3)];
        } else if (i % 3 == 1)
            mapQ[(int) (lstid[i] / 3)] = x[i];
        else
            mapU[(int) (lstid[i] / 3)] = x[i];
    }

    return 0;
}
