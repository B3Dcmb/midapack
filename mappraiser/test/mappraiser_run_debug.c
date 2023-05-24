//
// Created by sbiquard on 5/17/23.
//

#include "utils.h"

#include <libgen.h>
#include <mappraiser.h>
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    //____________________________________________________________
    // MPI initialization

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //____________________________________________________________
    // Parameters

    // Get directory of the executable

    char    buffer[1024];
    char   *directory;
    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (len != -1) {
        buffer[len] = '\0';
        directory   = dirname(buffer);
        if (rank == 0) printf("Executable Directory: %s\n", directory);
    } else {
        if (rank == 0) fprintf(stderr, "Failed to retrieve the executable directory.\n");
        MPI_Finalize();
        return 1;
    }

    // Define the data and the output directories

    char datapath[1024];
    strcpy(datapath, directory);
    strcat(datapath, "/../../../mappraiser/test/data_bin/");

    char outpath[1024];
    strcpy(outpath, directory);
    strcat(outpath, "/../../../mappraiser/test/out_ecg_debug");

    // proc 0 will create the output directory
    if (rank == 0) {
        printf("Checking output directory... ");
        if (createDirectory(outpath) != 0) {
            fprintf(stderr, "Cannot write mappraiser products.\n");
            MPI_Finalize();
            return 1;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Run parameters

    char         ref[1024]         = "run0";
    const int    solver            = 1;
    const int    precond           = 0;
    const int    Z_2lvl            = 0;
    const int    pointing_commflag = 6;
    const double tol               = 1e-6;
    const int    maxiter           = 3000;
    const int    enlFac            = 1;
    const int    ortho_alg         = 1;
    const int    bs_red            = 0;
    const int    nside             = 512;
    const int    Nnz               = 3;
    const int    lambda            = 8192;

    //____________________________________________________________
    // Load data

    if (rank == 0) printf("Loading data in %s...\n", datapath);

    MPI_Barrier(MPI_COMM_WORLD);

    // data distribution
    const int nb_blocks_proc[4] = {4, 2, 4, 4};
    int       data_size_proc[4] = {96000, 48000, 96000, 96000};
    const int nb_blocks_loc     = nb_blocks_proc[rank];
    const int nb_samp           = data_size_proc[rank];

    char tmp[1024];
    char filename[1024];

    // local_blocks_sizes
    int local_blocks_sizes[nb_blocks_loc];
    sprintf(filename, "local_blocks_sizes_%d.bin", rank);
    fillArrayFromFile(strcat(strcpy(tmp, datapath), filename), local_blocks_sizes, ARRAY_SIZE(local_blocks_sizes),
                      sizeof(local_blocks_sizes[0]));

    // pixels
    int pix[nb_samp * Nnz];
    sprintf(filename, "pixels_%d.bin", rank);
    fillArrayFromFile(strcat(strcpy(tmp, datapath), filename), pix, ARRAY_SIZE(pix), sizeof(pix[0]));

    // pixweights
    double pixweights[nb_samp * Nnz];
    sprintf(filename, "pixweights_%d.bin", rank);
    fillArrayFromFile(strcat(strcpy(tmp, datapath), filename), pixweights, ARRAY_SIZE(pixweights),
                      sizeof(pixweights[0]));

    // signal
    double signal[nb_samp];
    sprintf(filename, "signal_%d.bin", rank);
    fillArrayFromFile(strcat(strcpy(tmp, datapath), filename), signal, ARRAY_SIZE(signal), sizeof(signal[0]));

    // noise
    double noise[nb_samp];
    sprintf(filename, "noise_%d.bin", rank);
    fillArrayFromFile(strcat(strcpy(tmp, datapath), filename), noise, ARRAY_SIZE(noise), sizeof(noise[0]));

    // noiseless
    // for (int i = 0; i < nb_samp; ++i) noise[i] = 0;

    // invtt
    double invtt[lambda * nb_blocks_loc];
    sprintf(filename, "invtt_%d.bin", rank);
    fillArrayFromFile(strcat(strcpy(tmp, datapath), filename), invtt, ARRAY_SIZE(invtt), sizeof(invtt[0]));

    /*
        for (int proc = 0; proc < size; ++proc) {
            if (proc == rank) {
                printf("[proc %d] local_block_sizes: ", proc);
                printArray(local_blocks_sizes, ARRAY_SIZE(local_blocks_sizes), sizeof(local_blocks_sizes[0]), printInt);
                printf("\n");
                printf("[proc %d] pixels: ", proc);
                printArray(pixels, 15, sizeof(pixels[0]), printInt);
                printf("\n");
                printf("[proc %d] pixweights: ", proc);
                printArray(pixweights, 15, sizeof(pixweights[0]), printDouble);
                printf("\n");
                printf("[proc %d] signal: ", proc);
                printArray(signal, 15, sizeof(signal[0]), printDouble);
                printf("\n");
                printf("[proc %d] noise: ", proc);
                printArray(noise, 15, sizeof(noise[0]), printDouble);
                printf("\n");
                printf("[proc %d] invtt: ", proc);
                printArray(invtt, 10, sizeof(invtt[0]), printDouble);
                printf("\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    */

    //____________________________________________________________
    // Call MLmap

    MLmap(MPI_COMM_WORLD, outpath, ref, solver, precond, Z_2lvl, pointing_commflag, tol, maxiter, enlFac, ortho_alg,
          bs_red, nside, data_size_proc, nb_blocks_loc, local_blocks_sizes, Nnz, pix, pixweights, signal, noise, lambda,
          invtt);

    MPI_Finalize();
    return 0;
}
