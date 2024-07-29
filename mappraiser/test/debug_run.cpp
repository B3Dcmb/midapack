//
// Created by sbiquard on 5/24/23.
//

#include "utils.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <linux/limits.h>
#include <mappraiser.h>
#include <mpi.h>
#include <string>
#include <unistd.h>
#include <vector>

std::string getExecDirectory() {
    std::string directory;
    char buffer[PATH_MAX];

    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (len != -1) {
        buffer[len] = '\0';
    }

    char *lastSlash = strrchr(buffer, '/');
    if (lastSlash != nullptr) {
        *lastSlash = '\0';
        directory = buffer;
    }

    return directory;
}

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

    std::string directory = getExecDirectory();
    if (rank == 0)
        std::cout << "Running executable from " << directory << "\n";

    // Define the data and the output directories

    std::string data_path = directory + "/data_bin";
    std::string output_path = directory + "/debug_out";

    // proc 0 will create the output directory
    if (rank == 0) {
        std::cout << "Checking output directory... ";
        if (create_directory(output_path.c_str(), 0755) != 0) {
            std::cerr << "Cannot write mappraiser products :'(\n";
            MPI_Finalize();
            return 1;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Run parameters

    int solver = 0;
    int precond = 0;
    int Z_2lvl = 0;
    int pointing_commflag = 2;
    double tol = 1e-6;
    int maxiter = 100;
    int enlFac = 1;
    int ortho_alg = 1;
    int bs_red = 0;
    int nside = 512;
    int gap_stgy = 1;
    bool do_gap_filling = false;
    uint64_t realization = 0;
    int Nnz = 3;
    int lambda = 8192;
    double sample_rate = 200;

    char *outpath = (char *)output_path.c_str();
    char *ref = (gap_stgy == 0) ? (char *)"cond" : (char *)"marg";

    // bool to fill noise vector with zeros
    // (--> noiseless run but solver will still iterate)
    bool fill_noise_zero = false;

    // bool to trigger true noiseless mode
    // (--> set noise to zero + set lambda to 1)
    bool noiseless = false;

    // modify input indices to mimic the observation of a single sky pixel
    bool single_pixel = false;

    // only run mapmaker on a single data block
    // note: only possible with 1 MPI process
    bool single_block = false;

    //____________________________________________________________
    // Process options

    // noiseless mode activate
    if (noiseless) {
        fill_noise_zero = true;
        lambda = 1;
    }

    // force single process if single_det is activated
    if (single_block && size > 1) {
        if (rank == 0)
            std::cerr
                << "single_block mode only works with a single MPI process\n";
        exit(EXIT_FAILURE);
    }

    //____________________________________________________________
    // Load data

    std::string fname;

    if (rank == 0)
        std::cout << "Loading data from " << data_path << "\n";

    MPI_Barrier(MPI_COMM_WORLD);

    // data distribution

    int nb_blocks_loc;
    int nb_samp;
    std::vector<int> data_size_proc;

    if (!single_block) {
        const std::array<int, 4> nb_blocks_proc = {4, 2, 4, 4};
        data_size_proc = {412920, 206460, 412920, 412920};
        nb_blocks_loc = nb_blocks_proc[rank];
        nb_samp = data_size_proc[rank];
    } else {
        // we have a single MPI process
        data_size_proc = {103230};
        nb_blocks_loc = 1;
        nb_samp = data_size_proc[rank];
    }

    // local_blocks_sizes

    std::vector<int> local_blocks_sizes(nb_blocks_loc);
    fname = data_path + "/local_blocks_sizes_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), local_blocks_sizes.data(),
                      local_blocks_sizes.size(), sizeof(local_blocks_sizes[0]));

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded local_blocks_sizes\n";
#endif

    // pixels

    std::vector<int> pix(nb_samp * Nnz);
    fname = data_path + "/pixels_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), pix.data(), pix.size(), sizeof(pix[0]));

    if (single_pixel) {
        for (int i = 0; i < nb_samp; i++) {
            int innz = i * Nnz;
            for (int j = 0; j < Nnz; j++) {
                // don't remove the gaps
                if (pix[innz + j] >= 0) {
                    pix[innz + j] = pix[j];
                }
            }
        }
    }

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded pixels\n";
#endif

    // pixweights

    std::vector<double> pixweights(nb_samp * Nnz);
    fname = data_path + "/pixweights_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), pixweights.data(), pixweights.size(),
                      sizeof(pixweights[0]));

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded pixweights\n";
#endif

    // signal

    std::vector<double> signal(nb_samp);
    fname = data_path + "/signal_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), signal.data(), signal.size(),
                      sizeof(signal[0]));

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded signal\n";
#endif

    // noise

    std::vector<double> noise(nb_samp);
    fname = data_path + "/noise_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), noise.data(), noise.size(),
                      sizeof(noise[0]));

    if (fill_noise_zero) {
        std::fill(noise.begin(), noise.end(), 0);
    }

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded noise\n";
#endif

    // invtt

    std::vector<double> inv_tt(lambda * nb_blocks_loc);
    if (!noiseless) {
        fname = data_path + "/invtt_" + std::to_string(rank) + ".bin";
        fillArrayFromFile(fname.c_str(), inv_tt.data(), inv_tt.size(),
                          sizeof(inv_tt[0]));
    } else {
        std::fill(inv_tt.begin(), inv_tt.end(), 1);
    }

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded inv_tt\n";
#endif

    // tt

    std::vector<double> tt(lambda * nb_blocks_loc);
    if (!noiseless) {
        fname = data_path + "/tt_" + std::to_string(rank) + ".bin";
        fillArrayFromFile(fname.c_str(), tt.data(), tt.size(), sizeof(tt[0]));
    } else {
        std::fill(tt.begin(), tt.end(), 1);
    }

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded tt\n";
#endif

    // detindxs

    std::vector<uint64_t> detindxs(nb_blocks_loc);
    fname = data_path + "/detindxs_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), detindxs.data(), detindxs.size(),
                      sizeof(detindxs[0]));

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded detindxs\n";
#endif

    // obsindxs

    std::vector<uint64_t> obsindxs(nb_blocks_loc);
    fname = data_path + "/obsindxs_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), obsindxs.data(), obsindxs.size(),
                      sizeof(obsindxs[0]));

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded obsindxs\n";
#endif

    // telescopes

    std::vector<uint64_t> telescopes(nb_blocks_loc);
    fname = data_path + "/telescopes_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), telescopes.data(), telescopes.size(),
                      sizeof(telescopes[0]));

#if 0
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "  loaded telescopes\n";
#endif

    //____________________________________________________________
    // Call MLmap

    MPI_Barrier(MPI_COMM_WORLD);

#if 1
    MLmap(MPI_COMM_WORLD, outpath, ref, solver, precond, Z_2lvl,
          pointing_commflag, tol, maxiter, enlFac, ortho_alg, bs_red, nside,
          gap_stgy, do_gap_filling, realization, data_size_proc.data(),
          nb_blocks_loc, local_blocks_sizes.data(), sample_rate,
          detindxs.data(), obsindxs.data(), telescopes.data(), Nnz, pix.data(),
          pixweights.data(), signal.data(), noise.data(), lambda, inv_tt.data(),
          tt.data());
#else
    gap_filling(MPI_COMM_WORLD, data_size_proc.data(), nb_blocks_loc,
                local_blocks_sizes.data(), Nnz, tt.data(), inv_tt.data(),
                lambda, noise.data(), pix.data(), realization, detindxs.data(),
                obsindxs.data(), telescopes.data(), sample_rate);
#endif

    MPI_Finalize();
    return 0;
}
