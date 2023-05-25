//
// Created by sbiquard on 5/24/23.
//

#include "utils.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <linux/limits.h>
#include <mappraiser.h>
#include <mpi.h>
#include <string>
#include <unistd.h>
#include <vector>

std::string getExecDirectory() {
    std::string directory;
    char        buffer[PATH_MAX];

    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer) - 1);
    if (len != -1) { buffer[len] = '\0'; }

    char *lastSlash = strrchr(buffer, '/');
    if (lastSlash != nullptr) {
        *lastSlash = '\0';
        directory  = buffer;
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
    if (rank == 0) std::cout << "Running executable from " << directory << std::endl;

    // Define the data and the output directories

    std::string data_path   = directory + "/../../../mappraiser/test/data_bin";
    std::string output_path = directory + "/../../../mappraiser/test/debug_out";

    // proc 0 will create the output directory
    if (rank == 0) {
        std::cout << "Checking output directory... ";
        if (createDirectory(output_path.c_str()) != 0) {
            std::cerr << "Cannot write mappraiser products :'(" << std::endl;
            MPI_Finalize();
            return 1;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Run parameters

    const int       solver            = 0;
    const int       precond           = 0;
    const int       Z_2lvl            = 0;
    const int       pointing_commflag = 6;
    const double    tol               = 1e-6;
    const int       maxiter           = 3000;
    const int       enlFac            = 1;
    const int       ortho_alg         = 1;
    const int       bs_red            = 0;
    const int       nside             = 512;
    const int       Nnz               = 3;
    const int       lambda            = 8192;
    const int       gap_stgy          = 3;
    const u_int64_t realization       = 0;

    // bool to fill noise vector with zeros
    // (--> noiseless run but solver will still iterate)
    bool noiseless_run = false;

    //____________________________________________________________
    // Load data

    std::string fname;

    if (rank == 0) std::cout << "Loading data from " << data_path << "..." << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // data distribution

    const std::array<int, 4> nb_blocks_proc = {4, 2, 4, 4};
    std::array<int, 4>       data_size_proc = {412920, 206460, 412920, 412920};
    const int                nb_blocks_loc  = nb_blocks_proc[rank];
    const int                nb_samp        = data_size_proc[rank];

    // local_blocks_sizes

    std::vector<int> local_blocks_sizes(nb_blocks_loc);
    fname = data_path + "/local_blocks_sizes_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), local_blocks_sizes.data(), local_blocks_sizes.size(),
                      sizeof(local_blocks_sizes[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded local_blocks_sizes" << std::endl;

    // pixels

    std::vector<int> pix(nb_samp * Nnz);
    fname = data_path + "/pixels_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), pix.data(), pix.size(), sizeof(pix[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded pixels" << std::endl;

    // pixweights

    std::vector<double> pixweights(nb_samp * Nnz);
    fname = data_path + "/pixweights_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), pixweights.data(), pixweights.size(), sizeof(pixweights[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded pixweights" << std::endl;

    // signal

    std::vector<double> signal(nb_samp);
    fname = data_path + "/signal_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), signal.data(), signal.size(), sizeof(signal[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded signal" << std::endl;

    // noise

    std::vector<double> noise(nb_samp);
    fname = data_path + "/noise_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), noise.data(), noise.size(), sizeof(noise[0]));

    if (noiseless_run) std::fill(noise.begin(), noise.end(), 0);

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded noise" << std::endl;

    // invtt

    std::vector<double> inv_tt(lambda * nb_blocks_loc);
    fname = data_path + "/invtt_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), inv_tt.data(), inv_tt.size(), sizeof(inv_tt[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded inv_tt" << std::endl;

    // tt

    std::vector<double> tt(lambda * nb_blocks_loc);
    fname = data_path + "/tt_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), tt.data(), tt.size(), sizeof(tt[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded tt" << std::endl;

    // detindxs

    std::vector<u_int64_t> detindxs(nb_blocks_loc);
    fname = data_path + "/detindxs_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), detindxs.data(), detindxs.size(), sizeof(detindxs[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded detindxs" << std::endl;

    // obsindxs

    std::vector<u_int64_t> obsindxs(nb_blocks_loc);
    fname = data_path + "/obsindxs_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), obsindxs.data(), obsindxs.size(), sizeof(obsindxs[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded obsindxs" << std::endl;

    // telescopes

    std::vector<u_int64_t> telescopes(nb_blocks_loc);
    fname = data_path + "/telescopes_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), telescopes.data(), telescopes.size(), sizeof(telescopes[0]));

    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if (rank == 0) std::cout << "  loaded telescopes" << std::endl;

    //____________________________________________________________
    // Call MLmap

    char *outpath = (char *) output_path.c_str();
    char *ref     = (char *) "run0";

    MPI_Barrier(MPI_COMM_WORLD);

    MLmap(MPI_COMM_WORLD, outpath, ref, solver, precond, Z_2lvl, pointing_commflag, tol, maxiter, enlFac, ortho_alg,
          bs_red, nside, gap_stgy, realization, data_size_proc.data(), nb_blocks_loc, local_blocks_sizes.data(),
          detindxs.data(), obsindxs.data(), telescopes.data(), Nnz, pix.data(), pixweights.data(), signal.data(),
          noise.data(), lambda, inv_tt.data(), tt.data());

    MPI_Finalize();
    return 0;
}
