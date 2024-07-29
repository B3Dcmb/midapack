//
// Created by sbiquard on 5/24/23.
//

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <linux/limits.h>
#include <mappraiser.h>
#include <mpi.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "utils.h"

template <typename T>
bool write_vector_to_file(const std::string &filename,
                          const std::vector<T> &data) {
    // Open the file in binary write mode
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false; // Error opening file
    }

    // Write the data to the file
    file.write(reinterpret_cast<const char *>(data.data()),
               sizeof(T) * data.size());

    // Close the file
    file.close();
    return true; // Success
}

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

    if (size > 1) {
        exit(EXIT_FAILURE);
    }

    //____________________________________________________________
    // Parameters

    // Get directory of the executable
    std::string directory = getExecDirectory();
    std::cout << "Running executable from " << directory << "\n";

    // Define the data and the output directories
    std::string data_path = directory + "/../../../mappraiser/test/data_bin";
    std::string output_path = directory + "/../../../mappraiser/test/debug_out";

    // proc 0 will create the output directory
    std::cout << "Checking output directory... ";
    if (create_directory(output_path.c_str(), 0755) != 0) {
        std::cerr << "Cannot write mappraiser products :'(\n";
        MPI_Finalize();
        return 1;
    }

    //____________________________________________________________
    // Load data

    // data distribution
    const int nb_blocks_tot = 1;
    const int nb_blocks_loc = 1;
    std::vector<int> data_size_proc = {103230};
    const int nb_samp = data_size_proc[rank];
    const int64_t gif = 0;

    std::string fname;

    // local_blocks_sizes
    std::vector<int> local_blocks_sizes(nb_blocks_loc);
    fname = data_path + "/local_blocks_sizes_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), local_blocks_sizes.data(),
                      local_blocks_sizes.size(), sizeof(local_blocks_sizes[0]));

    // signal
    std::vector<double> signal(nb_samp);
    fname = data_path + "/signal_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), signal.data(), signal.size(),
                      sizeof(signal[0]));

    // noise
    std::vector<double> noise(nb_samp);
    fname = data_path + "/noise_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), noise.data(), noise.size(),
                      sizeof(noise[0]));

    // invtt
    int lambda = 8192;
    std::vector<double> inv_tt(lambda * nb_blocks_loc);
    fname = data_path + "/invtt_" + std::to_string(rank) + ".bin";
    fillArrayFromFile(fname.c_str(), inv_tt.data(), inv_tt.size(),
                      sizeof(inv_tt[0]));

    //____________________________________________________________
    // Compute a product using stbmm

    Block tpltzblocks_Nm1[nb_blocks_loc];
    Tpltz Nm1;

    // flags for Toeplitz product strategy
    Flag flag_stgy;
    flag_stgy_init_auto(&flag_stgy);
    flag_stgy.flag_skip_build_gappy_blocks = 1;

    // Block definition
    defineBlocks_avg(tpltzblocks_Nm1, inv_tt.data(), nb_blocks_loc,
                     local_blocks_sizes.data(), lambda, gif);

    // Matrix definition
    defineTpltz_avg(&Nm1, nb_samp, 1, 1, tpltzblocks_Nm1, nb_blocks_loc,
                    nb_blocks_tot, gif, nb_samp, flag_stgy, MPI_COMM_WORLD);

    // time product operation
    mappraiser::system_stopwatch stopwatch;
    stbmmProd(&Nm1, signal.data());
    auto end = std::chrono::high_resolution_clock::now();
    auto milliseconds =
        stopwatch.elapsed_time<double, std::chrono::milliseconds>();
    std::cout << "stbmmProd time (samples: " << nb_samp
              << ", lambda: " << lambda << ") -> " << milliseconds << " ms\n";

    // write result to disk
    write_vector_to_file(output_path + "/stbmm_result.bin", signal);

    MPI_Finalize();
    return 0;
}
