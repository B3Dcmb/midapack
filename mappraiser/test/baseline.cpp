//
// Created by sbiquard on 7/7/23.
//

#include "mappraiser/RunningSum.h"
#include "mappraiser/rng.h"
#include "mappraiser/stopwatch.h"

#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

// template<typename T>
// void printVector(const std::vector<T> &vec) {
//     for (const auto &element: vec) { std::cout << element << ", "; }
//     std::cout << std::endl;
// }

// template<typename T>
// void printVector(const std::vector<T> &vec, size_t n_to_print, size_t offset = 0) {
//     if (offset + n_to_print > vec.size()) {
//         offset     = 0;
//         n_to_print = vec.size();
//     }
//     for (size_t i = offset; i < offset + n_to_print; ++i) { std::cout << vec[i] << ", "; }
//     std::cout << std::endl;
// }

template<typename T>
void flagged_running_average_naive(int samples, const T *values, const uint8_t *valid, double *baseline,
                                   int bandwidth) {
    // half size of window
    int w0 = bandwidth / 2;

#pragma omp parallel for default(shared) schedule(static)
    for (int i = 0; i < samples; ++i) {
        double avg   = 0;
        int    count = 0;
        int    start = i - w0;
        int    end   = i + w0 - 1;

        // compute the moving average of (valid) samples over [i-w0, i+w0]
#pragma omp simd
        for (int j = start; j < end + 1; ++j) {
            if (0 <= j && j < samples && valid[j]) {
                avg += values[j];
                ++count;
            }
        }

        baseline[i] = avg / static_cast<double>(count);
    }
}

template<typename T>
void flagged_running_average_smart(int samples, const T *buf, const uint8_t *valid, double *baseline, int bandwidth) {
    RunningSum<T, uint8_t> sum(bandwidth);
    for (int i = 0; i < samples; ++i) {
        sum.process_index(samples, i, buf, valid);
        baseline[i] = sum.baseline_value();
    }
}

bool allclose(int n, const double *t1, const double *t2, double rtol = 1e-5, double atol = 1e-8) {
    for (int i = 0; i < n; ++i) {
        double a = t1[i];
        double b = t2[i];
        if (std::abs(a - b) > (atol + rtol * std::abs(b))) return false;
    }
    return true;
}

int main(int argc, char *argv[]) {
    //____________________________________________________________
    // Parse command line arguments

    if (argc < 3) {
        std::cerr << "Usage: ./baseline <nval> <lambda>" << std::endl;
        return (EXIT_FAILURE);
    }

    // Number of values
    auto nval = std::stoi(argv[1]);

    // Size of window
    auto lambda = std::stoi(argv[2]);

    //____________________________________________________________
    // Generate random data

    // Uniform randoms in [0,1]
    using ValType = double;
    std::vector<ValType> values(nval);
    mappraiser::rng_dist_uniform_01(nval, 0, 0, 0, 0, values.data());

    for (int i = 0; i < nval; ++i) { values[i] *= i; }

    // printVector(values);

    // Valid intervals
    std::vector<uint8_t> valid(nval, 1);

    int offset  = nval / 20;
    int lgap    = nval / 10;
    int spacing = nval / 2;
    for (int i = offset; i < nval - lgap; i += spacing) {
        for (int j = 0; j < lgap; ++j) { valid[i + j] = 0; }
    }

    // printVector(valid);

    //____________________________________________________________
    // Baseline computation

    std::vector<double> baseline(nval);

    mappraiser::system_stopwatch watch;
    flagged_running_average_naive<ValType>(nval, values.data(), valid.data(), baseline.data(), lambda);
    auto time_naive = watch.elapsed_time<double, std::chrono::milliseconds>();

    // printVector(baseline);

    //____________________________________________________________
    // Baseline computation (new version)

    std::vector<double> baseline_smart(nval);

    mappraiser::system_stopwatch watch_smart;
    flagged_running_average_smart<ValType>(nval, values.data(), valid.data(), baseline_smart.data(), lambda);
    auto time_smart = watch_smart.elapsed_time<double, std::chrono::milliseconds>();

    // printVector(baseline_smart);

    //____________________________________________________________
    // Results (timing + correctness)

    std::cout << "threads     : " << omp_get_max_threads() << std::endl
              << "nval        : " << nval << std::endl
              << "lambda      : " << lambda << std::endl
              << "time (naive): " << time_naive << " ms" << std::endl
              << "time (smart): " << time_smart << " ms" << std::endl;

    if (allclose(nval, baseline.data(), baseline_smart.data())) {
        std::cout << "Baselines are close (within relative distance of 1e-5)" << std::endl;
    } else {
        std::cout << "Baselines are not close" << std::endl;
    }

    return 0;
}