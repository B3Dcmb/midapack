//
// Created by sbiquard on 09/11/23.
//

#include "new_utils.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <vector>

#define NS 32
#define NP 4

void Px(Mat *P, std::vector<double> &x, std::vector<double> &y) {
    MatVecProd(P, x.data(), y.data());
}

int test_MatVecProd(const int nnz = 3, bool verbose = false) {
    //____________________________________________________________
    // MPI initialization

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size > 1) {
        std::cout << "this test is not meant to be run by several processes";
        return 1;
    }

    //____________________________________________________________
    // Create pointing matrix P

    Mat P;

    const int m = NS;

    // 4 pixels (2 extra)
    std::vector<int> indices = {-2, -1, -2, 0,  -2, 1,  0,  0, 1,  -1, -2,
                                -2, 1,  -1, -1, 0,  -2, 0,  1, -2, 0,  -2,
                                -1, 0,  1,  -2, -1, 1,  -1, 0, -2, -2};
    std::vector<double> values(m * nnz);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < nnz; j++) {
            values[i * nnz + j] = j == 0 ? 1 : 0.5;
        }
    }

    P.ignore_extra = false;
    MatInit(&P, m, nnz, indices.data(), values.data(), nullptr, MPI_COMM_WORLD);

    if (verbose)
        mat_print(&P, "P");

    //____________________________________________________________
    // sky map

    const int n = NP;

    std::vector<double> s(n * nnz);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nnz; j++) {
            s[i * nnz + j] = pow(2, i) * pow(3, j);
        }
    }

    if (verbose)
        vec_print(s, "s");

    //____________________________________________________________
    // observe the map

    std::vector<double> Ps(m);
    Px(&P, s, Ps);

    if (verbose)
        vec_print(Ps, "Ps");

    //____________________________________________________________
    // assert correct result

    std::vector<double> result = {7.,  14., 7.,  28., 7.,  56., 28., 28.,
                                  56., 14., 7.,  7.,  56., 14., 14., 28.,
                                  7.,  28., 56., 7.,  28., 7.,  14., 28.,
                                  56., 7.,  14., 56., 14., 28., 7.,  7.};

    assert(result == Ps);

    MatFree(&P);
    return 0;
}

void Ptx(Mat *P, std::vector<double> &y, std::vector<double> &x) {
    TrMatVecProd(P, y.data(), x.data());
}

int test_TrMatVecProd(const int nnz = 3, bool verbose = false) {
    //____________________________________________________________
    // MPI initialization

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size > 1) {
        std::cout << "this test is not meant to be run by several processes";
        return 1;
    }

    //____________________________________________________________
    // Create pointing matrix P

    Mat P;

    const int m = NS;

    // 4 pixels (2 extra)
    std::vector<int> indices = {-2, -1, -2, 0,  -2, 1,  0,  0, 1,  -1, -2,
                                -2, 1,  -1, -1, 0,  -2, 0,  1, -2, 0,  -2,
                                -1, 0,  1,  -2, -1, 1,  -1, 0, -2, -2};
    std::vector<double> values(m * nnz);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < nnz; j++) {
            values[i * nnz + j] = j == 0 ? 1 : 0.5;
        }
    }

    P.ignore_extra = false;
    MatInit(&P, m, nnz, indices.data(), values.data(), nullptr, MPI_COMM_WORLD);

    if (verbose)
        mat_print(&P, "P");

    //____________________________________________________________
    // time stream

    std::vector<double> d = {7.,  14., 7.,  28., 7.,  56., 28., 28.,
                             56., 14., 7.,  7.,  56., 14., 14., 28.,
                             7.,  28., 56., 7.,  28., 7.,  14., 28.,
                             56., 7.,  14., 56., 14., 28., 7.,  7.};

    if (verbose)
        vec_print(d, "d");

    //____________________________________________________________
    // project on the sky

    const int n = NP;
    std::vector<double> Ptd(n * nnz);

    Ptx(&P, d, Ptd);

    if (verbose)
        vec_print(Ptd, "P^t d");

    //____________________________________________________________
    // assert correct result

    std::vector<double> result = {77.,  38.5, 38.5, 98.,  49.,  49.,
                                  224., 112., 112., 336., 168., 168.};
    assert(result == Ptd);

    MatFree(&P);
    return 0;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    std::cout << "Testing MatVecProd" << std::endl;
    int info = test_MatVecProd(3);
    if (info != 0)
        std::cout << "test MatVecProd failed with status " << info << std::endl;
    else
        std::cout << "test MatVecProd successful" << std::endl;

    std::cout << "Testing TrMatVecProd" << std::endl;
    info = test_TrMatVecProd(3);
    if (info != 0)
        std::cout << "test TrMatVecProd failed with status " << info
                  << std::endl;
    else
        std::cout << "test TrMatVecProd successful" << std::endl;

    MPI_Finalize();
    return 0;
}
