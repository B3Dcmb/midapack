//
// Created by sbiquard on 09/11/23.
//

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <midapack.h>
#include <mpi.h>
#include <vector>

#define NS 32
#define NP 4

void mat_print(Mat *A, const std::string &name = "") {
    if (!name.empty())
        std::cout << name << " =\n";
    std::cout << "[";
    // loop over rows
    for (int t = 0; t < A->m; t++) {
        std::cout << ((t == 0) ? "[" : " [");
        int tnnz = t * A->nnz;
        int c = 0;
        // loop over columns
        for (int j = 0; j < A->lcount; j++) {
            if (c < A->nnz && j == A->indices[tnnz + c]) {
                std::cout << A->values[tnnz + c];
                c++;
            } else {
                std::cout << 0;
            }
            if (j + 1 < A->lcount)
                std::cout << ", ";
        }
        std::cout << "]";
        if (t + 1 < A->m)
            std::cout << ",\n";
    }
    std::cout << "]\n";
}

template <typename T>
void vec_print(std::vector<T> &v, const std::string &name = "") {
    if (!name.empty())
        std::cout << name << " = ";
    std::cout << "[";
    ulong n = v.size();
    for (ulong i = 0; i < n; i++) {
        std::cout << v[i];
        if (i + 1 < n)
            std::cout << ", ";
    }
    std::cout << "]\n";
}

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
    std::array<int, m> base_indices = {
        -2, -1, -2, 0,  -2, 1,  0,  0, 1, -1, -2, -2, 1,  -1, -1, 0,
        -2, 0,  1,  -2, 0,  -2, -1, 0, 1, -2, -1, 1,  -1, 0,  -2, -2};
    std::vector<int> indices(m * nnz);
    std::vector<double> values(m * nnz);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < nnz; j++) {
            indices[i * nnz + j] = base_indices[i] * nnz + j;
            values[i * nnz + j] = j == 0 ? 1 : 0.5;
        }
    }

    P.flag_ignore_extra = false;
    MatInit(&P, m, nnz, indices.data(), values.data(), 0, MPI_COMM_WORLD);

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
    std::array<int, m> base_indices = {
        -2, -1, -2, 0,  -2, 1,  0,  0, 1, -1, -2, -2, 1,  -1, -1, 0,
        -2, 0,  1,  -2, 0,  -2, -1, 0, 1, -2, -1, 1,  -1, 0,  -2, -2};
    std::vector<int> indices(m * nnz);
    std::vector<double> values(m * nnz);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < nnz; j++) {
            indices[i * nnz + j] = base_indices[i] * nnz + j;
            values[i * nnz + j] = j == 0 ? 1 : 0.5;
        }
    }

    P.flag_ignore_extra = false;
    MatInit(&P, m, nnz, indices.data(), values.data(), 0, MPI_COMM_WORLD);

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

    std::cout << "Testing MatVecProd\n";
    int info = test_MatVecProd(3);
    if (info != 0)
        std::cout << "test MatVecProd failed with status " << info << "\n";
    else
        std::cout << "test MatVecProd successful\n";

    std::cout << "Testing TrMatVecProd\n";
    info = test_TrMatVecProd(3);
    if (info != 0)
        std::cout << "test TrMatVecProd failed with status " << info << "\n";
    else
        std::cout << "test TrMatVecProd successful\n";

    MPI_Finalize();
    return 0;
}
