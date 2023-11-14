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

void mat_print(Mat *A, const std::string &name = "") {
    if (!name.empty())
        std::cout << name << " =" << std::endl;
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
            std::cout << "," << std::endl;
    }
    std::cout << "]" << std::endl;
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
    std::cout << "]" << std::endl;
}

void Px(Mat *P, std::vector<double> &x, std::vector<double> &y) {
    MatVecProd(P, x.data(), y.data(), 0);
}

int test_MatVecProd(bool verbose = false) {
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
    // Create pointing matrix P observing 3 pixels (one extra)

    const int m = 32;  // samples
    const int nnz = 3; // non zero entries per row

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

    //____________________________________________________________
    // pointing matrix construction

    Mat P;
    P.flag_ignore_extra = false;
    MatInit(&P, m, nnz, indices.data(), values.data(), 0, MPI_COMM_WORLD);

    if (verbose)
        mat_print(&P, "P");

    //____________________________________________________________
    // sky map

    const int n = 4; // sky pixels
    std::vector<double> sky(n * nnz);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nnz; j++) {
            sky[i * nnz + j] = pow(2, i) * pow(3, j);
        }
    }

    if (verbose)
        vec_print(sky, "m");

    //____________________________________________________________
    // observe the map

    std::vector<double> tod(m);
    Px(&P, sky, tod);

    if (verbose)
        vec_print(tod, "Pm");

    //____________________________________________________________
    // assert correct result

    std::vector<double> result = {7.,  14., 7.,  28., 7.,  56., 28., 28.,
                                  56., 14., 7.,  7.,  56., 14., 14., 28.,
                                  7.,  28., 56., 7.,  28., 7.,  14., 28.,
                                  56., 7.,  14., 56., 14., 28., 7.,  7.};

    assert(result == tod);

    MatFree(&P);
    return 0;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    std::cout << "Testing MatVecProd" << std::endl;
    int info = test_MatVecProd();
    if (info != 0)
        std::cout << "test MatVecProd failed with status " << info << std::endl;
    else
        std::cout << "test MatVecProd successful" << std::endl;

    MPI_Finalize();
    return 0;
}
