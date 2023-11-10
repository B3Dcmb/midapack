//
// Created by sbiquard on 09/11/23.
//

#include <array>
#include <iostream>
#include <midapack.h>
#include <mpi.h>
#include <vector>

void mat_print(Mat *A) {
    // loop over rows
    for (int t = 0; t < A->m; t++) {
        int tnnz = t * A->nnz;
        int c = 0;
        // loop over columns
        for (int j = 0; j < A->lcount; j++) {
            if (c < A->nnz && j == A->indices[tnnz + c]) {
                std::cout << A->values[tnnz + c] << " ";
                c++;
            } else {
                std::cout << 0 << " ";
            }
        }
        std::cout << std::endl;
    }
}

void Px(Mat *P, std::vector<double> &x, std::vector<double> &y) {
    MatVecProd(P, x.data(), y.data(), 0);
}

int main(int argc, char *argv[]) {
    //____________________________________________________________
    // MPI initialization

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //____________________________________________________________
    // Create pointing matrix P observing 3 pixels (one extra)

    const int m = 8;   // samples
    const int nnz = 1; // non zero entries per row

    std::array<int, m> base_indices = {0, 0, -1, 1, 0, -1, 1, 1}; // 3 pixels
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

    std::cout << "Pointing matrix:" << std::endl;
    mat_print(&P);

    //____________________________________________________________
    // sky map

    const int n = 3; // sky pixels
    std::vector<double> sky(n * nnz);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nnz; j++) {
            sky[i * nnz + j] = 2 * i + 1;
        }
    }

    std::cout << "Sky map" << std::endl;
    for (auto &s : sky) {
        std::cout << " " << s;
    }
    std::cout << std::endl;

    //____________________________________________________________
    // observe the map

    std::vector<double> tod(m);
    Px(&P, sky, tod);

    std::cout << "d" << std::endl;
    for (auto &d : tod) {
        std::cout << " " << d;
    }
    std::cout << std::endl;

    MatFree(&P);
    MPI_Finalize();
    return 0;
}
