//
// Created by sbiquard on 09/11/23.
//

#include <array>
#include <iostream>
#include <midapack.h>
#include <mpi.h>

void mat_print(Mat *A) {
    for (int i = 0; i < A->m; i++) {
        int innz = i * A->nnz;
        for (int j = 0; j < A->lcount; j++) {
            if (j == A->lindices[j])
                std::cout << A->values[innz + j] << " ";
            else
                std::cout << 0.0 << " ";
        }
        std::cout << std::endl;
    }
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
    std::array<int, m> indices = {0, 0, -1, 1, 0, -1, 1, 1};
    std::array<double, m> values{};
    values.fill(1.0);

    std::cout << values[7] << std::endl;

    Mat P;
    P.flag_ignore_extra = false;
    MatInit(&P, m, nnz, indices.data(), values.data(), 0, MPI_COMM_WORLD);

    std::cout << "P.trash_pix = " << P.trash_pix << std::endl;

    mat_print(&P);

    MPI_Finalize();
    return 0;
}
