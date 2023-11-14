//
// Created by sbiquard on 14/11/23.
//

#include <array>
#include <cassert>
#include <iostream>
#include <mappraiser.h>
#include <mpi.h>

#define NS 32
#define NP 4

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

void init_fake_mat(Mat *P, MPI_Comm comm) {
    const int m = NS;

    // 4 pixels (2 extra)
    std::array<int, m> base_indices = {
        -2, -1, -2, 0,  -2, 1,  0,  0, 1, -1, -2, -2, 1,  -1, -1, 0,
        -2, 0,  1,  -2, 0,  -2, -1, 0, 1, -2, -1, 1,  -1, 0,  -2, -2};

    auto indices = static_cast<int *>(malloc(sizeof(int) * m));
    auto values = static_cast<double *>(malloc(sizeof(double) * m));

    for (int i = 0; i < m; i++) {
        indices[i] = base_indices[i];
        values[i] = 1;
    }

    MatInit(P, m, 1, indices, values, 0, comm);

    assert(P->lcount == NP);
}

int my_test(bool verbose = false) {
    //____________________________________________________________
    // MPI initialization

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size != 1) {
        if (rank == 0)
            std::cerr << "this test is meant to be run by 1 process";
        return 1;
    }

    std::cout << "Running my_test" << std::endl;

    //____________________________________________________________
    // Create fake pointing matrix P and gaps G

    Mat P;
    init_fake_mat(&P, MPI_COMM_WORLD);

    mat_print(&P, "pointing matrix");

    Gap G;
    G.ngap = build_pixel_to_time_domain_mapping(&P);

    build_gap_struct(0, &G, &P);

    if (verbose)
        print_gap_info(&G);

    // free memory
    free(P.ll);
    free(P.id_last_pix);
    MatFree(&P);
    free(G.id0gap);
    free(G.lgap);

    return 0;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int info = my_test(true);
    MPI_Barrier(MPI_COMM_WORLD);
    if (info != 0)
        std::cout << "test failed with status " << info << std::endl;
    else
        std::cout << "test successful" << std::endl;

    MPI_Finalize();
    return 0;
}
