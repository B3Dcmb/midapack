//
// Created by sbiquard on 14/11/23.
//

#include "new_utils.h"

#include <array>
#include <cassert>
#include <iostream>
#include <mappraiser.h>
#include <mpi.h>

#define NS 32
#define NP 4

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

    MatInit(P, m, 1, indices, values, nullptr, comm);
    assert(P->lcount == NP);
}

int test_time_domain_mapping(bool verbose = false) {
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

    std::cout << "Running test_time_domain_mapping" << std::endl;

    //____________________________________________________________
    // Create fake pointing matrix P and gaps G

    Mat P;
    init_fake_mat(&P, MPI_COMM_WORLD);
#if 0
    mat_print(&P, "pointing matrix");
#endif

    Gap G;
    G.ngap = build_pixel_to_time_domain_mapping(&P);

    build_gap_struct(0, &G, &P);

    if (verbose)
        print_gap_info(&G);

    //____________________________________________________________
    // assert that the results are correct

    const int true_ngap = 10;
    int true_lgap[true_ngap] = {3, 1, 3, 2, 1, 1, 2, 2, 1, 2};
    int64_t true_id0gap[true_ngap] = {0, 4, 9, 13, 16, 19, 21, 25, 28, 30};

    assert(G.ngap == true_ngap);
    for (int i = 0; i < true_ngap; i++) {
        assert(G.lgap[i] == true_lgap[i]);
        assert(G.id0gap[i] == true_id0gap[i]);
    }

    // free memory
    free(P.ll);
    free(P.pix_to_last_samp);
    MatFree(&P);
    free(G.id0gap);
    free(G.lgap);

    return 0;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int info = test_time_domain_mapping();
    MPI_Barrier(MPI_COMM_WORLD);
    if (info != 0)
        std::cout << "test failed with status " << info << std::endl;
    else
        std::cout << "test successful" << std::endl;

    MPI_Finalize();
    return 0;
}
