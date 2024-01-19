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
void print_array(T *arr, int n, const std::string &name = "") {
    if (!name.empty())
        std::cout << name << " =\n";
    std::cout << "[";
    for (int i = 0; i < n; i++) {
        std::cout << arr[i];
        if (i + 1 < n)
            std::cout << ", ";
    }
    std::cout << "]\n";
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

    std::cout << "Running test_time_domain_mapping\n";

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
    free(P.id_last_pix);
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
        std::cout << "test failed with status " << info << "\n";
    else
        std::cout << "test successful\n";

    MPI_Finalize();
    return 0;
}
