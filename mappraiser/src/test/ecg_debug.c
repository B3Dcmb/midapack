//
// Created by sbiquard on 5/17/23.
//

#include <mappraiser.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    // MPI initialization
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Load data


    // Call MLmap
    // MLmap(MPI_COMM_WORLD, );

    printf("hello from rank %d\n", rank);

    // Finalize

    MPI_Finalize();
    return 0;
}
