//
// Created by sbiquard on 14/11/23.
//

#include <iostream>
#include <mappraiser.h>
#include <mpi.h>

int my_test() { return 0; }

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    std::cout << "Testing my_test" << std::endl;
    int info = my_test();
    if (info != 0)
        std::cout << "test failed with status " << info << std::endl;
    else
        std::cout << "test successful" << std::endl;

    MPI_Finalize();
    return 0;
}
