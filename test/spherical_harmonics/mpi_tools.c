#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
#include <unistd.h>

#include "s2hat.h"
#include "midapack.h"
#include "s2hat_tools.h"


int mpi_create_subset(int number_ranks_to_divive, MPI_Comm initcomm, MPI_Comm *subset_comm){
    
    int i;
    int *ranks_not_const;
    MPI_Group global_mpi_group, mpi_subset_group;

    int tag = 0;
    

    for (i=0; i<number_ranks_to_divive; i++){
        ranks_not_const[i] = i;
    }
    const int ranks_const[number_ranks_to_divive] = ranks_not_const;

    // MPI_Comm_split(initcomm, color, initrank, &s2hat_comm);
    
    // Get the group of the whole communicator
    MPI_Comm_group(initcomm, &global_mpi_group);

    // Construct group containing all ranks under number_ranks_to_divive
    MPI_Comm_group(initcomm, &mpi_subset_group);
    MPI_Group_incl(global_mpi_group, number_ranks_to_divive, ranks_const, &mpi_subset_group);
    
    MPI_Comm_create_group(initcomm, mpi_subset_group, tag, subset_comm);

    MPI_Group_free(&global_mpi_group);
    MPI_Group_free(&mpi_subset_group);

    return 0;
}
