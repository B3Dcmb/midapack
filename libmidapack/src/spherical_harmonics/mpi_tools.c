#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
#include <unistd.h>

// #include "s2hat.h"
// #include "midapack.h"
#include "s2hat_tools.h"


int mpi_create_subset(int number_ranks_to_divive, MPI_Comm initcomm, MPI_Comm *subset_comm){
    
    int i;
    int *ranks_not_const;
    MPI_Group global_mpi_group, mpi_subset_group;

    int tag = 0;
    
    ranks_not_const = (int *)malloc(number_ranks_to_divive*sizeof(int));
    for (i=0; i<number_ranks_to_divive; i++){
        ranks_not_const[i] = i;
    }
    const int* ranks_const = ranks_not_const;

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

int all_reduce_to_all_indices_mappraiser(int *indices_pixel_local, int number_pixel_local, int nside, int* all_sky_pixels_observed, int root, MPI_Comm world_comm)
{
    /* Will get the number_pixel_local first pixels of indices_pixel_local to build */
    int i;
    // int *copy_ell_indices;
    int *com_val;

    long int npix = 12*nside*nside;

    com_val=(int *) calloc( npix,sizeof(int)); // Npix or less because of trash_pix ? To check later

    for (i=0; i<number_pixel_local; i++)
    {
        com_val[indices_pixel_local[i]] = 1;
    }

    // s2m(com_val, indices_ones, indices_pixel_local, number_pixel_local); 
    MPI_Reduce(com_val, all_sky_pixels_observed, npix, MPI_INT, MPI_SUM, root, world_comm);	//maximum index

    // free(copy_ell_indices);
    free(com_val);
    return 0;
}