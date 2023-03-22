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




void mpi_broadcast_s2hat_global_struc(S2HAT_GLOBAL_parameters *Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
    /* Use s2hat routines to broadcast s2hat structures */
    if (Local_param_s2hat.gangrank != -1){
        MPI_pixelizationBcast( &(Global_param_s2hat->pixelization_scheme), Local_param_s2hat.gangroot, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
        MPI_scanBcast(Global_param_s2hat->pixelization_scheme, &(Global_param_s2hat->scan_sky_structure_pixel), Local_param_s2hat.gangroot, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->pixpar.par1), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->pixpar.par2), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->nlmax), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->nmmax), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->nside), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    }
}

int distribute_full_sky_map_into_local_maps_S2HAT(double* full_sky_map, double *local_map_s2hat, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat, int nstokes){
    /* Distribute full sky map in ring ordering, with convention [npix, nstokes] in column-wise order among procs, into local maps */
    distribute_map(Global_param_s2hat.pixelization_scheme, 1, 0, nstokes, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size, 
        local_map_s2hat, full_sky_map, Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    // 1 for the number of maps, 0 for the index of the current map

    return 0;
}


int collect_partial_map_from_pixels(double* local_map_s2hat, double *output_submap, int first_pix, int last_pix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat, int nstokes){
    // Collect specific pixels from all local_map_s2hat to form a submap given first and last pixels 
    int submap_size = last_pix - first_pix; // Submapsize given by pixel numbers

    collect_partialmap(Global_param_s2hat.pixelization_scheme, 1, 0, nstokes, first_pix, last_pix, 
        output_submap, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, submap_size, local_map_s2hat, 
        Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    // 1 map given, 0 is the number of the map

    return 0;
}




void free_covariance_matrix(double ** covariance_matrix_NxN, int lmax){
    int ell_value;
    for (ell_value=0; ell_value<lmax+1; ell_value++){
        free(covariance_matrix_NxN[ell_value]);
    }
    free(covariance_matrix_NxN);    
}

void free_s2hat_GLOBAL_parameters_struct(S2HAT_GLOBAL_parameters *Global_param_s2hat){
    destroy_pixelization(Global_param_s2hat->pixelization_scheme);
    destroy_scan(Global_param_s2hat->scan_sky_structure_pixel);
    free(Global_param_s2hat);
}

void free_s2hat_LOCAL_parameters_struct(S2HAT_LOCAL_parameters *Local_param_s2hat){
    free(Local_param_s2hat->mvals);
    free(Local_param_s2hat);
}

void free_s2hat_parameters_struct(S2HAT_parameters *S2HAT_params){

    free_s2hat_GLOBAL_parameters_struct(S2HAT_params->Global_param_s2hat);
    if (S2HAT_params->Local_param_s2hat->gangrank != -1)
        free_s2hat_LOCAL_parameters_struct(S2HAT_params->Local_param_s2hat);

    free(S2HAT_params);
}

