#ifdef W_MPI
 #include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mapmat.h"
#include "s2hat_tools.h"

/* Taken from greedyreduce in mapmat.c, case ALLREDUCE
    Reduce all pixel sky maps distributed among different procs to a single full sky map
*/
int all_reduce_to_single_map_mappraiser(Mat *A, double* x, int nside, double* out_val){
    int i;
    double *lvalues;
    lvalues = (double *) malloc((A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));	//allocate and set to 0.0 local values
    memcpy(lvalues, x, (A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));		//copy local values into result values
    double *com_val;


    long int npix = 12*nside*nside;

    // com_val=(double *) malloc( A->com_count *sizeof(double));
    // out_val=(double *) malloc( A->com_count *sizeof(double));
    com_val=(double *) malloc( 3*npix *sizeof(double));
    // out_val=(double *) malloc( 3*npix *sizeof(double));
    for(i=0; i < 3*npix; i++){
        com_val[i]=0.0;
        out_val[i]=0.0;
    }
    // s2m(com_val, lvalues, A->com_indices, A->lcount-(A->nnz)*(A->trash_pix));
    s2m(com_val, lvalues, A->lindices, A->lcount-(A->nnz)*(A->trash_pix)); 
    
    /*for(i=0; i < A->com_count; i++){
        printf("%lf ", com_val[i]);
    } */
    MPI_Allreduce(com_val, out_val, 3*npix, MPI_DOUBLE, MPI_SUM, A->comm);	//maximum index
    free(lvalues);

    ///// TO MODIFY : How are the 3 stokes parameters for the map included ?
}

int distribute_map_S2HAT_ordering(double* full_sky_map, double *local_map_s2hat, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){

    distribute_map(Global_param_s2hat.pixelization_scheme, 1, 0, 3, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size, 
        local_map_s2hat, full_sky_map, Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gang_comm);
    // 1 for the number of maps, 0 for the index of the current map, 3 for the number of Stokes parameters

    return 0;
}


int brute_force_transfer_local_maps(Mat *A, double* local_pixel_map_MAPPRAISER, double *local_pixel_map_s2hat, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){

    int nside = Global_param_s2hat.nside;
    double* full_sky_map;


    full_sky_map=(double *) malloc( 3*npix *sizeof(double));
    all_reduce_to_single_map_mappraiser(Mat *A, local_pixel_map_MAPPRAISER, nside, full_sky_map);

    // Change ordering of full_sky_map to have what S2HAT expects

    distribute_map_S2HAT_ordering(full_sky_map, local_pixel_map_s2hat, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);
}

