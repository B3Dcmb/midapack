#ifdef W_MPI
 #include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <chealpix.h>
#include "s2hat.h"
#include "midapack.h"
#include "s2hat_tools.h"
#include "domain_generalization.h"


int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix, s2hat_dcomplex *local_alm, int domain_PCG_computation, int bool_apply_filter, int nstokes, S2HAT_parameters *S2HAT_params)
{   // Initalize PCG_var structure with S2HAT parameters and flags, before attribution of the data

    PCG_variable->local_map_pix = local_map_pix; // Local map in pixel domain, with MAPPRAISER convention
    PCG_variable->local_alm = local_alm; // Local alm in harmonic domain, with S2HAT convention

    if ((local_alm == NULL) && ((domain_PCG_computation == 1) || (bool_apply_filter == 1)))
    {
        PCG_variable->local_alm = (s2hat_dcomplex *) calloc( nstokes * S2HAT_params->size_alm, sizeof(s2hat_dcomplex));
    }

    // Specifications of the computation
    PCG_variable->domain_PCG_computation = domain_PCG_computation; // Domain chosen for PCG computation : 0 pixel, 1 harmonic
    PCG_variable->bool_apply_filter = bool_apply_filter; // 0 no filter applied, 1 Wiener filter applied
    PCG_variable->nstokes = nstokes; // Number of Stokes parameters for the computation, 3 usually means T, Q, U

    // S2HAT structures necessary for harmonic domain
    PCG_variable->S2HAT_parameters = S2HAT_params;
    // Contains S2HAT_parameters->Global_param_s2hat and S2HAT_parameters->Local_param_s2hat

    // Initialize update flags to 0 "no update needed"
    PCG_variable->does_map_pixel_need_update = 0; // 0 no update needed, 1 need update from local_alm
    PCG_variable->does_local_alm_need_update = 0; // 0 no update needed, 1 need update from local_map_pix
    return 0;
}


int init_butterfly_communication(Butterfly_struct *Butterfly_struc, Mat *A, int *indices_in, int count_in, int *indices_out, int count_out, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    Butterfly_struc->steps = log_2(size);
    Butterfly_struc->S = (int **)malloc(Butterfly_struc->steps * sizeof(int *)); // allocate sending maps tab
    Butterfly_struc->R = (int **)malloc(Butterfly_struc->steps * sizeof(int *)); // allocate receiving maps tab
    Butterfly_struc->nS = (int *)malloc(Butterfly_struc->steps * sizeof(int));   // allocate sending map sizes tab
    Butterfly_struc->nR = (int *)malloc(Butterfly_struc->steps * sizeof(int));   // allocate receiving map size tab
    butterfly_init(A->lindices + (A->nnz) * (A->trash_pix), A->lcount - (A->nnz) * (A->trash_pix), Butterfly_struc->R, Butterfly_struc->nR, Butterfly_struc->S, Butterfly_struc->nS, &(A->com_indices), &(A->com_count), Butterfly_struc->steps, comm);
    butterfly_init(int *indices, int count, int **R, int *nR, int **S, int *nS, int **com_indices, int *com_count, int steps, MPI_Comm comm)

    butterfly_reshuffle_init(indices_in, count_in, indices_out, count_out, Butterfly_struc->R, Butterfly_struc->nR, Butterfly_struc->S, Butterfly_struc->nS, int **com_indices, int *com_count, Butterfly_struc->steps, comm);




}



/* Taken from greedyreduce in mapmat.c, case ALLREDUCE
    Reduce all pixel sky maps distributed among different procs to a single full sky map
*/
// int all_reduce_to_single_map_mappraiser(Mat *A, PCG_var *PCG_variable, int nside, double* output_fullsky_map, int root){
int all_reduce_to_single_map_mappraiser(Mat *A, double* local_map_pix, int nside, double* output_fullsky_map, int root){
    int i;
    double *lvalues;
    lvalues = (double *) malloc((A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));	//allocate and set to 0.0 local values
    // memcpy(lvalues, PCG_variable->local_map_pix, (A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));		//copy local values into result values
    memcpy(lvalues, local_map_pix, (A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));		//copy local values into result values
    double *com_val;

    long int npix = 12*nside*nside;

    com_val=(double *) calloc( 3*npix,sizeof(double)); // Npix or less because of trash_pix ? To check later

    s2m(com_val, lvalues, &(A->lindices[(A->nnz)*(A->trash_pix)]), A->lcount-(A->nnz)*(A->trash_pix)); 
    
    /*for(i=0; i < A->com_count; i++){
        printf("%lf ", com_val[i]);
    } */
    // MPI_Allreduce(com_val, output_fullsky_map, 3*npix, MPI_DOUBLE, MPI_SUM, A->comm);	//maximum index
    MPI_Reduce(com_val, output_fullsky_map, 3*npix, MPI_DOUBLE, MPI_SUM, root, A->comm);	//maximum index

    free(lvalues);
    free(com_val);

    // Here, the output_fullsky_map is in nest ordering, ordering as [nstokes, npix] in column-wise ordering
}



int brute_force_transfer_local_maps(Mat *A, double* local_pixel_map_MAPPRAISER, double *local_pixel_map_s2hat, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
    // Temporary transfer method, which compile every thing and redistribute everything

    int nside = Global_param_s2hat.nside;
    double* full_sky_map;

    int npix = 12*nside*nside;

    double *full_sky_map_ring = (double *) malloc( 3*npix *sizeof(double));

    if (Local_param_s2hat.gangroot == Local_param_s2hat.gangrank){
        double *full_sky_map_nest = (double *) malloc( 3*npix *sizeof(double));
        all_reduce_to_single_map_mappraiser(A, local_pixel_map_MAPPRAISER, nside, full_sky_map_nest, Local_param_s2hat.gangroot);

        
        compute_full_map_nest2ring(full_sky_map_nest, full_sky_map_ring, nside, A->nnz, npix);
        // Change ordering of full_sky_map to have what S2HAT expects, from MAPPRAISER nest [nstokes, npix] to S2HAT ring [nstokes, npix]
        free(full_sky_map_nest);
    }

    distribute_full_sky_map_into_local_maps_S2HAT(full_sky_map_ring, local_pixel_map_s2hat, Global_param_s2hat, Local_param_s2hat, A->nnz);
    // Distribute global full sky map into local maps in each processors
    free(full_sky_map_ring);
}


int global_map_2_harmonic(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, S2HAT_parameters *S2HAT_params){
    // Transform local_maps pixel distribution from MAPPRAISER into a harmonic S2HAT a_lm distribution

    S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters Local_param_s2hat = *(S2HAT_params->Local_param_s2hat);

    double *local_pixel_map_s2hat = (double *)malloc(Local_param_s2hat.map_size*sizeof(double));
    brute_force_transfer_local_maps(A, local_pixel_map_MAPPRAISER, local_pixel_map_s2hat, Global_param_s2hat, Local_param_s2hat);
    // Temporary transfer method

    apply_pix2alm(local_pixel_map_s2hat, local_alm_s2hat, A->nnz, S2HAT_params);
    free(local_alm_s2hat);
}




/* Maybe useless ? */
// int obtain_pixel_distrib_from_S2HAT(Mat *A, int nside, int *lindices_ring_ordering, int flag_to_sort_output)
// {   // 
//     int index, index_pixel_ring;
//     int size_list_indices = (A->lcount-(A->nnz)*(A->trash_pix))/(A->nnz);
//     for(index=0; index < size_list_indices; index++){
//         nest2ring( nside, A->lindices[(A->nnz)*(A->trash_pix)+index*(A->nnz)], &(lindices_ring_ordering[index]));
//         // Change pixel in nest ordering to ring ordering
//     }
//     if (flag_to_sort_output == 1){
//         int sflag = 3; // Choose to have counting_sort method to sort, which also allows to delete redundant elements
//         ssort(lindices_ring_ordering, size_list_indices, sflag);
//     }
// }

int global_harmonic_2_map(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A,  S2HAT_parameters *S2HAT_params){
    S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters Local_param_s2hat = *(S2HAT_params->Local_param_s2hat);
   return 0;
}



int update_PCG_var(PCG_var *PCG_variable, Mat *A)
{   // Update either local_map_pix from local_alm or the inverse, or do nothing, depending on the update flags
    // Reset the update flags to 0 afterwards
    if (PCG_variable->does_map_pixel_need_update == 1){
        global_harmonic_2_map(PCG_variable->local_map_pix, PCG_variable->local_alm, A, PCG_variable->S2HAT_parameters);
        PCG_variable->does_map_pixel_need_update = 0;
    }

    if ((PCG_variable->does_local_alm_need_update == 1) && (PCG_variable->local_alm != NULL)){
        global_map_2_harmonic(PCG_variable->local_map_pix, PCG_variable->local_alm, A, PCG_variable->S2HAT_parameters);
        PCG_variable->does_local_alm_need_update = 0;
    }
}


/* Old - not usable : idea happen to not be applicable */
// int global_harmonic_2_map(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat)
// {
//     int nside = Global_param_s2hat.nside;

//     double *local_map_pix = (double *) malloc(Local_param_s2hat.map_size*sizeof(double));

//     apply_alm2pix(local_alm_s2hat, local_map_pix, Global_param_s2hat, Local_param_s2hat);

//     int *lindices_ring_ordering = (int *)malloc(A->lcount-(A->nnz)*(A->trash_pix)*sizeof(int));

//     obtain_pixel_distrib_from_S2HAT(A, nside, lindices_ring_ordering, 1);

//     lindices_ring_ordering[0]

//     int index_proc;
//     for(index_proc=0; index_proc<Local_param_s2hat.gangsize; index_proc++){
//         Local_param_s2hat.gangroot = index_proc;

        
//         collect_partial_map_from_pixels(local_map_pix, double *output_submap, int first_pix, int last_pix, Global_param_s2hat, Local_param_s2hat, nstokes);

//         Local_param_s2hat.gangroot = 0;

//     }
// }