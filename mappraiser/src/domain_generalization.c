#ifdef W_MPI
 #include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
// #include <chealpix.h>
// #include "s2hat.h"
// #include "mappraiser.h"
// #include "s2hat_tools.h"
#include "domain_generalization.h"


int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix)
{   // Initalize PCG_var structure with S2HAT parameters and flags, before attribution of the data

    PCG_variable->local_map_pix = local_map_pix; // Local map in pixel domain, with MAPPRAISER convention

    // Specifications of the computation
    // PCG_variable->domain_PCG_computation = domain_PCG_computation; // Domain chosen for PCG computation : 0 pixel, 1 harmonic
    // PCG_variable->bool_apply_filter = bool_apply_filter; // 0 no filter applied, 1 Wiener filter applied
    // PCG_variable->nstokes = nstokes; // Number of Stokes parameters for the computation, 3 usually means T, Q, U

    // S2HAT structures necessary for harmonic domain
    // PCG_variable->S2HAT_parameters = S2HAT_params;
    // Contains S2HAT_parameters->Global_param_s2hat and S2HAT_parameters->Local_param_s2hat

    // Initialize update flags to 0 "no update needed"
    // PCG_variable->does_map_pixel_need_update = 0; // 0 no update needed, 1 need update from local_alm
    // PCG_variable->does_local_alm_need_update = 0; // 0 no update needed, 1 need update from local_map_pix
    return 0;
}


int init_harmonic_superstruct(Mat *A, Harmonic_superstruct *Harm_struct, double *mask_binary, int nside, int lmax, char *c_ell_path, int number_correlations, int iter_alm, float error_alm)
{
    /* Initalize all structures necessary for harmonic structures : S2HAT_params for S2HAT operations, and the 
       Relies on the fact that nstokes == A->nnz
    */

    int flag_reshuffle_butterfly = 1; // Put Butterfly in reshuffle mode

    int rank, size;
    // int number_pixels_local, number_pixel_total;
    // int *pixel_numbered_ring;
    MPI_Comm worldcomm = A->comm;
    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &size);

    Files_path_WIENER_FILTER *Files_WF_struct = &(Harm_struct->S2HAT_params.Files_WF_struct);
    // printf("%d ~~~ Initializing struct_WF \n", rank); fflush(stdout);
    init_files_struct_WF(Files_WF_struct, nside, lmax, c_ell_path, number_correlations);
    // printf("%d ~~~ Initializing superstruct S2HAT_params \n", rank); fflush(stdout);
    init_s2hat_parameters_superstruct(Files_WF_struct, mask_binary, A->nnz, iter_alm, error_alm, &(Harm_struct->S2HAT_params), worldcomm);
    // Initialize S2HAT structures

    S2HAT_parameters *S2HAT_params = &(Harm_struct->S2HAT_params);
    int nstokes = S2HAT_params->nstokes;

    int *indices_local_MAPPRAISER_nest = A->lindices + (A->nnz) * (A->trash_pix); // Indices local MAPPRAISER in nest distribution
    int number_pixels_MAPPRAISER = A->lcount - (A->nnz) * (A->trash_pix);

    int *ordered_indices_local_MAPPRAISER_ring = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));
    int *projector_ring2nest = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));
    int *projector_nest2ring = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));

    // printf("%d ~~~ Projecting indices from nest2ring \n", rank); fflush(stdout);
    // get_projectors_indices(indices_local_MAPPRAISER_nest, ordered_indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, nstokes, S2HAT_params->Global_param_s2hat.nside, projector_ring2nest, projector_nest2ring, rank);
    get_projectors_indices(indices_local_MAPPRAISER_nest, ordered_indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, nstokes, S2HAT_params->Global_param_s2hat.nside, projector_ring2nest, projector_nest2ring);
    // As we have nest distribution with MAPPRAISER convention, convert the indices into ring distribution with S2HAT convention

    S2HAT_params->local_projector_values_ring2nest = projector_ring2nest;
    S2HAT_params->local_projector_values_nest2ring = projector_nest2ring;
    
    int *indices_local_S2HAT_ring = S2HAT_params->Local_param_s2hat.pixel_numbered_ring; // Indices local S2HAT in ring distribution
    int number_pixels_S2HAT = nstokes * S2HAT_params->Local_param_s2hat.map_size; // Local map size

    // printf("%d ~~~ Preparing butterfly communication part 1 !!! \n", rank); fflush(stdout);
    prepare_butterfly_communication(ordered_indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, indices_local_S2HAT_ring, number_pixels_S2HAT, flag_reshuffle_butterfly, &(Harm_struct->MAPPRAISER_to_S2HAT_butterfly), worldcomm);
    // printf("%d ~~~ Preparing butterfly communication part 2 !!! \n", rank); fflush(stdout);
    prepare_butterfly_communication(indices_local_S2HAT_ring, number_pixels_S2HAT, ordered_indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, flag_reshuffle_butterfly, &(Harm_struct->S2HAT_to_MAPPRAISER_butterfly), worldcomm);

    /* Attributing the Harmonic superstruc variables*/
    // Harmonic_superstruct->S2HAT_parameters = S2HAT_params;
    // Attributing the S2HAT superstructure to Harmonic superstructure

    // Harmonic_superstruct->S2HAT_to_MAPPRAISER = ring2MAPP_butterfly;
    // Harmonic_superstruct->MAPPRAISER_to_S2HAT = MAPP2ring_butterfly;
    // Attributing the Butterfly structures to Harmonic superstructure
    // printf("%d ~~~ Done with initializing Harmonic struct !!! \n", rank); fflush(stdout);
    return 0;
}


// int get_projectors_indices_v2(int *indices_nest, int *ordered_indices_ring, int size_indices, int nstokes, int nside, int *projector_ring2nest, int *projector_nest2ring, int rank)
// {
//   /* Build projectors to convert indices_nest from the nest pixel distribution to ring pixel distribution,
//      then reorder them so that the indices in ring order are monotonous 
//      The indices_nest are expected to be not have any redundancy */

//   int i, j;
//   // printf("%d <<< Starting projecting indices_v2 !!! \n", rank); fflush(stdout);
//   int *indices_ring = (int *)malloc(size_indices*sizeof(int));

//   // printf("%d <<< Converting nest2ring !!! \n", rank); fflush(stdout);
//   convert_indices_nest2ring(indices_nest, indices_ring, size_indices, nstokes, nside);

//   // printf("%d <<< memcpy rings !!! \n", rank); fflush(stdout);
//   memcpy(ordered_indices_ring, indices_ring, size_indices*sizeof(int));

//   // printf("%d <<< ssort !!! \n", rank); fflush(stdout);
// //   ssort(ordered_indices_ring, size_indices, 0); 
//   for (i=0; i<size_indices; i++) 
//   {
//     projector_nest2ring[i] = i;
//   }
//   ssort_with_indices(ordered_indices_ring, projector_nest2ring, size_indices, 0);
//   // Argument flag=0 to use quicksort to sort the indices

//   // printf("%d <<< Finally for loops !!! \n", rank); fflush(stdout);
// //   for (i=0; i<size_indices; i++) 
// //   {
// //     for (j=0; j<size_indices; j++)
// //     {
// //       if (ordered_indices_ring[i] == indices_ring[j])
// //       {
// //         projector_ring2nest[i] = j;
// //         projector_nest2ring[j] = i;
// //         break ;
// //         // The indices_nest are expected to be not have any redundancy
// //       }
// //     }
// //   }
//   int index;
//   for (i=0; i<size_indices; i++) 
//   {
//     // index = elem_in_list_elem(ordered_indices_ring[i], indices_ring, size_indices);
//     // projector_ring2nest[i] = index;
//     projector_ring2nest[projector_nest2ring[i]] = i;
//   }
//   // printf("%d <<< Done projecting indices_v2 !!! \n", rank); fflush(stdout);
//   free(indices_ring);
//   // free(ordered_indices_ring);
//   return 0;
// }


int global_map_2_harmonic(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, Harmonic_superstruct *Harmonic_sup){
//     // Transform local_maps pixel distribution from MAPPRAISER into a harmonic S2HAT a_lm distribution

    // S2HAT_GLOBAL_parameters *Global_param_s2hat = &(Harmonic_sup->S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(Harmonic_sup->S2HAT_params.Local_param_s2hat);
    // Butterfly_superstruct *Butterfly_map2harmonic = &(Harmonic_sup->MAPPRAISER_to_S2HAT_butterfly);
    
    MPI_Comm worldcomm = A->comm;

    int rank, size;
    MPI_Comm_size(worldcomm, &size);
    MPI_Comm_rank(worldcomm, &rank);

    int nstokes = Harmonic_sup->S2HAT_params.nstokes;
    double *local_pixel_map_S2HAT;

    int number_pixels_MAPPRAISER = A->lcount-(A->nnz)*(A->trash_pix);
    double *local_pixel_map_MAPPRAISER_ring = (double *)calloc(number_pixels_MAPPRAISER,sizeof(double));
    int *local_pixel_indices_MAPPRAISER_ring = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));
    // double *local_pixel_map_MAPPRAISER_nest = (double *)calloc(true_size_local,sizeof(double));

    int max_size_test = 10;
    int index=0;
    int check_monotony = monotony(A->lindices + (A->nnz) * (A->trash_pix), number_pixels_MAPPRAISER);
    printf("%d »»»»»» Initial MAPP nest map with size %d monotony_check %d : From pixel 0 to %d - %f %d -", rank, number_pixels_MAPPRAISER, check_monotony, max_size_test, local_pixel_map_MAPPRAISER[index], A->lindices[index]);
    for (index=1;index<max_size_test;index++){
        printf("- %f %d -", local_pixel_map_MAPPRAISER[index], A->lindices[index]);
    }
    printf("\n"); fflush(stdout);

    project_values_into_different_scheme(local_pixel_map_MAPPRAISER, number_pixels_MAPPRAISER, Harmonic_sup->S2HAT_params.local_projector_values_nest2ring, local_pixel_map_MAPPRAISER_ring);
    int *local_pixel_indices_MAPPRAISER_ring_tmp = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));
    convert_indices_nest2ring(A->lindices + (A->nnz) * (A->trash_pix), local_pixel_indices_MAPPRAISER_ring_tmp, number_pixels_MAPPRAISER, nstokes, Harmonic_sup->S2HAT_params.Global_param_s2hat.nside);
    project_int_values_into_different_scheme(local_pixel_indices_MAPPRAISER_ring_tmp, number_pixels_MAPPRAISER, Harmonic_sup->S2HAT_params.local_projector_values_nest2ring, local_pixel_indices_MAPPRAISER_ring);
    free(local_pixel_indices_MAPPRAISER_ring_tmp);

    index = 0;
    check_monotony = monotony(local_pixel_indices_MAPPRAISER_ring, number_pixels_MAPPRAISER);
    printf("%d »»»»»» Transform MAPP local map into ring map with size %d monotony_check %d : From pixel 0 to %d - %f %d -", rank, number_pixels_MAPPRAISER, check_monotony, max_size_test, local_pixel_map_MAPPRAISER_ring[index], local_pixel_indices_MAPPRAISER_ring[index]);
    for (index=1;index<max_size_test;index++){
        printf("- %f %d -", local_pixel_map_MAPPRAISER_ring[index], local_pixel_indices_MAPPRAISER_ring[index]);
    }
    printf("\n"); fflush(stdout);

    int number_pixels_S2HAT = nstokes * Harmonic_sup->S2HAT_params.Local_param_s2hat.map_size; // Local map size
    if (number_pixels_S2HAT)
        local_pixel_map_S2HAT = (double *)calloc(number_pixels_S2HAT,sizeof(double));

    check_monotony = monotony(local_pixel_indices_MAPPRAISER_ring, number_pixels_MAPPRAISER);
    printf("%d »»»»»» ........ Performing butterfly communication ! monotony_check %d number_pixels_S2HAT %d \n",rank, check_monotony, number_pixels_S2HAT); fflush(stdout);
    // int number_of_nan = 0;
    // for (index=0;index<A->lcount-(A->nnz)*(A->trash_pix);index++){
    //     if (!(local_pixel_map_MAPPRAISER_ring[index] == local_pixel_map_MAPPRAISER_ring[index])){
    //         number_of_nan++;
    //         }
    //     if (!(local_pixel_indices_MAPPRAISER_ring[index] == local_pixel_indices_MAPPRAISER_ring[index])){
    //         number_of_nan++;
    //         }
    // }
    // if (number_of_nan)
    //     printf("%d «««««««« NUMBER OF NANS NON-ZERO FOR local_pixel_map_MAPPRAISER_ring before butterfly comm : %d \n", rank, number_of_nan); fflush(stdout);
    perform_butterfly_communication(local_pixel_map_MAPPRAISER_ring, local_pixel_indices_MAPPRAISER_ring, A->lcount-(A->nnz)*(A->trash_pix), local_pixel_map_S2HAT, Harmonic_sup->S2HAT_params.Local_param_s2hat.pixel_numbered_ring, number_pixels_S2HAT, &(Harmonic_sup->MAPPRAISER_to_S2HAT_butterfly), worldcomm);
    printf("%d »»»»»» ........ Butterfly communication performed !\n",rank); fflush(stdout);

    index = 0;
    printf("%d »»»»»» New check MAPP local map into ring map : From pixel 0 to %d - %f %d -", rank, max_size_test, local_pixel_map_MAPPRAISER_ring[index], local_pixel_indices_MAPPRAISER_ring[index]); fflush(stdout);
    for (index=1;index<max_size_test;index++){
        printf("- %f %d -", local_pixel_map_MAPPRAISER_ring[index], local_pixel_indices_MAPPRAISER_ring[index]);
    }
    printf("\n"); fflush(stdout);

    // if (number_pixels_S2HAT){
    if (Local_param_s2hat->gangrank >= 0){
        index = 0;
        check_monotony = monotony(Harmonic_sup->S2HAT_params.Local_param_s2hat.pixel_numbered_ring, number_pixels_S2HAT);
        // printf("%d »»»»»» Check obtained S2HAT map with size %d monotony_check %d : From pixel 0 to %d - %f %d -", rank, number_pixels_S2HAT, check_monotony, max_size_test, local_pixel_map_S2HAT[index], Harmonic_sup->S2HAT_params.Local_param_s2hat.pixel_numbered_ring[index]); fflush(stdout);
        for (index=1;index<max_size_test;index++){
            // printf("- %f %d -", local_pixel_map_S2HAT[index], Harmonic_sup->S2HAT_params.Local_param_s2hat.pixel_numbered_ring[index]);
        }
        // printf("\n"); fflush(stdout);
    }
    // printf("%d »»»»»» TEST S2HAT 0 ! \n", rank); fflush(stdout);
    // for (index=0;index<number_pixels_S2HAT;index++){
    //     if (!(local_pixel_map_S2HAT[index] == local_pixel_map_S2HAT[index])){
    //         number_of_nan++;
    //         }
    // }
    // printf("%d »»»»»» TEST S2HAT 1 ! \n", rank); fflush(stdout);
    // if (number_of_nan)
    //     printf("%d «««««««« NUMBER OF NANS NON-ZERO FOR number_pixels_S2HAT before apply_pix2alm : %d \n", rank, number_of_nan); fflush(stdout);
    
    printf("%d »»»»»» TEST S2HAT 2 ! \n", rank); fflush(stdout);
    if (Local_param_s2hat->gangrank >= 0)
        apply_pix2alm(local_pixel_map_S2HAT, local_alm_s2hat, &(Harmonic_sup->S2HAT_params));

    printf("%d »»»»»» TEST S2HAT 3 ! \n", rank); fflush(stdout);
    // for (index=0;index<Harmonic_sup->S2HAT_params.nstokes*Harmonic_sup->S2HAT_params.size_alm;index++){
    //     if (!(local_alm_s2hat[index].re == local_alm_s2hat[index].re)){
    //         number_of_nan++;
    //         }
    //     if (!(local_alm_s2hat[index].im == local_alm_s2hat[index].im)){
    //         number_of_nan++;
    //         }
    // }
    printf("%d »»»»»» TEST S2HAT 4 ! \n", rank); fflush(stdout);
    // if (number_of_nan)
    //     printf("%d «««««««« NUMBER OF NANS NON-ZERO within alms : %d \n", rank, number_of_nan); fflush(stdout);
    
    printf("%d »»»»»» TEST S2HAT 5 ! \n", rank); fflush(stdout);
    if (number_pixels_S2HAT)
        free(local_pixel_map_S2HAT);
    printf("%d »»»»»» Free step 0 ! \n", rank); fflush(stdout);
    free(local_pixel_indices_MAPPRAISER_ring);
    free(local_pixel_map_MAPPRAISER_ring);
    printf("%d »»»»»» Done ! \n", rank); fflush(stdout);
    return 0;
}

int global_harmonic_2_map(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, Harmonic_superstruct *Harmonic_sup){
//     // Transform harmonic S2HAT a_lm distribution into a local_maps pixel distribution from MAPPRAISER

    // S2HAT_GLOBAL_parameters *Global_param_s2hat = &(Harmonic_sup->S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(Harmonic_sup->S2HAT_params.Local_param_s2hat);
    // Butterfly_superstruct *Butterfly_map2harmonic = &(Harmonic_sup->MAPPRAISER_to_S2HAT_butterfly);
    
    MPI_Comm worldcomm = A->comm;

    int rank, size;
    MPI_Comm_size(worldcomm, &size);
    MPI_Comm_rank(worldcomm, &rank);

    int nstokes = Harmonic_sup->S2HAT_params.nstokes;
    double *local_pixel_map_S2HAT;
    
    int number_pixels_S2HAT = nstokes * Harmonic_sup->S2HAT_params.Local_param_s2hat.map_size; // Local map size
    if (number_pixels_S2HAT)
        local_pixel_map_S2HAT = (double *)calloc(number_pixels_S2HAT,sizeof(double));
    
    printf("%d ««««« Applying alm2pix \n", rank); fflush(stdout);
    if (Local_param_s2hat->gangrank >= 0)
        apply_alm2pix(local_pixel_map_S2HAT, local_alm_s2hat, &(Harmonic_sup->S2HAT_params));

    int number_pixels_MAPPRAISER = A->lcount-(A->nnz)*(A->trash_pix);
    double *local_pixel_map_MAPPRAISER_ring = (double *)calloc(number_pixels_MAPPRAISER,sizeof(double));
    int *local_pixel_indices_MAPPRAISER_ring = (int *)calloc(number_pixels_MAPPRAISER,sizeof(int));
    // double *local_pixel_map_MAPPRAISER_nest = (double *)calloc(true_size_local,sizeof(double));

    int *local_pixel_indices_MAPPRAISER_ring_tmp = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));
    printf("%d ««««« Projecting indices nest2ring \n", rank); fflush(stdout);
    convert_indices_nest2ring(A->lindices + (A->nnz) * (A->trash_pix), local_pixel_indices_MAPPRAISER_ring_tmp, number_pixels_MAPPRAISER, nstokes, Harmonic_sup->S2HAT_params.Global_param_s2hat.nside);
    project_int_values_into_different_scheme(local_pixel_indices_MAPPRAISER_ring_tmp, number_pixels_MAPPRAISER, Harmonic_sup->S2HAT_params.local_projector_values_nest2ring, local_pixel_indices_MAPPRAISER_ring);
    free(local_pixel_indices_MAPPRAISER_ring_tmp);
    
    printf("%d ««««« Performing butterfly communication \n", rank); fflush(stdout);
    perform_butterfly_communication(local_pixel_map_S2HAT, Harmonic_sup->S2HAT_params.Local_param_s2hat.pixel_numbered_ring, number_pixels_S2HAT, local_pixel_map_MAPPRAISER_ring, local_pixel_indices_MAPPRAISER_ring, number_pixels_MAPPRAISER, &(Harmonic_sup->S2HAT_to_MAPPRAISER_butterfly), worldcomm);

    printf("%d ««««« End of butterfly communication \n", rank); fflush(stdout);
    project_values_into_different_scheme(local_pixel_map_MAPPRAISER_ring, number_pixels_MAPPRAISER, Harmonic_sup->S2HAT_params.local_projector_values_ring2nest, local_pixel_map_MAPPRAISER);
    // local_pixel_indices_MAPPRAISER_ring
    // m2m(local_pixel_map_MAPPRAISER_ring, , number_pixels_MAPPRAISER, local_pixel_map_MAPPRAISER, A->lindices + (A->nnz) * (A->trash_pix), number_pixels_MAPPRAISER);

    printf("%d ««««« Free step 0 !\n", rank); fflush(stdout);
    if (number_pixels_S2HAT)
        free(local_pixel_map_S2HAT);

    printf("%d ««««« Free step 1 !\n", rank); fflush(stdout);
    free(local_pixel_indices_MAPPRAISER_ring);
    printf("%d ««««« Free step 2 !\n", rank); fflush(stdout);
    free(local_pixel_map_MAPPRAISER_ring);
    printf("%d ««««« Done \n", rank); fflush(stdout);
    return 0;
}

int get_mask_from_indices(Mat *A, int *mask_binary, int nside, int root)
{
    /* Get binary mask from indices observed, given by Mat */
    int i;
    long int number_pixels_total = 12*nside*nside;
    
    // int *all_sky_pixels_observed = (int *)malloc(number_pixels_total*sizeof(int));

    // A->lindices + (A->nnz) * (A->trash_pix)
    // A->lcount-(A->nnz)*(A->trash_pix)

    int *indices_ring = (int *)malloc((A->lcount-(A->nnz)*(A->trash_pix))*sizeof(int));

    convert_indices_nest2ring(A->lindices + (A->nnz) * (A->trash_pix), indices_ring, A->lcount-(A->nnz)*(A->trash_pix), A->nnz, nside);

    all_reduce_to_all_indices_mappraiser(indices_ring, (A->lcount-(A->nnz)*(A->trash_pix))/(A->nnz), nside, mask_binary, root, A->comm);

    free(indices_ring);

    for (i=0; i<number_pixels_total; i++)
    {
        if (mask_binary[i] != 0)
            mask_binary[i] = 1;
    }
    return 0;
}


int free_harmonic_superstruct(Harmonic_superstruct *Harmonic_sup, int rank)
{
    free_butterfly_superstruct(&(Harmonic_sup->S2HAT_to_MAPPRAISER_butterfly), rank);
    free_butterfly_superstruct(&(Harmonic_sup->MAPPRAISER_to_S2HAT_butterfly), rank);

    free_s2hat_parameters_struct(&(Harmonic_sup->S2HAT_params));
    return 0;
}

// int free_PCG_var(PCG_var *PCG_var_obj)
// {
//     free(PCG_var_obj->local_map_pix);
//     // if (PCG_var_obj->S2HAT_parameters != NULL)
//     //     free_s2hat_parameters_struct(PCG_var_obj->S2HAT_parameters);
//     free(PCG_var_obj);
// }


