#ifdef W_MPI
 #include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <chealpix.h>
#include "s2hat.h"
// #include "midapack.h"
#include "s2hat_tools.h"
// #include "domain_generalization.h"


int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix, int domain_PCG_computation, int bool_apply_filter)
{   // Initalize PCG_var structure with S2HAT parameters and flags, before attribution of the data

    PCG_variable->local_map_pix = local_map_pix; // Local map in pixel domain, with MAPPRAISER convention

    // Specifications of the computation
    PCG_variable->domain_PCG_computation = domain_PCG_computation; // Domain chosen for PCG computation : 0 pixel, 1 harmonic
    PCG_variable->bool_apply_filter = bool_apply_filter; // 0 no filter applied, 1 Wiener filter applied
    // PCG_variable->nstokes = nstokes; // Number of Stokes parameters for the computation, 3 usually means T, Q, U

    // S2HAT structures necessary for harmonic domain
    // PCG_variable->S2HAT_parameters = S2HAT_params;
    // Contains S2HAT_parameters->Global_param_s2hat and S2HAT_parameters->Local_param_s2hat

    // Initialize update flags to 0 "no update needed"
    // PCG_variable->does_map_pixel_need_update = 0; // 0 no update needed, 1 need update from local_alm
    // PCG_variable->does_local_alm_need_update = 0; // 0 no update needed, 1 need update from local_map_pix
    return 0;
}


int init_harmonic_superstruct(Mat *A, Harmonic_superstruct *Harm_struct, int *mask_binary)
{
    /* Initalize all structures necessary for harmonic structures : S2HAT_params for S2HAT operations, and the 
       Relies on the fact that nstokes == A->nnz
    */
    
    int i, size;
    int flag_reshuffle_butterfly = 1;
    int number_pixels_local, number_pixel_total;
    int *pixel_numbered_ring;
    MPI_comm worldcomm = A->comm;
    MPI_Comm_size(worldcomm, &size);

    init_s2hat_parameters_superstruct(Files_WF_struct, mask_binary, A->nnz, &(Harm_struct->S2HAT_params), world_comm);
    // Initialize S2HAT structures
    S2HAT_parameters *S2HAT_params = &(Harm_struct->S2HAT_params);
    int nstokes = S2HAT_params->nstokes;

    int *indices_local_MAPPRAISER_nest = A->lindices + (A->nnz) * (A->trash_pix); // Indices local MAPPRAISER in nest distribution
    int number_pixels_MAPPRAISER = A->lcount - (A->nnz) * (A->trash_pix);

    int ordered_indices_local_MAPPRAISER_ring = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));
    int *projector_ring2nest = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));
    int *projector_nest2ring = (int *)malloc(number_pixels_MAPPRAISER*sizeof(int));

    get_projectors_indices(indices_local_MAPPRAISER_nest, ordered_indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, nstokes, S2HAT_params->Global_param_s2hat.nside, projector_ring2nest, projector_nest2ring);
    // As we have nest distribution with MAPPRAISER convention, convert the indices into ring distribution with S2HAT convention

    S2HAT_params->local_projector_values_ring2nest = projector_ring2nest;
    S2HAT_params->local_projector_values_nest2ring = projector_nest2ring;
    
    int *indices_local_S2HAT_ring = S2HAT_params->Local_param_s2hat->pixel_numbered_ring; // Indices local S2HAT in ring distribution
    int number_pixels_S2HAT = nstokes*S2HAT_params->Local_param_s2hat->map_size; // Local map size

    prepare_butterfly_communication(ordered_indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, indices_local_S2HAT_ring, number_pixels_S2HAT, flag_reshuffle_butterfly, &(Harm_struct->MAPPRAISER_to_S2HAT_butterfly), worldcomm);
    prepare_butterfly_communication(indices_local_S2HAT_ring, number_pixels_S2HAT, ordered_indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, flag_reshuffle_butterfly, &(Harm_struct->S2HAT_to_MAPPRAISER_butterfly), worldcomm);

    /* Attributing the Harmonic superstruc variables*/
    // Harmonic_superstruct->S2HAT_parameters = S2HAT_params;
    // Attributing the S2HAT superstructure to Harmonic superstructure

    // Harmonic_superstruct->S2HAT_to_MAPPRAISER = ring2MAPP_butterfly;
    // Harmonic_superstruct->MAPPRAISER_to_S2HAT = MAPP2ring_butterfly;
    // Attributing the Butterfly structures to Harmonic superstructure
}

int global_map_2_harmonic(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, Harmonic_superstruct *Harmonic_sup){
//     // Transform local_maps pixel distribution from MAPPRAISER into a harmonic S2HAT a_lm distribution

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(Harmonic_sup->S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(Harmonic_sup->S2HAT_params.Local_param_s2hat);
    Butterfly_superstruct *Butterfly_map2harmonic = &(Harmonic_sup->MAPPRAISER_to_S2HAT_butterfly);
    
    MPI_Comm worldcomm = A->comm;

    int i, rank, size;
    MPI_Comm_size(worldcomm, &size);
    MPI_Comm_rank(worldcomm, &rank);

    double *local_pixel_map_S2HAT;

    int number_pixels_MAPPRAISER = A->lcount-(A->nnz)*(A->trash_pix);
    double *local_pixel_map_MAPPRAISER_ring = (double *)calloc(number_pixels_MAPPRAISER,sizeof(double));
    int *local_pixel_indices_MAPPRAISER_ring = (int *)calloc(number_pixels_MAPPRAISER,sizeof(int));
    // double *local_pixel_map_MAPPRAISER_nest = (double *)calloc(true_size_local,sizeof(double));

    
    project_values_into_different_scheme(local_pixel_map_MAPPRAISER, number_pixels_MAPPRAISER, Harmonic_sup->S2HAT_params->local_projector_values_nest2ring, local_pixel_map_MAPPRAISER_ring);
    project_values_into_different_scheme(A->lindices + (A->nnz) * (A->trash_pix), number_pixels_MAPPRAISER, Harmonic_sup->S2HAT_params->local_projector_values_nest2ring, local_pixel_indices_MAPPRAISER_ring);

    int number_pixels_S2HAT = nstokes*S2HAT_params->Local_param_s2hat->map_size; // Local map size
    if (number_pixels_S2HAT)
        local_pixel_map_S2HAT = (double *)calloc(number_pixels_S2HAT*sizeof(double));

    perform_butterfly_communication(local_pixel_map_MAPPRAISER_ring, local_pixel_indices_MAPPRAISER_ring, A->lcount-(A->nnz)*(A->trash_pix), local_pixel_map_S2HAT, Harmonic_sup->S2HAT_params->Local_param_s2hat->pixel_numbered_ring, number_pixels_S2HAT, &(Harm_struct->MAPPRAISER_to_S2HAT_butterfly), worldcomm);
    
    if (Local_param_s2hat->gangrank >= 0)
        apply_pix2alm(local_pixel_map_S2HAT, local_alm_s2hat, &(Harmonic_sup->S2HAT_params));
    
    if (number_pixels_S2HAT)
        free(local_pixel_map_S2HAT);

    free(local_pixel_indices_MAPPRAISER_ring);
    free(local_pixel_map_MAPPRAISER_ring);
    return 0;
}

int global_harmonic_2_map(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, Harmonic_superstruct *Harmonic_sup){
//     // Transform harmonic S2HAT a_lm distribution into a local_maps pixel distribution from MAPPRAISER

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(Harmonic_sup->S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(Harmonic_sup->S2HAT_params.Local_param_s2hat);
    Butterfly_superstruct *Butterfly_map2harmonic = &(Harmonic_sup->MAPPRAISER_to_S2HAT_butterfly);
    
    MPI_Comm worldcomm = A->comm;

    int i, rank, size;
    MPI_Comm_size(worldcomm, &size);
    MPI_Comm_rank(worldcomm, &rank);

    double *local_pixel_map_S2HAT;
    
    int number_pixels_S2HAT = nstokes*S2HAT_params->Local_param_s2hat->map_size; // Local map size
    if (number_pixels_S2HAT)
        local_pixel_map_S2HAT = (double *)calloc(number_pixels_S2HAT*sizeof(double));

    if (Local_param_s2hat->gangrank >= 0)
        apply_alm2pix(local_alm_s2hat, local_pixel_map_S2HAT, &(Harmonic_sup->S2HAT_params));

    int number_pixels_MAPPRAISER = A->lcount-(A->nnz)*(A->trash_pix);
    double *local_pixel_map_MAPPRAISER_ring = (double *)calloc(number_pixels_MAPPRAISER,sizeof(double));
    int *local_pixel_indices_MAPPRAISER_ring = (int *)calloc(number_pixels_MAPPRAISER,sizeof(int));
    // double *local_pixel_map_MAPPRAISER_nest = (double *)calloc(true_size_local,sizeof(double));

    project_values_into_different_scheme(A->lindices + (A->nnz) * (A->trash_pix), number_pixels_MAPPRAISER, Harmonic_sup->S2HAT_params->local_projector_values_nest2ring, local_pixel_indices_MAPPRAISER_ring);
    
    perform_butterfly_communication(local_pixel_map_S2HAT, Harmonic_sup->S2HAT_params->Local_param_s2hat->pixel_numbered_ring, number_pixels_S2HAT, local_pixel_map_MAPPRAISER_ring, local_pixel_indices_MAPPRAISER_ring, number_pixels_MAPPRAISER, &(Harm_struct->S2HAT_to_MAPPRAISER_butterfly), worldcomm);

    project_values_into_different_scheme(local_pixel_map_MAPPRAISER_ring, number_pixels_MAPPRAISER, Harmonic_sup->S2HAT_params->local_projector_values_ring2nest, local_pixel_map_MAPPRAISER);
    
    
    if (number_pixels_S2HAT)
        free(local_pixel_map_S2HAT);

    free(local_pixel_indices_MAPPRAISER_ring);
    free(local_pixel_map_MAPPRAISER_ring);
    return 0;
}

int get_mask_from_indices(Mat *A, int *mask_binary, int nside, int root)
{
    /* Get binary mask from indices observed, given by Mat */
    int i;
    long int number_pixels_total = 12*nside*nside;
    
    int *all_sky_pixels_observed = (int *)malloc(number_pixels_total*sizeof(int));

    A->lindices + (A->nnz) * (A->trash_pix)
    A->lcount-(A->nnz)*(A->trash_pix)

    int *indices_ring = (int *)malloc((A->lcount-(A->nnz)*(A->trash_pix))*sizeof);

    convert_indices_nest2ring(A->lindices + (A->nnz) * (A->trash_pix), indices_ring, A->lcount-(A->nnz)*(A->trash_pix), A->nnz, nside)

    all_reduce_to_all_indices_mappraiser(indices_ring, (A->lcount-(A->nnz)*(A->trash_pix))/(A->nnz), nside, mask_binary, root, A->comm);

    free(indices_ring);

    for (i=0; i<number_pixels_total; i++)
    {
        if (mask_binary[i] != 0)
            mask_binary = 1;
    }
}


int free_harmonic_superstruct(Harmonic_superstruct *Harmonic_sup, int rank)
{
    free_butterfly_superstruct(&(Harmonic_sup->S2HAT_to_MAPPRAISER), rank);
    free_butterfly_superstruct(&(Harmonic_sup->MAPPRAISER_to_S2HAT), rank);

    free_s2hat_parameters_struct(&(Harmonic_sup->S2HAT_parameters));
}

int free_PCG_var(PCG_var *PCG_var_obj)
{
    free(PCG_var_obj->local_map_pix);
    // if (PCG_var_obj->S2HAT_parameters != NULL)
    //     free_s2hat_parameters_struct(PCG_var_obj->S2HAT_parameters);
    free(PCG_var_obj);
}


