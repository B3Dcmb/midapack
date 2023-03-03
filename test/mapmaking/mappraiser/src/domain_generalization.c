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


int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix, int domain_PCG_computation, int bool_apply_filter, int nstokes, S2HAT_parameters *S2HAT_params)
{   // Initalize PCG_var structure with S2HAT parameters and flags, before attribution of the data

    PCG_variable->local_map_pix = local_map_pix; // Local map in pixel domain, with MAPPRAISER convention
    // PCG_variable->local_alm = local_alm; // Local alm in harmonic domain, with S2HAT convention

    // if ((local_alm == NULL) && ((domain_PCG_computation == 1) || (bool_apply_filter == 1)))
    // {
    //     PCG_variable->local_alm = (s2hat_dcomplex *) calloc( nstokes * S2HAT_params->size_alm, sizeof(s2hat_dcomplex));
    // }

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

int init_butterfly_communication(Butterfly_struct *Butterfly_obj, int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, int do_we_need_to_project_into_different_scheme, MPI_Comm comm)
{
    // Initialize the butterfly communication
    int size;
    MPI_Comm_size(comm, &size);

    Butterfly_obj->classic_or_reshuffle_butterfly = flag_classic_or_reshuffle_butterfly;
    // 0 if classic butterfly, e.g. if pixel distributions are the same before/after ; 1 if reshuffle butterfly

    Butterfly_obj->steps = log_2(size);
    Butterfly_obj->S = (int **)malloc(Butterfly_obj->steps * sizeof(int *)); // allocate sending maps tab
    Butterfly_obj->R = (int **)malloc(Butterfly_obj->steps * sizeof(int *)); // allocate receiving maps tab
    Butterfly_obj->nS = (int *)malloc(Butterfly_obj->steps * sizeof(int));   // allocate sending map sizes tab
    Butterfly_obj->nR = (int *)malloc(Butterfly_obj->steps * sizeof(int));   // allocate receiving map size tab
    // butterfly_init(A->lindices + (A->nnz) * (A->trash_pix), A->lcount - (A->nnz) * (A->trash_pix), A->R, A->nR, A->S, A->nS, &(A->com_indices), &(A->com_count), A->steps, comm);

    switch(flag_classic_or_reshuffle_butterfly)
    {
        case 0: // Classic butterfly
            butterfly_reduce_init(indices_in, count_in, Butterfly_obj->R, Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS, &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count), Butterfly_obj->steps, comm);

        case 1: // Reshuffle butterfly
            butterfly_reshuffle_init(indices_in, count_in, indices_out, count_out, Butterfly_obj->R, Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS, &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count), Butterfly_obj->steps, comm);
    }

    Butterfly_obj->do_we_need_to_project_into_different_scheme = do_we_need_to_project_into_different_scheme;
    return 0;
}


int init_harmonic_superstruct(int is_pixel_scheme_MAPPRAISER_ring, Mat *A, Files_path_WIENER_FILTER *Files_WF_struct, S2HAT_parameters *S2HAT_params, Butterfly_struct *MAPP2ring_butterfly, Butterfly_struct *ring2MAPP_butterfly, int *mask_binary)
{
    /* Initalize all structures necessary for harmonic structures : S2HAT_params for S2HAT operations, and the 2 Butterfly_struct MAPP2ring_butterfly and ring2MAPP_butterfly for communication purposes
       indices_out and count_out correspond to the pixel distribution with the scheme prefered for TOD operation, which we call the "PCG pixel distribution scheme"
       If the pixel distribution is not in ring distribution, then it is expected to be in nest pixel distribution and we need to change it

       Relies on the fact that nstokes == A->nnz
    */
    int i;
    int flag_classic_butterfly = 0, flag_reshuffle_butterfly = 1;
    long int number_pixels_local, number_pixel_total;
    long int *indices_local_MAPPRAISER_ring; // New numbering of pixels used by MAPPRAISER in ring scheme, contained in local proc
    long int *pixel_numbered_ring;
    
    int *indices_local_MAPPRAISER = A->lindices + (A->nnz) * (A->trash_pix); // Indices 
    int number_pixels_MAPPRAISER = A->lcount - (A->nnz) * (A->trash_pix);
    // Will pixels be co-added ? lindices won't

    init_s2hat_parameters_superstruct(Files_WF_struct, mask_binary, A->nnz, S2HAT_params, A->comm);
    // Initialize S2HAT structures

    map_size = S2HAT_params->Local_param_s2hat->map_size; // Local map size
    
    if (is_pixel_scheme_MAPPRAISER_ring != 0)
    {
        // We work with a ring distribution
        number_pixel_total = 12*nside*nside;
        // indices_local_MAPPRAISER_ring = indices_local_MAPPRAISER;
        for (i=0; i<number_pixels_MAPPRAISER; i++)
            indices_local_MAPPRAISER_ring[i] = (indices_local_MAPPRAISER[i]%(S2HAT_params->nstokes))*number_pixel_total + indices_local_MAPPRAISER[i]/(S2HAT_params->nstokes)
    }
    else // If the pixel distribution is not in ring distribution, then it is expected to be in nest pixel distribution and we need to change it
    {   
        // We work with a nest distribution
        indices_local_MAPPRAISER_ring = (long int *)malloc(map_size*sizeof(long int));

        get_projectors_ring_and_nest(indices_local_MAPPRAISER, indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, projector_ring2nest, projector_nest2ring, S2HAT_params->nstokes)
        // If nest distribution with MAPPRAISER convention, convert the indices into ring distribution with S2HAT convention
        
        ring2MAPP_butterfly->projector_values = projector_ring2nest;
        MAPP2ring_butterfly->projector_values = projector_nest2ring;
        // Define the projectors

        ring2MAPP_butterfly->ordered_indices = indices_local_MAPPRAISER_ring; // Get ordered indices ring for future purpose
        MAPP2ring_butterfly->ordered_indices = indices_local_MAPPRAISER_ring; // Get ordered indices ring for future purpose
    }
    
    MAPP2ring_butterfly->indices_local_MAPPRAISER_ring = indices_local_MAPPRAISER_ring;

    pixel_numbered_ring = Local_param_s2hat->pixel_numbered_ring;

    // init_butterfly_communication(TOD2alm_butterfly, pixel_numbered_ring, number_pixels_local, NULL, 0, flag_classic_butterfly, world_comm);
    // // Communication from TOD into pixel map in RING scheme, for harmonic operation purposes

    init_butterfly_communication(ring2MAPP_butterfly, pixel_numbered_ring, map_size, indices_local_MAPPRAISER, number_pixels_MAPPRAISER, flag_reshuffle_butterfly, A->comm);
    // Communication from ring scheme pixel map in the PCG pixel distribution scheme
    

    init_butterfly_communication(MAPP2ring_butterfly, indices_local_MAPPRAISER, number_pixels_MAPPRAISER, pixel_numbered_ring, map_size, flag_reshuffle_butterfly, A->comm);
    // Communication from ring scheme pixel map in the PCG pixel distribution scheme
}

int butterfly_communication(double *values_to_communicate, int *indices_in, int count_in, double *values_out, double *indices_out, int count_out, Butterfly_struct *Butterfly_obj, MPI_Comm comm)
{
    // Perform the butterfly communication
    int k;
    int nSmax, nRmax;
    double *com_val, *values_ordered;

    for (k = 0; k < Butterfly_obj->steps; k++) /*compute max communication buffer size*/
    {
        if (Butterfly_obj->nR[k] > nRmax)
            nRmax = Butterfly_obj->nR[k];
        if (Butterfly_obj->nS[k] > nSmax)
            nSmax = Butterfly_obj->nS[k];
    }

    /* Copy value */
    com_val = (double *)malloc(Butterfly_obj->com_count * sizeof(double));
    for (i = 0; i < Butterfly_obj->com_count; i++)
        com_val[i] = 0.0;

    if (Butterfly_obj->do_we_need_to_project_into_different_scheme)
    {
        // We work with a ring distribution
        m2m(values_to_communicate, indices_in, count_in, com_val, Butterfly_struct->com_indices, Butterfly_struct->com_count);
    }
    else // If the pixel distribution is not in ring distribution, then it is expected to be in nest pixel distribution and we need to change it
    {   
        // We work with a nest distribution
        values_ordered = (double *)malloc(count_in*sizeof(double));
        // Project the nest distribution into rings
        project_values_into_different_scheme(values_to_communicate, count_in, Butterfly_struct->projector_values, values_ordered);
        
        m2m(values_ordered, Butterfly_obj->ordered_indices, count_in, com_val, Butterfly_struct->com_indices, Butterfly_struct->com_count);
    }

    switch(flag_classic_or_reshuffle_butterfly)
    {
        case 0: // Classic butterfly
            modified_butterfly_reduce(Butterfly_obj->R, Butterfly_obj->nR, nRmax, Butterfly_obj->S, Butterfly_obj->nS, nSmax, com_val, Butterfly_obj->steps, comm);

        case 1: // Reshuffle butterfly
            butterfly_reshuffle(Butterfly_obj->R, Butterfly_obj->nR, nRmax, Butterfly_obj->S, Butterfly_obj->nS, nSmax, com_val, Butterfly_obj->steps, comm);
    }
    m2m(com_val, Butterfly_struct->com_indices, Butterfly_struct->com_count, values_out, indices_out, count_out);
    
    if (Butterfly_obj->do_we_need_to_project_into_different_scheme == 0)
        free(values_ordered);
    
    free(com_val);
}


int global_map_2_harmonic(Butterfly_struct *Butterfly_map2harmonic, double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, S2HAT_parameters *S2HAT_params){
//     // Transform local_maps pixel distribution from MAPPRAISER into a harmonic S2HAT a_lm distribution

    S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters Local_param_s2hat = *(S2HAT_params->Local_param_s2hat);

    double *local_pixel_map_ring = (double *)calloc(Local_param_s2hat->map_size, sizeof(double));


    butterfly_communication(local_pixel_map_MAPPRAISER, A->lindices + (A->nnz) * (A->trash_pix), A->lcount - (A->nnz) * (A->trash_pix), 
                            local_pixel_map_ring, Local_param_s2hat->pixel_numbered_ring, Local_param_s2hat->map_size, 
                            Butterfly_map2harmonic, A->comm);

    apply_pix2alm(local_pixel_map_ring, local_alm_s2hat, S2HAT_params);

    free(local_pixel_map_ring);

    return 0;
}

int global_harmonic_2_map(Butterfly_struct *Butterfly_harmonic2map, double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A,  S2HAT_parameters *S2HAT_params){
    S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters Local_param_s2hat = *(S2HAT_params->Local_param_s2hat);
    
    double *local_map_pix_ring = (double *)calloc(Local_param_s2hat->map_size, sizeof(double));
    double *local_pixel_map_MAPPRAISER_ring = (double *)calloc((A->lcount - (A->nnz) * (A->trash_pix)),sizeof(double));    

    apply_alm2pix(local_alm_s2hat, local_map_pix_ring, S2HAT_params);
    
    butterfly_communication(local_pixel_map_ring, Local_param_s2hat->pixel_numbered_ring, Local_param_s2hat->map_size, 
                            local_pixel_map_MAPPRAISER_ring, Butterfly_harmonic2map->ordered_indices, A->lcount - (A->nnz) * (A->trash_pix), 
                            Butterfly_harmonic2map, A->comm);

    project_values_into_different_scheme(local_pixel_map_MAPPRAISER_ring, A->lcount - (A->nnz) * (A->trash_pix), Butterfly_harmonic2map->projector_values, local_pixel_map_MAPPRAISER);

    free(local_map_pix_ring);
    free(local_pixel_map_MAPPRAISER_ring);

    return 0;
}

int free_Butterfly_struct(Butterfly_struct *Butterfly_obj)
{
    free(Butterfly_obj->com_indices);
    free(Butterfly_obj->R);
    free(Butterfly_obj->nR);
    free(Butterfly_obj->S); 
    free(Butterfly_obj->nS);

    if (Butterfly_obj->do_we_need_to_project_into_different_scheme)
    {
        free(Butterfly_obj->projector_values);
        free(Butterfly_obj->ordered_indices);
    }

    free(Butterfly_obj);
}

int free_PCG_var(PCG_var *PCG_var_obj)
{
    free(PCG_var_obj->local_map_pix);

    if (PCG_var_obj->S2HAT_parameters != NULL)
        free_s2hat_parameters_struct(PCG_var_obj->S2HAT_parameters);
    

    free(PCG_var_obj);
}


// /* Taken from greedyreduce in mapmat.c, case ALLREDUCE
//     Reduce all pixel sky maps distributed among different procs to a single full sky map
// */
// int all_reduce_to_single_map_mappraiser(Mat *A, PCG_var *PCG_variable, int nside, double* output_fullsky_map, int root){
// int all_reduce_to_single_map_mappraiser(Mat *A, double* local_map_pix, int nside, double* output_fullsky_map, int root){
//     int i;
//     double *lvalues;
//     lvalues = (double *) malloc((A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));	//allocate and set to 0.0 local values
//     // memcpy(lvalues, PCG_variable->local_map_pix, (A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));		//copy local values into result values
//     memcpy(lvalues, local_map_pix, (A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));		//copy local values into result values
//     double *com_val;

//     long int npix = 12*nside*nside;

//     com_val=(double *) calloc( 3*npix,sizeof(double)); // Npix or less because of trash_pix ? To check later

//     s2m(com_val, lvalues, &(A->lindices[(A->nnz)*(A->trash_pix)]), A->lcount-(A->nnz)*(A->trash_pix)); 
    
//     /*for(i=0; i < A->com_count; i++){
//         printf("%lf ", com_val[i]);
//     } */
//     // MPI_Allreduce(com_val, output_fullsky_map, 3*npix, MPI_DOUBLE, MPI_SUM, A->comm);	//maximum index
//     MPI_Reduce(com_val, output_fullsky_map, 3*npix, MPI_DOUBLE, MPI_SUM, root, A->comm);	//maximum index

//     free(lvalues);
//     free(com_val);

//     // Here, the output_fullsky_map is in nest ordering, ordering as [nstokes, npix] in column-wise ordering
// }

// int brute_force_transfer_local_maps(Mat *A, double* local_pixel_map_MAPPRAISER, double *local_pixel_map_s2hat, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
//     // Temporary transfer method, which compile every thing and redistribute everything

//     int nside = Global_param_s2hat.nside;
//     double* full_sky_map;

//     int npix = 12*nside*nside;

//     double *full_sky_map_ring = (double *) malloc( 3*npix *sizeof(double));

//     if (Local_param_s2hat.gangroot == Local_param_s2hat.gangrank){
//         double *full_sky_map_nest = (double *) malloc( 3*npix *sizeof(double));
//         all_reduce_to_single_map_mappraiser(A, local_pixel_map_MAPPRAISER, nside, full_sky_map_nest, Local_param_s2hat.gangroot);

        
//         compute_full_map_nest2ring(full_sky_map_nest, full_sky_map_ring, nside, A->nnz, npix);
//         // Change ordering of full_sky_map to have what S2HAT expects, from MAPPRAISER nest [nstokes, npix] to S2HAT ring [nstokes, npix]
//         free(full_sky_map_nest);
//     }

//     distribute_full_sky_map_into_local_maps_S2HAT(full_sky_map_ring, local_pixel_map_s2hat, Global_param_s2hat, Local_param_s2hat, A->nnz);
//     // Distribute global full sky map into local maps in each processors
//     free(full_sky_map_ring);
// }

// /* Taken from greedyreduce in mapmat.c, case ALLREDUCE
//     Reduce all pixel sky maps distributed among different procs to a single full sky map
// */
// int all_reduce_to_single_map_mappraiser(Mat *A, PCG_var *PCG_variable, int nside, double* output_fullsky_map, int root){


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
