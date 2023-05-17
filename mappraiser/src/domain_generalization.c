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


int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix, int domain_PCG_computation, int bool_apply_filter, int nstokes)
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
    // PCG_variable->S2HAT_parameters = S2HAT_params;
    // Contains S2HAT_parameters->Global_param_s2hat and S2HAT_parameters->Local_param_s2hat

    // Initialize update flags to 0 "no update needed"
    // PCG_variable->does_map_pixel_need_update = 0; // 0 no update needed, 1 need update from local_alm
    // PCG_variable->does_local_alm_need_update = 0; // 0 no update needed, 1 need update from local_map_pix
    return 0;
}


int init_harmonic_superstruct(int is_pixel_scheme_MAPPRAISER_ring, Mat *A, Harmonic_superstruct *Harm_struct, int *mask_binary)
{
    /* Initalize all structures necessary for harmonic structures : S2HAT_params for S2HAT operations, and the 2 Butterfly_struct MAPP2ring_butterfly and ring2MAPP_butterfly for communication purposes
       indices_out and count_out correspond to the pixel distribution with the scheme prefered for TOD operation, which we call the "PCG pixel distribution scheme"
       If the pixel distribution is not in ring distribution, then it is expected to be in nest pixel distribution and we need to change it

       Relies on the fact that nstokes == A->nnz
    */
    
    int i, size;
    int flag_classic_butterfly = 0, flag_reshuffle_butterfly = 1;
    long int number_pixels_local, number_pixel_total;
    long int *indices_local_MAPPRAISER_ring; // New numbering of pixels used by MAPPRAISER in ring scheme, contained in local proc
    long int *pixel_numbered_ring;
    MPI_comm world_comm = A->comm;
    MPI_Comm_size(world_comm, &size);

    Butterfly_struct *MAPP2ring_butterfly = (Butterfly_struct *)malloc(1*sizeof(Butterfly_struct));
    Butterfly_struct *ring2MAPP_butterfly = (Butterfly_struct *)malloc(1*sizeof(Butterfly_struct)):
    
    int *indices_local_MAPPRAISER = A->lindices + (A->nnz) * (A->trash_pix); // Indices 
    int number_pixels_MAPPRAISER = A->lcount - (A->nnz) * (A->trash_pix);
    // Will pixels be co-added ? lindices won't

    init_s2hat_parameters_superstruct(Files_WF_struct, mask_binary, A->nnz, S2HAT_params, world_comm);
    // Initialize S2HAT structures
   

    map_size = S2HAT_params->Local_param_s2hat->map_size; // Local map size
    
    // if (is_pixel_scheme_MAPPRAISER_ring != 0)
    // {
    //     // We work with a ring distribution
    //     number_pixel_total = 12*nside*nside;
    //     // indices_local_MAPPRAISER_ring = indices_local_MAPPRAISER;
    //     for (i=0; i<number_pixels_MAPPRAISER; i++)
    //         indices_local_MAPPRAISER_ring[i] = (indices_local_MAPPRAISER[i]%(S2HAT_params->nstokes))*number_pixel_total + indices_local_MAPPRAISER[i]/(S2HAT_params->nstokes)
    // }
    // else // If the pixel distribution is not in ring distribution, then it is expected to be in nest pixel distribution and we need to change it
    // {   
    //     // We work with a nest distribution
    //     indices_local_MAPPRAISER_ring = (long int *)malloc(map_size*sizeof(long int));

    //     get_projectors_indices(indices_local_MAPPRAISER, indices_local_MAPPRAISER_ring, number_pixels_MAPPRAISER, projector_ring2MAPP, projector_MAPP2ring, S2HAT_params->nstokes)
    //     // If nest distribution with MAPPRAISER convention, convert the indices into ring distribution with S2HAT convention
        
    //     ring2MAPP_butterfly->projector_values = projector_ring2MAPP;
    //     MAPP2ring_butterfly->projector_values = projector_MAPP2ring;
    //     // Define the projectors

    //     ring2MAPP_butterfly->ordered_indices = indices_local_MAPPRAISER_ring; // Get ordered indices ring for future purpose
    //     MAPP2ring_butterfly->ordered_indices = indices_local_MAPPRAISER_ring; // Get ordered indices ring for future purpose
    // }

    // Mirroring the sizes of the indices and data which will be
    int new_number_pixels_local_MAPPRAISER = 0;
    mirror_data_butterfly(NULL, NULL, number_pixels_MAPPRAISER, NULL, NULL, &new_number_pixels_local_MAPPRAISER, MIRROR_SIZE, world_comm);
    new_number_pixels_local_MAPPRAISER += number_pixels_MAPPRAISER;

    int *new_indices_local_MAPPRAISER
    mirror_data_butterfly(NULL, NULL, number_pixels_MAPPRAISER, NULL, NULL, &new_number_pixels_local_MAPPRAISER, MIRROR_INDICES, world_comm);

    // MPI_Comm *comm_butterfly;
    int nb_rank_Butterfly = pow(2,log2(size));
    // mpi_create_subset(nb_rank_Butterfly, world_comm, comm_butterfly);
    construct_butterfly_struct(Butterfly_struct *Butterfly_obj, int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, int do_we_need_to_project_into_different_scheme, MPI_Comm worlcomm)

    indices_local_MAPPRAISER_ring = (long int *)malloc(map_size*sizeof(long int));

    get_projectors_indices(new_indices_local_MAPPRAISER, indices_local_MAPPRAISER_ring, new_number_pixels_local_MAPPRAISER, projector_ring2MAPP, projector_MAPP2ring, S2HAT_params->nstokes)
    // If nest distribution with MAPPRAISER convention, convert the indices into ring distribution with S2HAT convention
    
    ring2MAPP_butterfly->projector_values = projector_ring2MAPP;
    MAPP2ring_butterfly->projector_values = projector_MAPP2ring;
    // Define the projectors

    ring2MAPP_butterfly->ordered_indices = indices_local_MAPPRAISER_ring; // Get ordered indices ring for future purpose
    MAPP2ring_butterfly->ordered_indices = indices_local_MAPPRAISER_ring; // Get ordered indices ring for future purpose


    MAPP2ring_butterfly->indices_local_MAPPRAISER_ring = indices_local_MAPPRAISER_ring;

    pixel_numbered_ring = Local_param_s2hat->pixel_numbered_ring;

    // init_butterfly_communication(TOD2alm_butterfly, pixel_numbered_ring, number_pixels_local, NULL, 0, flag_classic_butterfly, world_comm);
    // // Communication from TOD into pixel map in RING scheme, for harmonic operation purposes

    init_butterfly_communication(ring2MAPP_butterfly, pixel_numbered_ring, map_size, new_indices_local_MAPPRAISER, new_number_pixels_local_MAPPRAISER, flag_reshuffle_butterfly, A->comm);
    // Communication from ring scheme pixel map in the PCG pixel distribution scheme
    

    init_butterfly_communication(MAPP2ring_butterfly, new_indices_local_MAPPRAISER, new_number_pixels_local_MAPPRAISER, pixel_numbered_ring, map_size, flag_reshuffle_butterfly, A->comm);
    // Communication from ring scheme pixel map in the PCG pixel distribution scheme


    
    /* Attributing the Harmonic superstruc variables*/
    Harmonic_superstruct->S2HAT_parameters = S2HAT_params;
    // Attributing the S2HAT superstructure to Harmonic superstructure

    Harmonic_superstruct->S2HAT_to_MAPPRAISER = ring2MAPP_butterfly;
    Harmonic_superstruct->MAPPRAISER_to_S2HAT = MAPP2ring_butterfly;
    // Attributing the Butterfly structures to Harmonic superstructure
}

int global_map_2_harmonic(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, Harmonic_superstruct *Harmonic_sup){
//     // Transform local_maps pixel distribution from MAPPRAISER into a harmonic S2HAT a_lm distribution

    S2HAT_GLOBAL_parameters *Global_param_s2hat = Harmonic_sup->S2HAT_params->Global_param_s2hat;
    S2HAT_LOCAL_parameters *Local_param_s2hat = Harmonic_sup->S2HAT_params->Local_param_s2hat;
    Butterfly_struct *Butterfly_map2harmonic = Harmonic_sup->MAPPRAISER_to_S2HAT;
    MPI_Comm world_comm = A->comm;

    int i, rank, size;
    MPI_Comm_size(worldcomm, &size);
    MPI_Comm_rank(worldcomm, &rank);
    int true_size_local = Harmonic_sup->Butterfly_struct_supplement->new_size_local;

    double *local_pixel_map_ring;
    
    double *local_pixel_map_MAPPRAISER_ring = (double *)calloc(true_size_local,sizeof(double));
    // double *local_pixel_map_MAPPRAISER_nest = (double *)calloc(true_size_local,sizeof(double));


    project_values_into_different_scheme(local_pixel_map_MAPPRAISER, A->lcount-(A->nnz)*(A->trash_pix), Harmonic_sup->Butterfly_struct_supplement->projector_values, local_pixel_map_MAPPRAISER_ring);


    if (Harmonic_sup->Butterfly_struct_supplement->size_from_mirror != 0){
        // If data was mirrored previously, unmirror data

        int size_mirror = Harmonic_sup->Butterfly_struct_supplement->size_from_mirror;
        int *indices_from_mirror = Harmonic_sup->Butterfly_struct_supplement->indices_from_mirror;
        int size_butterfly = pow(2, log2(size));

        double *unmirrored_local_pixel_map_ring = (double *)malloc(size_mirror*sizeof(double));
        mirror_data_butterfly(local_pixel_map_MAPPRAISER_ring, NULL, A->lcount - (A->nnz) * (A->trash_pix), unmirrored_local_pixel_map_ring, NULL, &(size_mirror), MIRROR_DATA, world_comm);

        if (rank < size_butterfly){
            double *new_map = (double *)malloc(true_size_local*sizeof(double));
            
            m2m(local_pixel_map_MAPPRAISER, A->lindices + (A->nnz) * (A->trash_pix), A->lcount - (A->nnz) * (A->trash_pix), new_map, Harmonic_sup->Butterfly_struct_supplement->ordered_indices, true_size_local);
            m2m(unmirrored_local_pixel_map_ring, indices_from_mirror, size_mirror, new_map, Harmonic_sup->Butterfly_struct_supplement->ordered_indices, true_size_local);
            // Reassemble the new map values after mirroring in a single map
            
            memcpy(local_pixel_map_MAPPRAISER_ring, new_map, true_size_local*sizeof(double));
        }
        free(unmirrored_local_pixel_map_ring);
    }
    // free(local_pixel_map_MAPPRAISER_nest);
    
    local_pixel_map_ring = (double *)calloc(Local_param_s2hat->map_size, sizeof(double));
    butterfly_communication(local_pixel_map_MAPPRAISER_ring, Harmonic_sup->Butterfly_struct_supplement->ordered_indices, true_size_local, 
                            local_pixel_map_ring, Local_param_s2hat->pixel_numbered_ring, Local_param_s2hat->map_size, 
                            Butterfly_map2harmonic);
    // butterfly_communication(local_pixel_map_MAPPRAISER, A->lindices + (A->nnz) * (A->trash_pix), A->lcount - (A->nnz) * (A->trash_pix), 
    //                         local_pixel_map_ring, Local_param_s2hat->pixel_numbered_ring, Local_param_s2hat->map_size, 
    //                         Butterfly_map2harmonic, A->comm);

    if (Local_param_s2hat->gangrank >= 0){
        apply_pix2alm(local_pixel_map_ring, local_alm_s2hat, Harmonic_sup->S2HAT_params);
    }

    free(local_pixel_map_ring);
    free(local_pixel_map_MAPPRAISER_ring);

    return 0;
}

int global_harmonic_2_map(s2hat_dcomplex *local_alm_s2hat, double* local_pixel_map_MAPPRAISER, Mat *A, Harmonic_superstruct *Harmonic_sup){
    S2HAT_GLOBAL_parameters *Global_param_s2hat = Harmonic_sup->S2HAT_params->Global_param_s2hat;
    S2HAT_LOCAL_parameters *Local_param_s2hat = Harmonic_sup->S2HAT_params->Local_param_s2hat;
    Butterfly_struct *Butterfly_harmonic2map = Harmonic_sup->MAPPRAISER_to_S2HAT;
    int true_size_local = Harmonic_sup->Butterfly_struct_supplement->new_size_local;

    MPI_Comm world_comm = A->comm;

    int rank, size;
    MPI_Comm_size(worldcomm, &size);
    MPI_Comm_rank(worldcomm, &rank);

    double *local_map_pix_ring = (double *)calloc(Local_param_s2hat->map_size, sizeof(double));
    double *local_pixel_map_MAPPRAISER_ring = (double *)calloc(Harmonic_sup->Butterfly_struct_supplement->new_size_local,sizeof(double));

    if (Local_param_s2hat->gangrank >= 0){
        apply_alm2pix(local_alm_s2hat, local_map_pix_ring, Harmonic_sup->S2HAT_params);
    }

    butterfly_communication(local_map_pix_ring, Local_param_s2hat->pixel_numbered_ring, Local_param_s2hat->map_size, 
                            local_pixel_map_MAPPRAISER_ring, Harmonic_sup->Butterfly_struct_supplement->ordered_indices, true_size_local, 
                            Butterfly_harmonic2map);


    if (Harmonic_sup->Butterfly_struct_supplement->size_from_mirror != 0){
        // If data was mirrored previously, unmirror data

        int size_mirror = Harmonic_sup->Butterfly_struct_supplement->size_from_mirror;
        int *indices_from_mirror = Harmonic_sup->Butterfly_struct_supplement->indices_from_mirror;
        int size_butterfly = pow(2, log2(size));

        double *true_local_pixel_map_ring = (double *)calloc(size_from_mirror, sizeof(double));
        if (rank < size_butterfly)
            m2m(local_pixel_map_MAPPRAISER_ring, Harmonic_sup->Butterfly_struct_supplement->ordered_indices, true_size_local, true_local_pixel_map_ring, indices_from_mirror, size_mirror);
        
        double *unmirrored_local_pixel_map_ring = (double *)malloc((A->lcount - (A->nnz) * (A->trash_pix))*sizeof(double));
        mirror_data_butterfly(true_local_pixel_map_ring, NULL, size_mirror, unmirrored_local_pixel_map_ring, NULL, &(A->lcount - (A->nnz) * (A->trash_pix)), UNMIRROR_DATA, world_comm);

        if (rank >= size_butterfly){
            memcpy(local_pixel_map_MAPPRAISER_ring, unmirrored_local_pixel_map_ring, (A->lcount - (A->nnz) * (A->trash_pix))*sizeof(double));
        }
        free(unmirrored_local_pixel_map_ring);
        free(true_local_pixel_map_ring);
    }

    project_values_into_different_scheme(local_pixel_map_MAPPRAISER_ring, A->lcount - (A->nnz) * (A->trash_pix), Harmonic_sup->Butterfly_struct_supplement->projector_values, local_pixel_map_MAPPRAISER);

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
    // if (Butterfly_obj->do_we_need_to_project_into_different_scheme)
    // {
    //     free(Butterfly_obj->projector_values);
    //     free(Butterfly_obj->ordered_indices);
    // }
    free(Butterfly_obj->projector_values);
    free(Butterfly_obj->ordered_indices);

    free(Butterfly_obj);
}

int free_harmonic_superstruct(Harmonic_superstruct *Harmonic_sup)
{
    if (Harmonic_sup->S2HAT_to_MAPPRAISER != NULL){
        free_Butterfly_struct(Harmonic_sup->S2HAT_to_MAPPRAISER);
        free_Butterfly_struct(Harmonic_sup->MAPPRAISER_to_S2HAT);
    }
    if (Harmonic_sup->S2HAT_parameters != NULL)
        free_s2hat_parameters_struct(Harmonic_sup->S2HAT_parameters);
}

int free_PCG_var(PCG_var *PCG_var_obj)
{
    free(PCG_var_obj->local_map_pix);
    // if (PCG_var_obj->S2HAT_parameters != NULL)
    //     free_s2hat_parameters_struct(PCG_var_obj->S2HAT_parameters);
    free(PCG_var_obj);
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
