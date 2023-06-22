#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
// #include <mpi.h>
#include <unistd.h>
// #include "s2hat.h"
// #include <chealpix.h>

#include "midapack.h"
// #include "s2hat.h"



int main_0(int argc, char** argv){
// int main(int argc, char** argv){
    s2hat_pixparameters pixpar;
    pixpar.par1 = 16;
    int pixchoice = PIXCHOICE_HEALPIX;
    s2hat_pixeltype pixelization_scheme;
    set_pixelization(pixchoice, pixpar, &pixelization_scheme);
    printf("Hi ! \n"); fflush(stdout);
    return 0;
}

// int main_indices_monotony(int argc, char** argv){
int main(int argc, char** argv){

    char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_band_limited_1024_0.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits";
    int rank, nprocs;
    int i;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int check_monotony;
    // int nside = 512;
    int nside = 2;
    int lmax = 2*nside+2;
    int domain_PCG_computation = 0;
    int bool_apply_filter = 0;
    int nstokes = 3;
    int number_correlations = 4;


    int npix = 12*nside*nside;
    int index = 0;
    int max_size_test = 10;

    int lcount = (npix/nprocs + 10)*nstokes;
    

    if (rank == nprocs-1)
        lcount += (npix%nprocs - 10)*nstokes;
    printf("%d --- Producing lcount !!! %d over %d\n", rank, lcount, nstokes*npix); fflush(stdout);

    int number_pixels_MAPP = lcount;
    
    int *lindices = (int *)malloc(number_pixels_MAPP*sizeof(int));

    int first_pixel = rank*(npix/nprocs)*nstokes;
    printf("%d --- Producing lindices !!! first_pixel %d --- lindices -", rank, first_pixel);
    for (i=0; i<lcount; i++){
        lindices[i] = first_pixel + i;
        printf("- %d -", lindices[i]);
    }
    printf("\n"); fflush(stdout);

    check_monotony = monotony_v2(lindices, lcount);

    printf("%d --- FIRST check_monotony lindices : %d \n", rank, check_monotony);
    
    Files_path_WIENER_FILTER Files_WF_struct;
    printf("%d --- Initializing struct_WF \n", rank); fflush(stdout);
    init_files_struct_WF(&Files_WF_struct, nside, lmax, c_ell_path, number_correlations);
    int *mask_binary = NULL;
    S2HAT_parameters S2HAT_params;
    init_s2hat_parameters_superstruct(&Files_WF_struct, mask_binary, nstokes, &S2HAT_params, worldcomm);

    int *projector_ring2nest = (int *)malloc(number_pixels_MAPP*sizeof(int));
    int *projector_nest2ring = (int *)malloc(number_pixels_MAPP*sizeof(int));
    int *ordered_indices_local_ring = (int *)malloc(number_pixels_MAPP*sizeof(int));
    int *unordered_indices_local_ring = (int *)malloc(number_pixels_MAPP*sizeof(int));
    int *ordered_indices_local_ring_bis = (int *)malloc(number_pixels_MAPP*sizeof(int));

    int *indices_ring = (int *)malloc(number_pixels_MAPP*sizeof(int));
    convert_indices_nest2ring(lindices, indices_ring, number_pixels_MAPP, nstokes, nside);
    memcpy(unordered_indices_local_ring, indices_ring, number_pixels_MAPP*sizeof(int));
    
    free(indices_ring);

    get_projectors_indices(lindices, ordered_indices_local_ring, number_pixels_MAPP, nstokes, nside, projector_ring2nest, projector_nest2ring);

    int *local_pixel_indices_MAPPRAISER_ring = (int *)malloc(number_pixels_MAPP*sizeof(int));
    project_int_values_into_different_scheme(lindices, number_pixels_MAPP, projector_nest2ring, local_pixel_indices_MAPPRAISER_ring);
    // project_int_values_into_different_scheme(lindices, number_pixels_MAPP, projector_ring2nest, local_pixel_indices_MAPPRAISER_ring);
    project_int_values_into_different_scheme(unordered_indices_local_ring, number_pixels_MAPP, projector_nest2ring, ordered_indices_local_ring_bis);

    check_monotony = monotony_v2(ordered_indices_local_ring, lcount);
    printf("%d --- SECOND check_monotony ordered_indices_local_ring : %d \n", rank, check_monotony);
    check_monotony = monotony_v2(local_pixel_indices_MAPPRAISER_ring, lcount);
    printf("%d --- THIRD check_monotony local_pixel_indices_MAPPRAISER_ring -- projected : %d \n", rank, check_monotony);

    check_monotony = monotony_v2(ordered_indices_local_ring_bis, lcount);
    printf("%d --- THIRDbis check_monotony ordered_indices_local_ring_bis -- projected : %d \n", rank, check_monotony);

    printf("%d --- local_pixel_indices_MAPPRAISER_ring complete : size %d !!! -", rank, lcount);
    for (i=0; i<lcount/nstokes; i++){
        printf("- *%d* %d %d %d %d -", i, ordered_indices_local_ring[i], ordered_indices_local_ring_bis[i], local_pixel_indices_MAPPRAISER_ring[i], unordered_indices_local_ring[i]);
    }
    printf("\n"); fflush(stdout);

    printf("%d --- projectors complete : size %d !!! -", rank, lcount);
    for (i=0; i<lcount/nstokes; i++){
        // index = i*nstokes;
        printf("- %i pr2n %d pn2r %d -", i, projector_ring2nest[i], projector_nest2ring[i]);
    }
    printf("\n"); fflush(stdout);


    int *s2hat_indices = S2HAT_params.Local_param_s2hat.pixel_numbered_ring;
    check_monotony = monotony_v2(s2hat_indices, nstokes * S2HAT_params.Local_param_s2hat.map_size);
    printf("%d --- FOURTH check_monotony s2hat_indices : %d with size %d true_size %d \n", rank, check_monotony, S2HAT_params.Local_param_s2hat.map_size*nstokes, S2HAT_params.Local_param_s2hat.map_size);

    printf("%d --- s2hat indices !!! first : %d ; last map_size/2-1: %d ; last map_size/2: %d ; last map_size: %d ; last nstokes*map_size %d \n", rank, s2hat_indices[0], s2hat_indices[S2HAT_params.Local_param_s2hat.map_size/2-1], s2hat_indices[S2HAT_params.Local_param_s2hat.map_size/2], s2hat_indices[S2HAT_params.Local_param_s2hat.map_size-1], s2hat_indices[S2HAT_params.Local_param_s2hat.map_size*nstokes-1]);

    // if (rank==nprocs-1){
    //     printf("%d --- Last proc s2hat indices : size %d !!! -", rank, nstokes * S2HAT_params.Local_param_s2hat.map_size);
    //     for (i=0; i<nstokes * S2HAT_params.Local_param_s2hat.map_size; i++){
    //         printf("- %d -", s2hat_indices[i]);
    //     }
    //     printf("\n"); fflush(stdout);
    // }

    // free_PCG_var(PCG_variable);
    printf("%d --- Free step \n", rank); fflush(stdout);

    printf("%d --- Done !!! \n", rank); fflush(stdout);
    free(local_pixel_indices_MAPPRAISER_ring);
    free(projector_nest2ring);
    free(projector_ring2nest);
    free(lindices);
    free(ordered_indices_local_ring);
    free(unordered_indices_local_ring);
    free(ordered_indices_local_ring_bis);
    free_s2hat_parameters_struct(&S2HAT_params);

    return 0;
}

