#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
#include <lapacke.h>
// #include "fitsio.h"
// #include <mpi.h>
#include <unistd.h>
// #include "s2hat.h"
#include <fitsio.h>
#include <chealpix.h>

#include "mappraiser.h"

void apply_Wiener_filter_pixel(int nside, int lmax, int nstokes, double *CMB_map_ring, double *CMB_map_output, double *c_ells, int number_correlations, double *mask_binary, MPI_Comm worldcomm){
    /* Apply Wiener filter in pixel space to CMB map */

    printf("############################################################ START !!! ############################################################ \n");
    // char *path_CMB_map = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/Map_test_band_limited.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE_128.fits";
    char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_SO.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_cut_latS.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_cut_equator.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_pt_src.fits";
    

    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTE_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_onlyBB_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_onlyEEBB_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_band_limited_1024_0.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTE.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTTTE.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTTTE_v2.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128_woTE.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128_only_BB.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128_only_EEBB.fits";
    // char *path_mask = "/pscratch/sd/m/mag/Masks_files/mask_SO_64.fits";
    // char *path_mask = "/pscratch/sd/m/mag/Masks_files/mask_cut_latS_64.fits";
    // char *path_mask = "/pscratch/sd/m/mag/Masks_files/mask_cut_equator_64.fits";
    char *path_mask = "/pscratch/sd/m/mag/Masks_files/mask_cut_pt_src_64.fits";

    int rank, size;

    // MPI_Init();
    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &size);

    int i, j, ell_value;

    // int nside = 512;
    // int nside_2 = 64;
    // int lmax_2 = 2*nside_2;
    int npix = 12*nside*nside;
    
    double *mask_binary_=NULL;
    // double *mask_binary= (double *)malloc(npix*sizeof(double));
    // read_fits_mask(nside, mask_binary, path_mask, 1);
    
    int nstokes_2 = 3; // T, Q, U
    int number_correlations_2 = 4; // TT, EE, BB, TE
    // int nstokes = 2; //Q, U
    // int number_correlations = 2; // EE, BB


    
    int index, index_2;
    int order_matrix = nstokes;

    // double *CMB_map_ring_2 = (double *) malloc( 3*npix*sizeof(double));
    // read_TQU_maps(nside_2, CMB_map_ring_2, path_CMB_map, 3);


    int max_size_pixels = 40;
    // printf("%d --- Reading first elems of map \n", rank);
    // for (i=0; i<max_size_pixels; i++)
    //     printf("- %f %f %f -", CMB_map_ring_2[i], CMB_map_ring_2[i+npix], CMB_map_ring_2[i+2*npix]);
    // printf("\n");

    printf("%d --- Initializing S2HAT_param \n", rank); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations_2);
    
    
    
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary_, nstokes_2, &S2HAT_params, worldcomm);
    printf("%d --- End of initializing S2HAT_param \n", rank); fflush(stdout);

    double *local_map_pix = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("%d --- Distributing S2HAT_param \n", rank); fflush(stdout);
    
    double *CMB_map_ring_2 = (double *)malloc(nstokes*npix*sizeof(double));
    memcpy(CMB_map_ring_2, CMB_map_ring, nstokes*npix*sizeof(double));
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring_2, local_map_pix, &S2HAT_params);
    // free(CMB_map_ring_2);
    printf("%d --- End of distributing S2HAT_param \n", rank); fflush(stdout);
    

    s2hat_dcomplex *local_alm_s2hat;

    local_alm_s2hat = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);
    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);


    printf("%d ----!!!! ***local_map_output*** nstokes %d ; map_size %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    double *local_map_output = (double *)calloc(nstokes*S2HAT_params.Local_param_s2hat.map_size,sizeof(double));


    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);


    memcpy(CMB_map_output, local_map_output, nstokes*npix*sizeof(double));
    // memcpy(CMB_map_output, CMB_map_ring, nstokes*npix*sizeof(double));

    printf("%d --- Free steps 0 ! \n", rank); fflush(stdout);
    // free(local_map_pix);
    printf("%d --- Free steps 1 ! \n", rank); fflush(stdout);
    free(local_map_output);

    printf("%d --- Done !!! \n", rank); fflush(stdout);
}


void apply_Wiener_filter_pixel_test(int nside, int lmax, int nstokes, double *CMB_map_ring, double *CMB_map_output, double *c_ells, int number_correlations, double *mask_binary, MPI_Comm worldcomm){
    /* Apply Wiener filter in pixel space to CMB map */

    char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";

    
    printf("############################################################ START !!! ############################################################ \n"); fflush(stdout);
    int rank, size;

    int npix = 12*nside*nside;
    // MPI_Init();
    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &size);

    // memcpy(CMB_map_output, CMB_map_ring, nstokes*npix*sizeof(double));

    printf("%d --- Initializing S2HAT_param -- nside %d lmax %d nstokes %d \n", rank, nside, lmax, nstokes); fflush(stdout);
    char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128.fits";
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);

    double *mask_binary_ = NULL;
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary_, nstokes, &S2HAT_params, worldcomm);

    double *local_map_pix = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));

    double *CMB_map_ring_2 = (double *)malloc(nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring_2, path_CMB_map, 3);
    // printf("%d --- Test map_size %d npix %d \n", rank, S2HAT_params.Local_param_s2hat.map_size, npix);
    // distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring, local_map_pix, &S2HAT_params);
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring_2, local_map_pix, &S2HAT_params);
    free(CMB_map_ring_2);
    // memcpy(CMB_map_output, local_map_pix, nstokes*npix*sizeof(double));

    s2hat_dcomplex *local_alm_s2hat = (s2hat_dcomplex *)calloc(nstokes*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);
    printf("%d --- map2harmonic-end --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);


    printf("%d ----!!!! ***local_map_output*** nstokes %d ; map_size %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));


    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);
    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);

    memcpy(CMB_map_output, local_map_output, nstokes*npix*sizeof(double));

    printf("%d --- Free steps 0 ! \n", rank); fflush(stdout);
    free(local_map_pix);
    printf("%d --- Free steps 1 ! \n", rank); fflush(stdout);
    free(local_map_output);

    printf("%d --- Done !!! \n", rank); fflush(stdout);
}


void apply_Wiener_filter_pixel_true(int nside, int lmax, int nstokes, double *CMB_map_ring, double *CMB_map_output, double *c_ells, int number_correlations, double *mask_binary, MPI_Comm worldcomm){
    /* Apply Wiener filter in pixel space to CMB map */

    printf("############################################################ START !!! ############################################################ \n"); fflush(stdout);
    int rank, size;

    // MPI_Init();
    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &size);

    int ell_value;
    // int i, j;

    int npix = 12*nside*nside;
    // int index, index_2;
    // int order_matrix = nstokes;

    printf("%d --- Initializing S2HAT_param -- nside %d lmax %d nstokes %d \n", rank, nside, lmax, nstokes); fflush(stdout);
    char *c_ell_path = "";
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);

    double *local_map_pix = (double *)calloc(nstokes*S2HAT_params.Local_param_s2hat.map_size,sizeof(double));

    double *CMB_map_output_copy_0 = (double *)calloc(nstokes*npix,sizeof(double));
    memcpy(CMB_map_output_copy_0, CMB_map_ring, nstokes*npix*sizeof(double));
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_output_copy_0, local_map_pix, &S2HAT_params);
    free(CMB_map_output_copy_0);

    s2hat_dcomplex *local_alm_s2hat;

    local_alm_s2hat = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    int i, max_size_pixels=40;
    printf("%d --- Reading first elems of local map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", local_map_pix[i], CMB_map_ring[i]);
    printf("\n");

    // S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    // S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);
    printf("%d --- map2harmonic-end --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    // free(local_map_pix);

    double **inverse_covariance_matrix;
    inverse_covariance_matrix = malloc((lmax+1)*sizeof(double *));
    for(ell_value=0; ell_value<lmax+1; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);

    get_covariance_matrix_block_diagonal(c_ells, number_correlations, inverse_covariance_matrix, &S2HAT_params);

    int index, index_2, order_matrix=nstokes, size_ell_max=10;
    for (ell_value=0;ell_value<size_ell_max; ell_value++){
        printf("\n #### ell_init= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
        printf("\n");
    }
    }
    for(ell_value=0; ell_value<lmax+1; ell_value++){
        get_cholesky_decomposition_inverted(nstokes, inverse_covariance_matrix[ell_value], 'L');
    }

    for (ell_value=0;ell_value<size_ell_max; ell_value++){
        printf("\n #### ell_Inverted= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
        printf("\n");
    }
    }
    // double normalization_factor = 0;
    // for (ell_value=0; ell_value<lmax; ell_value++){
    //     normalization_factor += (2*ell_value + 1)/(4*M_PI);
    // }

    // printf("Value of normalization factor %f \n", normalization_factor);
    // normalization_factor = normalization_factor;
    // char cholesky_part = 'L';
    // int info;

    // for(ell_value=0; ell_value<lmax; ell_value++){
    //     for (index=0; index < order_matrix; index++){
    //         for (index_2=0; index_2<order_matrix; index_2++){
    //             inverse_covariance_matrix[ell_value][index*order_matrix + index_2] *= 1/normalization_factor;
    //         }
    //     }
    //     dpotrf_(&cholesky_part, &nstokes, inverse_covariance_matrix[ell_value], &nstokes, &info);

    //     inverse_covariance_matrix[ell_value][nstokes] = 0;
    //     if (nstokes == 3){
    //         inverse_covariance_matrix[ell_value][6] = 0;
    //         inverse_covariance_matrix[ell_value][7] = 0;
    //     }
    // }

    printf("%d ----!!!! ***local_map_output*** nstokes %d ; map_size %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));

    s2hat_dcomplex *local_alm_s2hat_inverted;

    local_alm_s2hat_inverted = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);
    // apply_inv_block_diag_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, &S2HAT_params);

    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    
    // apply_alm2pix(local_map_output, local_alm_s2hat_inverted, &S2HAT_params);
    apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);

    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);



    double *CMB_map_output_copy = (double *)calloc(nstokes*npix,sizeof(double));
    gather_map(local_map_output, CMB_map_output_copy, nstokes, &S2HAT_params);
    memcpy(CMB_map_output, CMB_map_output_copy, nstokes*npix*sizeof(double));
    free(CMB_map_output_copy);
    
    printf("%d --- Reading first elems of local map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", CMB_map_output[i], CMB_map_ring[i]);
    printf("\n");
    printf("%d --- Reading first elems of local map + npix \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", CMB_map_output[i+npix], CMB_map_ring[i+npix]);
    printf("\n");
    printf("%d --- Reading first elems of local map + 2*npix \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", CMB_map_output[i+2*npix], CMB_map_ring[i+2*npix]);
    printf("\n");


    double average_relative_error_T_0_1 = 0;
    double average_relative_error_Q_0_1 = 0;
    double average_relative_error_U_0_1 = 0;
    double max_average_relative_error_T_0_1 = 0;
    double max_average_relative_error_Q_0_1 = 0;
    double max_average_relative_error_U_0_1 = 0;

    int size_total_probed = npix;
    // size_total_probed = 10;

    int init_index = 0;
    int number_0 = 0;

    index = init_index;
    max_average_relative_error_T_0_1 = fabs((CMB_map_ring[index]-CMB_map_output[index])/CMB_map_ring[index]);
    max_average_relative_error_Q_0_1 = fabs((CMB_map_ring[index+npix]-CMB_map_output[index+npix])/CMB_map_ring[index+npix]);
    max_average_relative_error_U_0_1 = fabs((CMB_map_ring[index+2*npix]-CMB_map_output[index+2*npix])/CMB_map_ring[index+2*npix]);
    for(index=1;index<size_total_probed+init_index;index++){
        average_relative_error_T_0_1 += fabs((CMB_map_ring[index]-CMB_map_output[index])/CMB_map_ring[index]);
        if (max_average_relative_error_T_0_1 < fabs((CMB_map_ring[index]-CMB_map_output[index])/CMB_map_ring[index]))
            max_average_relative_error_T_0_1 = fabs((CMB_map_ring[index]-CMB_map_output[index])/CMB_map_ring[index]);

        if(CMB_map_ring[index] == 0.000000){
            number_0++;
        }
        if(nstokes>1){
            average_relative_error_Q_0_1 += fabs((CMB_map_ring[index+npix]-CMB_map_output[index+npix])/CMB_map_ring[index+npix]);
            if (max_average_relative_error_Q_0_1 < fabs((CMB_map_ring[index+npix]-CMB_map_output[index+npix])/CMB_map_ring[index+npix]))
                max_average_relative_error_Q_0_1 = fabs((CMB_map_ring[index+npix]-CMB_map_output[index+npix])/CMB_map_ring[index+npix]);

            average_relative_error_U_0_1 += fabs((CMB_map_ring[index+2*npix]-CMB_map_output[index+2*npix])/CMB_map_ring[index+2*npix]);
            if (max_average_relative_error_U_0_1 < fabs((CMB_map_ring[index+2*npix]-CMB_map_output[index+2*npix])/CMB_map_ring[index+2*npix]))
                max_average_relative_error_U_0_1 = fabs((CMB_map_ring[index+2*npix]-CMB_map_output[index+2*npix])/CMB_map_ring[index+2*npix]);
        }
    }
    printf("\n");
    printf("Average 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            average_relative_error_T_0_1/(size_total_probed), average_relative_error_Q_0_1/(size_total_probed), average_relative_error_U_0_1/(size_total_probed));
    printf("Maximum 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            max_average_relative_error_T_0_1, max_average_relative_error_Q_0_1, max_average_relative_error_U_0_1);

    printf("Number of 0 : %d \n", number_0);
    fflush(stdout);
    

    printf("%d --- Free step 2 ! \n", rank); fflush(stdout);
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0){
        printf("%d --- Free step 1b ! \n", rank); fflush(stdout);
        free(local_alm_s2hat_inverted);
    }
    printf("%d --- Free step 2b ! \n", rank); fflush(stdout);
    // free(local_map_pix);
    free(local_map_output);
    printf("%d --- Free step 3 ! \n", rank); fflush(stdout);
    free_covariance_matrix(inverse_covariance_matrix, lmax);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    free_s2hat_parameters_struct(&S2HAT_params);

    printf("%d --- Done !!! \n", rank); fflush(stdout);
}
