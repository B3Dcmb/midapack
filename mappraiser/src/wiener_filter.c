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


void apply_Wiener_filter_pixel(int nside, int lmax, int nstokes, double *CMB_map, double *CMB_map_output, double *c_ells, int number_correlations, double *mask_binary, MPI_Comm worldcomm)
{
    // MPI_Init( &argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &nprocs);

    int npix = 12*nside*nside;
    int index, index_2;

    char *c_ell_path = "";
    

    // if (rank==0)
        // printf("Wiener-filter pixel application ! ############################################################################################################################ \n");
    // printf("Initializing MPI %d %d \n", rank, nprocs); fflush(stdout);

    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);

    // printf("Initializing S2HAT_param \n"); fflush(stdout);
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);
    // printf("Finish initializing S2HAT_param !! \n"); fflush(stdout);

    S2HAT_LOCAL_parameters Local_param_s2hat = S2HAT_params.Local_param_s2hat;

    if (Local_param_s2hat.gangrank >= 0){
        // printf("Entering Wiener-filter with rank %d \n", rank); fflush(stdout);
        
        double *CMB_map_temp  = (double *)malloc(3*npix*sizeof(double));
        // double CMB_map_temp[3*npix];
        memcpy(CMB_map_temp, CMB_map, 3*npix*sizeof(double));

        // double *local_map_pix = (double *) calloc( 3*Local_param_s2hat.map_size, sizeof(double));
        // MPI_Barrier(Local_param_s2hat.gangcomm);
        double local_map_pix[3*Local_param_s2hat.map_size];
        distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_temp, local_map_pix, &S2HAT_params);
        free(CMB_map_temp);
        // MPI_Barrier(Local_param_s2hat.gangcomm);
        
        
        // s2hat_dcomplex *local_alm, *local_alm_out;
        // local_alm = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        // local_alm_out = (s2hat_dcomplex *) malloc( (3*S2HAT_params.size_alm)*sizeof(s2hat_dcomplex));
        // s2hat_dcomplex *local_alm_inverted = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        s2hat_dcomplex local_alm[3*S2HAT_params.size_alm];
        s2hat_dcomplex local_alm_inverted[3*S2HAT_params.size_alm];

        
        for (index=0; index<3*S2HAT_params.size_alm; index++){
            local_alm[index].re = 0;
            local_alm[index].im = 0;
            local_alm_inverted[index].re = 0;
            local_alm_inverted[index].im = 0;
        }
        
        // Doing Pix2Alm
        apply_pix2alm(local_map_pix, local_alm, &S2HAT_params);
        // free(local_map_pix);

        // Getting covariance matrix
        int ell_value;
        double **inverse_covariance_matrix;
        inverse_covariance_matrix = malloc((lmax+1)*sizeof(double *));
        // double **inverse_covariance_matrix_2;
        // inverse_covariance_matrix_2 = malloc((lmax+1)*sizeof(double *));
        
        for(ell_value=0; ell_value<lmax+1; ell_value++){
            inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
            // inverse_covariance_matrix_2[ell_value] = calloc(nstokes*nstokes,sizeof(double));
            // inverse_covariance_matrix_2[ell_value][0] = 1;
            // inverse_covariance_matrix_2[ell_value][4] = 1;
            // inverse_covariance_matrix_2[ell_value][8] = 1;
        }
        get_covariance_matrix_block_diagonal(c_ells, number_correlations, inverse_covariance_matrix, &S2HAT_params);

        char cholesky_part = 'L';
        int info;
        for(ell_value=0; ell_value<lmax+1; ell_value++){
            get_cholesky_decomposition_inverted(nstokes, inverse_covariance_matrix[ell_value], cholesky_part);
        }
        // double normalization_factor = 0;
        // for (ell_value=0; ell_value<lmax+1; ell_value++){
        //     normalization_factor += (2*ell_value + 1)/(4*M_PI);
        // }

        for(ell_value=0; ell_value<lmax+1; ell_value++){
            // for (index=0; index < nstokes; index++){
            //     for (index_2=0; index_2<nstokes; index_2++){
            //         inverse_covariance_matrix[ell_value][index*nstokes + index_2] *= 1/normalization_factor;
            //     }
            // }
            dpotrf_(&cholesky_part, &nstokes, inverse_covariance_matrix[ell_value], &nstokes, &info);

            inverse_covariance_matrix[ell_value][nstokes] = 0;
            if (nstokes == 3){
                inverse_covariance_matrix[ell_value][6] = 0;
                inverse_covariance_matrix[ell_value][7] = 0;
            }
        }

        apply_inv_block_diag_covariance_matrix_to_alm(local_alm, local_alm_inverted, inverse_covariance_matrix, &S2HAT_params);


        free_covariance_matrix(inverse_covariance_matrix, lmax);
        // free_covariance_matrix(inverse_covariance_matrix_2, lmax);



        // double *local_map_pix_out = (double *)calloc(3*Local_param_s2hat.map_size, sizeof(double));
        double local_map_pix_out[3*Local_param_s2hat.map_size];
        
        // apply_alm2pix(local_map_pix_out, local_alm, &S2HAT_params);
        apply_alm2pix(local_map_pix_out, local_alm_inverted, &S2HAT_params);
        // free(local_alm);
        // free(local_alm_out);
        

        double *full_sky_map_2;
        full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
        // MPI_Barrier(Local_param_s2hat.gangcomm);
        gather_map(local_map_pix_out, full_sky_map_2, nstokes, &S2HAT_params);
        // free(local_map_pix_out);
        // MPI_Barrier(Local_param_s2hat.gangcomm);

        if (rank==0)
            memcpy(CMB_map_output, full_sky_map_2, nstokes*npix*sizeof(double));

        // if (rank == 0){
        //     printf("//////////////////////// Comparison new map vs old \n");
        //     fflush(stdout);

        //     double average_relative_error_T_0_1 = 0;
        //     double average_relative_error_Q_0_1 = 0;
        //     double average_relative_error_U_0_1 = 0;
        //     double max_average_relative_error_T_0_1 = 0;
        //     double max_average_relative_error_Q_0_1 = 0;
        //     double max_average_relative_error_U_0_1 = 0;

        //     int size_total_probed = npix;

        //     int init_index = 0;
        //     int number_0 = 0;

        //     double pixel_value_T;
        //     double pixel_value_Q;
        //     double pixel_value_U;

        //     int index = init_index;

        //     max_average_relative_error_T_0_1 = fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
        //     max_average_relative_error_Q_0_1 = fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
        //     max_average_relative_error_U_0_1 = fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);
        //     for(index=init_index;index<size_total_probed+init_index;index++){
        //         pixel_value_T = fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
        //         average_relative_error_T_0_1 += pixel_value_T;
        //         if (max_average_relative_error_T_0_1 < pixel_value_T)
        //             max_average_relative_error_T_0_1 = pixel_value_T;
        //         if(CMB_map[index] == 0.000000){
        //             number_0++;
        //         }
        //         if(nstokes>1){
        //             pixel_value_Q = fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
        //             average_relative_error_Q_0_1 += pixel_value_Q;
        //             if (max_average_relative_error_Q_0_1 < pixel_value_Q)
        //                 max_average_relative_error_Q_0_1 = pixel_value_Q;
        //             if(CMB_map[index+npix] == 0.000000){
        //                 number_0++;
        //             }
        //             pixel_value_U = fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);
        //             average_relative_error_U_0_1 += pixel_value_U;
        //             if (max_average_relative_error_U_0_1 < pixel_value_U)
        //                 max_average_relative_error_U_0_1 = pixel_value_U;
        //             if(CMB_map[index+2*npix] == 0.000000){
        //                 number_0++;
        //             }
        //         }
        //     }
        //     printf("\n");
        //     printf("Average 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
        //             average_relative_error_T_0_1/(size_total_probed), average_relative_error_Q_0_1/(size_total_probed), average_relative_error_U_0_1/(size_total_probed));
        //     printf("Maximum 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
        //             max_average_relative_error_T_0_1, max_average_relative_error_Q_0_1, max_average_relative_error_U_0_1);
        //     printf("Number of 0 : %d \n", number_0);
        //     fflush(stdout);
        // }


        free(full_sky_map_2);
        // free(local_alm_inverted);

        free_s2hat_parameters_struct(&S2HAT_params);

    // MPI_Finalize();
    }
    // printf("Wiener-filter pixel done for rank %d ! \n", rank);
    fflush(stdout);
}


// NSIDE=64
// Average 0-1 relative error on all pixels with 3 49152 : T 0.025333 ; Q 0.003341 ; U 0.003203 
// Maximum 0-1 relative error on all pixels with 3 49152 : T 431.783585 ; Q 39.388716 ; U 38.044393 
// Average 0-1 relative error on all pixels with 3 49152 : T 0.002534 ; Q 0.003341 ; U 0.003203 
// Maximum 0-1 relative error on all pixels with 3 49152 : T 27.571147 ; Q 39.388716 ; U 38.044393 

// NSIDE=2
// Average 0-1 relative error on all pixels with 3 48 : T 1.239642 ; Q 0.190214 ; U 0.190214 // Without inv_cov = 1, 1 MPI_Task
// Maximum 0-1 relative error on all pixels with 3 48 : T 3.795595 ; Q 0.453338 ; U 0.453338 
// Average 0-1 relative error on all pixels with 3 48 : T 0.059915 ; Q 0.190214 ; U 0.190214 // Without inv_cov = 1, 2 MPI_Task
// Maximum 0-1 relative error on all pixels with 3 48 : T 0.112457 ; Q 0.453338 ; U 0.453338 
// Average 0-1 relative error on all pixels with 3 48 : T 0.059915 ; Q 0.190214 ; U 0.190214 // With inv_cov = 1, 1 MPI_Task
// Maximum 0-1 relative error on all pixels with 3 48 : T 0.112457 ; Q 0.453338 ; U 0.453338 
// Average 0-1 relative error on all pixels with 3 48 : T 0.059915 ; Q 0.190214 ; U 0.178022 // With inv_cov = 1, 2 MPI_Task
// Maximum 0-1 relative error on all pixels with 3 48 : T 0.112457 ; Q 0.453338 ; U 0.481687 
// NSIDE=512
// Good one with, for 1 MPI_task :
// Average 0-1 relative error on all pixels with 3 3145728 : T 0.000510 ; Q 0.000166 ; U 0.000156 
// Maximum 0-1 relative error on all pixels with 3 3145728 : T 99.646403 ; Q 15.353078 ; U 54.5316
// For more :
// Average 0-1 relative error on all pixels with 3 3145728 : T 0.000075 ; Q 0.000166 ; U 0.000156 
// Maximum 0-1 relative error on all pixels with 3 3145728 : T 18.512676 ; Q 15.353078 ; U 54.531675 

void apply_red_matrix_x_alm(int nside, int lmax, int nstokes, double *CMB_map, double *CMB_map_output, double **red_matrix, double *mask_binary, MPI_Comm worldcomm)
{
    // MPI_Init( &argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &nprocs);
    // printf("Test -1 !!! \n"); fflush(stdout);
    int npix = 12*nside*nside;
    int index, index_2;

    char *c_ell_path = "";
    int number_correlations = ceil((nstokes*nstokes)/2) + floor(nstokes/2) + nstokes%2;
    if (nstokes == 1)
        number_correlations = 1;
    // printf("Test -1a !!! nb_correl %d \n", number_correlations); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    // printf("Test -1ab !!! \n"); fflush(stdout);
    // printf("Initializing S2HAT_param \n"); fflush(stdout);
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);
    // printf("Finish initializing S2HAT_param !! \n"); fflush(stdout);
    // printf("Test -1bb !!! \n"); fflush(stdout);

    S2HAT_LOCAL_parameters Local_param_s2hat = S2HAT_params.Local_param_s2hat;
    // printf("Test -1b !!! \n"); fflush(stdout);
    if (Local_param_s2hat.gangrank >= 0){
        // printf("Entering Wiener-filter with rank %d \n", rank); fflush(stdout);
        
        double *CMB_map_temp  = (double *)malloc(3*npix*sizeof(double));
        // double CMB_map_temp[3*npix];
        memcpy(CMB_map_temp, CMB_map, 3*npix*sizeof(double));
        // printf("Test 0 !!! \n"); fflush(stdout);

        // double *local_map_pix = (double *) calloc( 3*Local_param_s2hat.map_size, sizeof(double));
        // MPI_Barrier(Local_param_s2hat.gangcomm);
        double local_map_pix[3*Local_param_s2hat.map_size];
        distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_temp, local_map_pix, &S2HAT_params);
        // printf("Test 1a !!! \n"); fflush(stdout);
        free(CMB_map_temp);
        // printf("Test 1b !!! \n"); fflush(stdout);
        
        // s2hat_dcomplex *local_alm, *local_alm_out;
        // local_alm = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        // local_alm_out = (s2hat_dcomplex *) malloc( (3*S2HAT_params.size_alm)*sizeof(s2hat_dcomplex));
        // s2hat_dcomplex *local_alm_inverted = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        s2hat_dcomplex local_alm[3*S2HAT_params.size_alm];
        s2hat_dcomplex local_alm_modified[3*S2HAT_params.size_alm];
        
        for (index=0; index<3*S2HAT_params.size_alm; index++){
            local_alm[index].re = 0;
            local_alm[index].im = 0;
            local_alm_modified[index].re = 0;
            local_alm_modified[index].im = 0;
        }
        // printf("Test 1c !!! \n"); fflush(stdout);
        // Doing Pix2Alm
        apply_pix2alm(local_map_pix, local_alm, &S2HAT_params);
        // free(local_map_pix);
        // printf("Test 1d !!! \n"); fflush(stdout);

        // Getting covariance matrix
        int ell_value;
        double **red_matrix_copy;
        red_matrix_copy = malloc((lmax+1)*sizeof(double *));

        // printf("Test 1e !!! \n"); fflush(stdout);
        for(ell_value=0; ell_value<lmax+1; ell_value++){
            red_matrix_copy[ell_value] = calloc(nstokes*nstokes,sizeof(double));
            // printf("- %d -", ell_value); fflush(stdout);
            memcpy(red_matrix_copy[ell_value], red_matrix+ell_value*nstokes*nstokes, nstokes*nstokes*sizeof(double));
            
        }
        // printf("Test 1f !!! \n"); fflush(stdout);

        // printf("Start applying red matrix \n"); fflush(stdout);
        apply_inv_block_diag_covariance_matrix_to_alm(local_alm, local_alm_modified, red_matrix_copy, &S2HAT_params);
        // printf("End applying red matrix \n"); fflush(stdout);

        free_covariance_matrix(red_matrix_copy, lmax);
        // printf("End freeing red matrix \n"); fflush(stdout);

        // double *local_map_pix_out = (double *)calloc(3*Local_param_s2hat.map_size, sizeof(double));
        double local_map_pix_out[3*Local_param_s2hat.map_size];
        
        // apply_alm2pix(local_map_pix_out, local_alm, &S2HAT_params);
        apply_alm2pix(local_map_pix_out, local_alm_modified, &S2HAT_params);
        // free(local_alm);
        // free(local_alm_out);

        double *full_sky_map_2;
        full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
        gather_map(local_map_pix_out, full_sky_map_2, nstokes, &S2HAT_params);
        // free(local_map_pix_out);
        // MPI_Barrier(Local_param_s2hat.gangcomm);

        if (rank==0)
            memcpy(CMB_map_output, full_sky_map_2, nstokes*npix*sizeof(double));


        free(full_sky_map_2);
        // free(local_alm_modified);

        free_s2hat_parameters_struct(&S2HAT_params);

    // MPI_Finalize();
    }
    // printf("Wiener-filter pixel done for rank %d ! \n", rank);
    fflush(stdout);
}
