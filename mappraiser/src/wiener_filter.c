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

    // char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits"; //// TO PUT !!!!
    char *c_ell_path = ""; //// TO PUT !!!!
    
    printf("NEW SIM ################################################################################################################################################ \n");
    printf("Initializing MPI %d %d \n", rank, nprocs); fflush(stdout);

    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);

    printf("Initializing S2HAT_param \n"); fflush(stdout);
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);
    printf("Finish initializing S2HAT_param !! \n"); fflush(stdout);

    S2HAT_LOCAL_parameters Local_param_s2hat = S2HAT_params.Local_param_s2hat;

    if (Local_param_s2hat.gangrank >= 0){
        
        printf("Entering Wiener-filter with rank %d \n", rank); fflush(stdout);
        
        double *CMB_map_temp  = (double *)malloc(3*npix*sizeof(double));
        // double CMB_map_temp[3*npix];
        memcpy(CMB_map_temp, CMB_map, 3*npix*sizeof(double));

        // double *local_map_pix = (double *) calloc( 3*Local_param_s2hat.map_size, sizeof(double));
        // MPI_Barrier(Local_param_s2hat.gangcomm);
        double local_map_pix[3*Local_param_s2hat.map_size];
        distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_temp, local_map_pix, &S2HAT_params);
        free(CMB_map_temp);
        // MPI_Barrier(Local_param_s2hat.gangcomm);
        
        
        s2hat_dcomplex *local_alm, *local_alm_out;
        local_alm = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        // local_alm_out = (s2hat_dcomplex *) malloc( (3*S2HAT_params.size_alm)*sizeof(s2hat_dcomplex));
        s2hat_dcomplex *local_alm_inverted = (s2hat_dcomplex *) calloc( (3*S2HAT_params.size_alm),sizeof(s2hat_dcomplex));
        // s2hat_dcomplex local_alm[3*S2HAT_params.size_alm];
        // s2hat_dcomplex local_alm_inverted[3*S2HAT_params.size_alm];

        int index;
        for (index=0; index<3*S2HAT_params.size_alm; index++){
            local_alm[index].re = 0;
            local_alm[index].im = 0;
            // local_alm_inverted[index].re = 0;
            // local_alm_inverted[index].im = 0;
        }
        // Doing Pix2Alm
        MPI_Barrier(Local_param_s2hat.gangcomm);
        apply_pix2alm(local_map_pix, local_alm, &S2HAT_params);
        // free(local_map_pix);
        MPI_Barrier(Local_param_s2hat.gangcomm);
        // Getting covariance matrix
        int ell_value;
        double **inverse_covariance_matrix;
        inverse_covariance_matrix = malloc((lmax+1)*sizeof(double *));
        double **inverse_covariance_matrix_2;
        inverse_covariance_matrix_2 = malloc((lmax+1)*sizeof(double *));
        
        for(ell_value=0; ell_value<lmax+1; ell_value++){
            inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
            inverse_covariance_matrix_2[ell_value] = calloc(nstokes*nstokes,sizeof(double));
            inverse_covariance_matrix_2[ell_value][0] = 1;
            inverse_covariance_matrix_2[ell_value][4] = 1;
            inverse_covariance_matrix_2[ell_value][8] = 1;
        }
        get_covariance_matrix_block_diagonal(c_ells, number_correlations, inverse_covariance_matrix, &S2HAT_params);
        // inverse_covariance_matrix[lmax][8] = 1;

        int index_2, size_ell_max=10;
        int counter_diff=0;
        for(ell_value=0; ell_value<lmax+1; ell_value++){
            for (index=0; index<nstokes*nstokes; index++){
                if (inverse_covariance_matrix[ell_value][index] - inverse_covariance_matrix_2[ell_value][index] != 0){
                    // printf("- %d %d %f %f -", ell_value, index, inverse_covariance_matrix[ell_value][index], inverse_covariance_matrix_2[ell_value][index]);
                    counter_diff++;
                    }
            }
        }
        printf("Count : %d \n", counter_diff);
        
        // for(ell_value=0; ell_value<lmax; ell_value++){
        // for(ell_value=0; ell_value<size_ell_max; ell_value++){
        //     if (ell_value<size_ell_max){
        //     // if (ell_value<lmax){
        //         printf("\n #### cov mat = %d \n", ell_value);
        //         for (index=0; index < nstokes; index++){
        //             for (index_2=0; index_2<nstokes; index_2++){
        //                 printf("%.7f \t", inverse_covariance_matrix[ell_value][index*nstokes + index_2]);
        //             }
        //             printf("\n");
        //         }
        //     }
        // }

        apply_inv_block_diag_covariance_matrix_to_alm(local_alm, local_alm_inverted, inverse_covariance_matrix_2, &S2HAT_params);

        int count=0;
        for (index=0;index<3*S2HAT_params.size_alm;index++){
            if (fabs(local_alm[index].re - local_alm_inverted[index].re) > 0.005){
                count++;
                printf("- re %d -", index);
            }
            if (fabs(local_alm[index].im - local_alm_inverted[index].im) > 0.005){
                count++;
                printf("- im %d -", index);
            }
        }
        printf("\n");

        printf("Count alms : %d / %d \n", count, nstokes*S2HAT_params.size_alm); fflush(stdout);
        // int max_size_test = 20;
        // printf("%d --- Local_alm diff - %f %f | %f %f -", rank, local_alm[0].re, local_alm[0].im, local_alm_inverted[0].re, local_alm_inverted[0].im);
        // for (index=1;index<max_size_test;index++){
        //     printf("- %f %f | %f %f -", local_alm[index].re, local_alm[index].im, local_alm_inverted[index].re, local_alm_inverted[index].im);
        // }
        // printf(" \n");
        // free_covariance_matrix(inverse_covariance_matrix, lmax+1);
        // free_covariance_matrix(inverse_covariance_matrix_2, lmax+1);



        // double *local_map_pix_out = (double *)calloc(3*Local_param_s2hat.map_size, sizeof(double));
        double local_map_pix_out[3*Local_param_s2hat.map_size];
        MPI_Barrier(Local_param_s2hat.gangcomm);
        // apply_alm2pix(local_map_pix_out, local_alm, &S2HAT_params);
        apply_alm2pix(local_map_pix_out, local_alm_inverted, &S2HAT_params);
        // free(local_alm);
        // free(local_alm_out);
        

        double *full_sky_map_2;
        full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
        MPI_Barrier(Local_param_s2hat.gangcomm);
        gather_map(local_map_pix_out, full_sky_map_2, nstokes, &S2HAT_params);
        // free(local_map_pix_out);
        // MPI_Barrier(Local_param_s2hat.gangcomm);

        if (rank==0)
            memcpy(CMB_map_output, full_sky_map_2, nstokes*npix*sizeof(double));

        if (rank == 0){
            printf("//////////////////////// Comparison new map vs old \n");
            fflush(stdout);

            double average_relative_error_T_0_1 = 0;
            double average_relative_error_Q_0_1 = 0;
            double average_relative_error_U_0_1 = 0;
            double max_average_relative_error_T_0_1 = 0;
            double max_average_relative_error_Q_0_1 = 0;
            double max_average_relative_error_U_0_1 = 0;

            int size_total_probed = npix;

            int init_index = 0;
            int number_0 = 0;

            double pixel_value_T;
            double pixel_value_Q;
            double pixel_value_U;

            int index = init_index;

            max_average_relative_error_T_0_1 = fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
            max_average_relative_error_Q_0_1 = fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
            max_average_relative_error_U_0_1 = fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);
            for(index=init_index;index<size_total_probed+init_index;index++){
                pixel_value_T = fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
                average_relative_error_T_0_1 += pixel_value_T;
                if (max_average_relative_error_T_0_1 < pixel_value_T)
                    max_average_relative_error_T_0_1 = pixel_value_T;
                if(CMB_map[index] == 0.000000){
                    number_0++;
                }
                if(nstokes>1){
                    pixel_value_Q = fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
                    average_relative_error_Q_0_1 += pixel_value_Q;
                    if (max_average_relative_error_Q_0_1 < pixel_value_Q)
                        max_average_relative_error_Q_0_1 = pixel_value_Q;
                    if(CMB_map[index+npix] == 0.000000){
                        number_0++;
                    }
                    pixel_value_U = fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);
                    average_relative_error_U_0_1 += pixel_value_U;
                    if (max_average_relative_error_U_0_1 < pixel_value_U)
                        max_average_relative_error_U_0_1 = pixel_value_U;
                    if(CMB_map[index+2*npix] == 0.000000){
                        number_0++;
                    }
                }
            }
            printf("\n");
            printf("Average 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
                    average_relative_error_T_0_1/(size_total_probed), average_relative_error_Q_0_1/(size_total_probed), average_relative_error_U_0_1/(size_total_probed));
            printf("Maximum 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
                    max_average_relative_error_T_0_1, max_average_relative_error_Q_0_1, max_average_relative_error_U_0_1);
            printf("Number of 0 : %d \n", number_0);
            fflush(stdout);
        }


        free(full_sky_map_2);
        free(local_alm_inverted);
        free_covariance_matrix(inverse_covariance_matrix, lmax);
        free_covariance_matrix(inverse_covariance_matrix_2, lmax);
        free_s2hat_parameters_struct(&S2HAT_params);

    // MPI_Finalize();
    }
    printf("Wiener-filter pixel done for rank %d ! \n", rank);
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
// Average 0-1 relative error on all pixels with 3 48 : T 0.341294 ; Q 0.190214 ; U 0.190214 // With inv_cov = 1, 1 MPI_Task
// Maximum 0-1 relative error on all pixels with 3 48 : T 0.865423 ; Q 0.453338 ; U 0.453338 
// Average 0-1 relative error on all pixels with 3 48 : T 0.059915 ; Q 0.190214 ; U 0.178022 // With inv_cov = 1, 2 MPI_Task
// Maximum 0-1 relative error on all pixels with 3 48 : T 0.112457 ; Q 0.453338 ; U 0.481687 
// NSIDE=512
// Good one with, for 1 MPI_task :
// Average 0-1 relative error on all pixels with 3 3145728 : T 0.000510 ; Q 0.000166 ; U 0.000156 
// Maximum 0-1 relative error on all pixels with 3 3145728 : T 99.646403 ; Q 15.353078 ; U 54.5316
// For more :
// Average 0-1 relative error on all pixels with 3 3145728 : T 0.000075 ; Q 0.000166 ; U 0.000156 
// Maximum 0-1 relative error on all pixels with 3 3145728 : T 18.512676 ; Q 15.353078 ; U 54.531675 
// void apply_Wiener_filter_pixel_good_one(int nside_, int lmax_, int nstokes_, double *CMB_map_ring_, double *CMB_map_output_, double *c_ells_, int number_correlations_, double *mask_binary_, MPI_Comm worldcomm)
// {

//     char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_band_limited_1024_0.fits";
//     // char *path_CMB_map = "/Users/mag/Documents/PHD1Y/Space_Work/Inference_Sampling/map_files/Map_band_limited_1024_0.fits";
//     double *CMB_map;

//     int f_sky, npix;
//     int i_index, ncomp;
//     int index;

//     int nside = 512;
//     int lmax = 2*nside; //+2; //+10; //3*nside-1 ;//1500;//025;
//     npix = 12*nside*nside;

//     char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits"; //// TO PUT !!!!
//     int number_correlations = 4; // 6; // TO MODIFY LATER !!!!!!!
//     int nstokes = 3;
//     // int number_correlations = 1; // 6; // TO MODIFY LATER !!!!!!!
//     // int nstokes = 1;    

//     int rank, nprocs;
//     MPI_Comm gangcomm;
    
//     // MPI_Init( &argc, &argv);
//     MPI_Comm_rank( worldcomm, &rank);
//     MPI_Comm_size( worldcomm, &nprocs);

//     printf("NEW SIM ################################################################################################################################################ \n");
//     printf("Initializing MPI %d %d \n", rank, nprocs);
//     fflush(stdout);

//     int gangrank = rank;
//     int gangsize = nprocs;
//     gangcomm = MPI_COMM_WORLD;


//     printf("Initializing Global_param_s2hat \n");

//     S2HAT_parameters S2HAT_params;
//     init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
//     printf("--- Test init %d %d \n", (S2HAT_params.Files_WF_struct).lmax_Wiener_Filter, lmax); fflush(stdout);
//     int i;

//     double *mask_binary=NULL;// = calloc(12*nside*nside, sizeof(int));

//     printf("Initializing S2HAT_param \n"); fflush(stdout);
//     init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, gangcomm);
//     printf("Finish initializing S2HAT_param !! \n"); fflush(stdout);

//     S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
//     S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);
//     int first_ring = Local_param_s2hat->first_ring;
//     // printf("Test verif3 : %d %d \n", first_ring, Global_param_s2hat->pixelization_scheme.fpix[first_ring-1]); fflush(stdout);
//     // printf("###### Test4 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     // printf("###### Test12 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     printf("Getting CMB maps !!!\n");
//     fflush(stdout);
//     // printf("###### Test10 - %d \n", Local_param_s2shat->pixel_numbered_ring[0]);
//     npix = 12*nside*nside;
//     printf("###### Test nside - %d %d \n", nside, npix);
//     fflush(stdout);
//     CMB_map = (double *) malloc( 3*npix*sizeof(double));
//     printf("###### ??? \n");
//     fflush(stdout);
//     printf("###### Test13 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     fflush(stdout);
//     read_TQU_maps(nside, CMB_map, path_CMB_map, 3);
//     printf("###### Test14 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     // printf("Reading map - rank %d \n", gangrank);
//     // fflush(stdout);
//     printf("CMB map read !!!\n");
//     fflush(stdout);
//     printf("###### Test11 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     // double *CMB_map_temp  = CMB_map+0; // &(CMB_map[npix]);
//     double *CMB_map_temp  = (double *)malloc(3*npix*sizeof(double));
//     memcpy(CMB_map_temp, CMB_map, 3*npix*sizeof(double));

//     // printf("###### Test9 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     // printf("Minute test1 - rank %d - %f %f %f \n", gangrank, CMB_map[0], CMB_map[npix], CMB_map[2*npix]); fflush(stdout);
//     // printf("Minute test1_bis - rank %d - %d \n", gangrank, Global_param_s2hat->nlmax); fflush(stdout);


//     double *local_map_pix;//, *local_map_pix_E, *local_map_pix_B;
//     local_map_pix = (double *) calloc( 3*Local_param_s2hat->map_size, sizeof(double));
//     // printf("###### Test8 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]); fflush(stdout);
//     // printf("###### TEST2 --- MPI struct - %d \n", Local_param_s2hat->gangcomm); fflush(stdout);
//     // printf("###### CMB Temp - %f \n", CMB_map_temp[0]); fflush(stdout);

//     distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_temp, local_map_pix, &S2HAT_params);
    
//     // printf("###### Test5 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     // printf("Tiny test - rank %d - %f \n", gangrank, CMB_map_temp[0]); fflush(stdout);

//     s2hat_dcomplex *local_alm;
    
//     local_alm = (s2hat_dcomplex *) calloc( ((3*(Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax)),sizeof(s2hat_dcomplex));

//     // printf("## Test dims local_alm : %d, %d, %d \n", Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, nstokes);

//     // printf("Little test before pix2alm - rank %d - %f %f \n", gangrank, local_map_pix[0], CMB_map_temp[0]);
//     // printf("Little test before pix2alm - rank %d - %f %f \n", gangrank, local_map_pix[0], CMB_map[0]);
//     fflush(stdout);


//     int max_size_test = 20;

//     // printf("Map_pix after distribute_map - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
//     // printf("Map_pix after distribute_map - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map_temp[0]);
//     // for (index=1;index<max_size_test;index++){
//     //     printf("- %f %f -", local_map_pix[index], CMB_map_temp[index]);
//     // }
//     // printf(" \n");
//     // printf("Map_pix -2- after distribute_map - rank %d - %f %f -", gangrank, local_map_pix[Local_param_s2hat->map_size], CMB_map_temp[npix]);
//     // for (index=Local_param_s2hat->map_size;index<max_size_test+Local_param_s2hat->map_size;index++){
//     //     printf("- %f %f -", local_map_pix[index], CMB_map_temp[index]);
//     // }
//     // printf(" \n");
//     // fflush(stdout);
//     // printf("###### Test6 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     int nmaps = 1; // We only provide 1 input set of alm coefficient
//     // int lda = nstokes;
//     int lda = Global_param_s2hat->nlmax;
    
//     double *local_w8ring;
//     int i_ring;

//     int spin;

//     // printf("----- Test nstokes : %d \n", S2HAT_params.nstokes);
//     apply_pix2alm(local_map_pix, local_alm, &S2HAT_params);

//     // printf("###### Test7 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     // printf("Test1 between the 2 pix2alm - rank %d \n", gangrank);
//     // fflush(stdout);
//     // printf("Test2 between the 2 pix2alm - rank %d \n", gangrank);
//     // fflush(stdout);

//     // printf("###### Test3 - %d \n", Local_param_s2hat->pixel_numbered_ring[0]);
//     // printf("Map_pix after first apply_pix2alm - rank %d - %f %f - \n", gangrank, local_map_pix[0], CMB_map_temp[0]); fflush(stdout);


//     // max_size_test = 50;
    
//     // printf("%d --- Local_alm long after 1st pix2alm - %f %f -", rank, local_alm[0].re, local_alm[0].im);
//     // for (index=1;index<max_size_test;index++){
//     //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
//     // }
//     // printf(" \n");
//     // printf("%d --- Local_alm long after 1st pix2alm + lmax*lmax/2 - %f %f -", rank, local_alm[lmax*lmax/20].re, local_alm[lmax*lmax/20].im);
//     // for (index=1;index<max_size_test;index++){
//     //     printf("- %f %f -", local_alm[index+lmax*lmax/20].re, local_alm[index+lmax*lmax/20].im);
//     // }
//     // printf(" \n");

//     // printf("Tset 0 -- %d \n", Local_param_s2hat->map_size); fflush(stdout);
//     printf("Tset 0F \n"); fflush(stdout);
//     double *local_map_pix_out = (double *)calloc(3*Local_param_s2hat->map_size, sizeof(double));
//     apply_alm2pix(local_map_pix_out, local_alm, &S2HAT_params);
//     printf("Map_pix after first apply_alm2pix - rank %d - %f %f - \n", gangrank, local_map_pix_out[0], CMB_map_temp[0]); fflush(stdout);

//     fflush(stdout);

//     double *full_sky_map_2;
//     full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
//     gather_map(local_map_pix_out, full_sky_map_2, nstokes, &S2HAT_params);
//     // memcpy(full_sky_map_2, local_map_pix, nstokes*npix*sizeof(double));

//     if (rank==0)
//         memcpy(CMB_map_output_, full_sky_map_2, nstokes*npix*sizeof(double));

//     if (Local_param_s2hat->gangrank == 0){
//     printf("//////////////////////// Comparison new map vs old \n");
//     fflush(stdout);

//     printf("From pixel 0 to %d - %f %f -", max_size_test, full_sky_map_2[0], CMB_map_temp[0]);
//     for (index=0;index<max_size_test;index++){
//         printf("- %f %f -", full_sky_map_2[index], CMB_map_temp[index]);
//         }
//     printf(" \n");
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

//     index = init_index;
    
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
//     }

//     free(CMB_map);
//     free(CMB_map_temp);

    
//     // printf("Test 9 !\n"); fflush(stdout);
//     free_s2hat_parameters_struct(&S2HAT_params);

//     printf("Test 10 !\n"); fflush(stdout);
//     // MPI_Finalize();

//     printf("Test finish ! - %d \n", gangrank);
//     fflush(stdout);

// }

