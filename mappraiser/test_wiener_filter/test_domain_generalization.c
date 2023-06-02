
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

// #include "mappraiser.h"
#include "domain_generalization.h"


// /* Initalize PCG_var structure */
// int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix, int domain_PCG_computation, int bool_apply_filter);

// /* Initialize harmonic superstructure */
// int init_harmonic_superstruct(Mat *A, Harmonic_superstruct *Harm_struct, int *mask_binary);

// /* Transforms local map pixels into local alm in harmonic domain and vice-versa*/
// int global_map_2_harmonic(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, Harmonic_superstruct *Harmonic_sup);
// int global_harmonic_2_map(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A,  Harmonic_superstruct *Harmonic_sup);

// int get_mask_from_indices(Mat *A, int *mask_binary, int nside, int root);

// int free_harmonic_superstruct(Harmonic_superstruct *Harmonic_sup, int rank);
// int free_PCG_var(PCG_var *PCG_var_obj);


int main_firsts_tests_init(int argc, char** argv){
// int main(int argc, char** argv){

    char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_band_limited_1024_0.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits";
    int rank, nprocs;
    int i;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int nside = 512;
    // int nside = 4;
    int lmax = 2*nside;
    int domain_PCG_computation = 0;
    int bool_apply_filter = 0;
    int nstokes = 3;
    int number_correlations = 4;


    int npix = 12*nside*nside;
    int index = 0;
    int max_size_test = 10;

    double *CMB_map_ring = (double *) malloc( nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, nstokes);
    // printf("%d --- Map read !!! %f %f %f\n", rank, CMB_map_ring[0], CMB_map_ring[1570048], CMB_map_ring[npix+1570048]); fflush(stdout);
    // printf("%d ---- Initial ring CMB map : From pixel 0 to %d - %f -", rank, max_size_test, CMB_map_ring[index]);
    // for (index=0;index<max_size_test;index++){
    //     printf("- %f -", CMB_map_ring[index]);
    // }
    
    printf(" \n"); fflush(stdout);
    double *CMB_map_nest = (double *) malloc( nstokes*npix*sizeof(double));
    convert_full_map_ring2nest(CMB_map_ring, CMB_map_nest, nside, nstokes);
    printf("%d --- Map transformed !!! %f %f \n", rank, CMB_map_nest[0], CMB_map_nest[1]); fflush(stdout);
    free(CMB_map_ring);


    int overlap = 0; // 10;
    Mat Obj_pointing_matrix;
    Obj_pointing_matrix.trash_pix = 0;
    Obj_pointing_matrix.nnz = nstokes;
    Obj_pointing_matrix.lcount = (npix/nprocs + overlap)*nstokes;
    Obj_pointing_matrix.comm = worldcomm;

    if (rank == nprocs-1)
        Obj_pointing_matrix.lcount += (npix%nprocs - overlap)*nstokes;
    printf("%d --- Producing lcount !!! %d over %d\n", rank, Obj_pointing_matrix.lcount, nstokes*npix); fflush(stdout);

    int number_pixels_MAPP = Obj_pointing_matrix.lcount-(Obj_pointing_matrix.nnz)*(Obj_pointing_matrix.trash_pix);
    
    Obj_pointing_matrix.lindices = (int *)malloc(number_pixels_MAPP*sizeof(int));

    int first_pixel = rank*(npix/nprocs)*nstokes;
    for (i=0; i<Obj_pointing_matrix.lcount; i++){
        Obj_pointing_matrix.lindices[i] = first_pixel + i;
    }
    printf("%d --- Producing lindices !!! first_pixel %d -- 0: %d last: %d \n", rank, first_pixel, Obj_pointing_matrix.lindices[0], Obj_pointing_matrix.lindices[Obj_pointing_matrix.lcount-1]); 
    fflush(stdout);

    double *local_map_MAPPRAISER, *local_map_MAPPRAISER_output;
    s2hat_dcomplex *local_alm_s2hat;
    
    local_map_MAPPRAISER = (double *)malloc(Obj_pointing_matrix.lcount*sizeof(double));

    printf("%d -- Getting min, max of lindices \n", rank);
    int min_val = Obj_pointing_matrix.lindices[0];
    int max_val = Obj_pointing_matrix.lindices[0];
    for (i=1; i<Obj_pointing_matrix.lcount; i++){
        if (min_val > Obj_pointing_matrix.lindices[i])
            min_val = Obj_pointing_matrix.lindices[i];
        if (max_val < Obj_pointing_matrix.lindices[i])
            max_val = Obj_pointing_matrix.lindices[i];
    }
    printf("%d -- Min indice : %d ; Max indice : %d \n", rank, min_val, max_val); fflush(stdout);
    printf("%d -- Test last value CMB_map_nest %f \n", rank, CMB_map_nest[nstokes*npix - 1]); fflush(stdout);
    printf("%d --- Test before setting local_map_MAPPRAISER !!! \n", rank); fflush(stdout);
    for (i=0; i<Obj_pointing_matrix.lcount; i++){
        local_map_MAPPRAISER[i] = CMB_map_nest[Obj_pointing_matrix.lindices[i]];
    }
    printf("%d --- Producing local_map_MAPPRAISER !!! first %f last: %f \n", rank, local_map_MAPPRAISER[0], local_map_MAPPRAISER[Obj_pointing_matrix.lcount-1]); fflush(stdout);
    printf("%d --- From CMB_map_nest !!! first %f last: %f \n", rank, CMB_map_nest[Obj_pointing_matrix.lindices[0]], CMB_map_nest[Obj_pointing_matrix.lindices[Obj_pointing_matrix.lcount-1]]); 
    fflush(stdout);


    // PCG_var PCG_variable;
    // initialize_PCG_var_struct(&(PCG_variable), local_map_MAPPRAISER);

    int *mask_binary = NULL;
    Harmonic_superstruct Harmonic_struct;
    printf("%d --- initializing harmonic superstruct !!! \n", rank); fflush(stdout);
    init_harmonic_superstruct(&Obj_pointing_matrix, &(Harmonic_struct), mask_binary, nside, lmax, c_ell_path, number_correlations);
    printf("%d --- finishing producing harmonic superstruct !!!!!!! \n", rank); fflush(stdout);
    S2HAT_parameters *S2HAT_params = &(Harmonic_struct.S2HAT_params);

    if (S2HAT_params->Local_param_s2hat.gangrank >= 0)
        local_alm_s2hat = (s2hat_dcomplex *)malloc(S2HAT_params->nstokes*S2HAT_params->size_alm*sizeof(s2hat_dcomplex));

    printf("%d --- map2harmonic !!! \n", rank); fflush(stdout);
    global_map_2_harmonic(local_map_MAPPRAISER, local_alm_s2hat, &Obj_pointing_matrix, &(Harmonic_struct));

    max_size_test = 50;
    printf("%d ---- ***alms*** from pixel 0 to %d - %f %f -", rank, max_size_test, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    int number_of_nan = 0;
    for (index=0;index<S2HAT_params->nstokes*S2HAT_params->size_alm;index++){
        printf("- %f %f -", local_alm_s2hat[index].re, local_alm_s2hat[index].im);
        if (!(local_alm_s2hat[index].re == local_alm_s2hat[index].re)){
            local_alm_s2hat[index].re = 0;
            number_of_nan++;
            }
        if (!(local_alm_s2hat[index].im == local_alm_s2hat[index].im)){
            local_alm_s2hat[index].im = 0;
            number_of_nan++;
            }
    }
    printf(" \n"); fflush(stdout);

    printf("%d ----!!!! ***alms*** Number of nans %d \n", rank, number_of_nan); fflush(stdout);
    // double *local_map_pix = (double *)malloc(nstokes*S2HAT_params->Local_param_s2hat.map_size*sizeof(double));
    // apply_alm2pix(local_map_pix, local_alm_s2hat, S2HAT_params);

    //  max_size_test = 10;
    // printf("%d ---- INTERMEDIARY STEP alm2pix IN S2HAT : From pixel 0 to %d - %f -", rank, max_size_test, local_map_pix[index]);
    // for (index=0;index<max_size_test;index++){
    //     printf("- %f -", local_map_pix[index]);
    //     }
    // printf(" \n"); fflush(stdout);
    // free(local_map_pix);

    local_map_MAPPRAISER_output = (double *)malloc(number_pixels_MAPP*sizeof(double));
    printf("%d --- harmonic2map !!! \n", rank); fflush(stdout);
    global_harmonic_2_map(local_map_MAPPRAISER_output, local_alm_s2hat, &Obj_pointing_matrix, &(Harmonic_struct));

    printf("%d ---- //////////////////////// Comparison new map vs old \n", rank); fflush(stdout);

    max_size_test = 10;
    printf("%d ---- From pixel 0 to %d - %f %f -", rank, max_size_test, local_map_MAPPRAISER_output[0], local_map_MAPPRAISER[0]);
    for (index=0;index<max_size_test;index++){
        printf("- %f %f -", local_map_MAPPRAISER_output[index*nstokes], local_map_MAPPRAISER[index*nstokes]);
        }
    printf(" \n");
    fflush(stdout);

    // int new_pixel = npix;
    int new_pixel = 0;
    printf("%d ---2- From pixel %d to %d - %f %f -", rank, new_pixel, new_pixel+max_size_test, local_map_MAPPRAISER_output[new_pixel*nstokes+1], local_map_MAPPRAISER[new_pixel*nstokes+1]);
    for (index=new_pixel;index<new_pixel+max_size_test;index++){
        printf("- %f %f -", local_map_MAPPRAISER_output[index*nstokes+1], local_map_MAPPRAISER[index*nstokes+1]);
            }
    printf(" \n");
    fflush(stdout);

    // new_pixel = 2*npix;
    printf("%d ---3- From pixel %d to %d - %f %f -", rank, new_pixel, new_pixel+max_size_test, local_map_MAPPRAISER_output[new_pixel*nstokes+2], local_map_MAPPRAISER[new_pixel*nstokes+2]);
    for (index=new_pixel;index<new_pixel+max_size_test;index++){
        printf("- %f %f -", local_map_MAPPRAISER_output[index*nstokes+2], local_map_MAPPRAISER[index*nstokes+2]);
        }
    printf(" \n");

    double average_relative_error_T_0_1 = 0;
    double average_relative_error_Q_0_1 = 0;
    double average_relative_error_U_0_1 = 0;
    double max_average_relative_error_T_0_1 = 0;
    double max_average_relative_error_Q_0_1 = 0;
    double max_average_relative_error_U_0_1 = 0;
    double average_relative_error_0_2 = 0;
    double average_relative_error_0_2b = 0;
    double max_average_relative_error_0_2 = 0;

    int size_total_probed = Obj_pointing_matrix.lcount/nstokes;
    // size_total_probed = 10;
    int index_complete_map;
    int init_index = first_pixel;
    int index_nest = Obj_pointing_matrix.lindices[i];
    int number_0 = 0;
    int index_U_max = 0;

    index = init_index;
    printf("%d --- calculating statistics!!! \n", rank); fflush(stdout);
    max_average_relative_error_T_0_1 = fabs((local_map_MAPPRAISER_output[0]-CMB_map_nest[0])/CMB_map_nest[index]);
    max_average_relative_error_Q_0_1 = fabs((local_map_MAPPRAISER_output[0+1]-CMB_map_nest[index+1])/CMB_map_nest[index+1]);
    max_average_relative_error_U_0_1 = fabs((local_map_MAPPRAISER_output[0+2]-CMB_map_nest[index+2])/CMB_map_nest[index+2]);
    for(index=0;index<size_total_probed;index++){
        index_complete_map = first_pixel + index*nstokes;
        average_relative_error_T_0_1 += fabs((local_map_MAPPRAISER_output[index*nstokes]-CMB_map_nest[index_complete_map])/CMB_map_nest[index_complete_map]);
        if (max_average_relative_error_T_0_1 < fabs((local_map_MAPPRAISER_output[index*nstokes]-CMB_map_nest[index_complete_map])/CMB_map_nest[index_complete_map]));
            max_average_relative_error_T_0_1 = fabs((local_map_MAPPRAISER_output[index*nstokes]-CMB_map_nest[index_complete_map])/CMB_map_nest[index_complete_map]);

        if(CMB_map_nest[index] == 0.000000){
            number_0++;
        }
        if(nstokes>1){
            average_relative_error_Q_0_1 += fabs((local_map_MAPPRAISER_output[index*nstokes+1]-CMB_map_nest[index_complete_map+1])/CMB_map_nest[index_complete_map+1]);
            if (max_average_relative_error_Q_0_1 < fabs((local_map_MAPPRAISER_output[index*nstokes+1]-CMB_map_nest[index_complete_map+1])/CMB_map_nest[index_complete_map+1]));
                max_average_relative_error_Q_0_1 = fabs((local_map_MAPPRAISER_output[index*nstokes+1]-CMB_map_nest[index_complete_map+1])/CMB_map_nest[index_complete_map+1]);


            average_relative_error_U_0_1 += fabs((local_map_MAPPRAISER_output[index*nstokes+2]-CMB_map_nest[index_complete_map+2])/CMB_map_nest[index_complete_map+2]);
            if (max_average_relative_error_U_0_1 < fabs((local_map_MAPPRAISER_output[index*nstokes+2]-CMB_map_nest[index_complete_map+2])/CMB_map_nest[index_complete_map+2])){
                max_average_relative_error_U_0_1 = fabs((local_map_MAPPRAISER_output[index*nstokes+2]-CMB_map_nest[index_complete_map+2])/CMB_map_nest[index_complete_map+2]);
                index_U_max = index;
                }

        }
    }
    printf("\n");
    printf("%d --- Average 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", rank, nstokes, npix, 
            average_relative_error_T_0_1/(size_total_probed), average_relative_error_Q_0_1/(size_total_probed), average_relative_error_U_0_1/(size_total_probed));
    printf("%d --- Maximum 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", rank, nstokes, npix, 
            max_average_relative_error_T_0_1, max_average_relative_error_Q_0_1, max_average_relative_error_U_0_1);
    printf("%d --- Number of 0 : %d ; Pixel for max value U : %d \n", rank, number_0, index_U_max); fflush(stdout);
    
    printf("%d --- calculating statistics!!! \n", rank); fflush(stdout);
    max_average_relative_error_0_2 = fabs((local_map_MAPPRAISER_output[0]-local_map_MAPPRAISER[0])/local_map_MAPPRAISER[0]);

    average_relative_error_0_2b = fabs((local_map_MAPPRAISER_output[0]-local_map_MAPPRAISER[0])/local_map_MAPPRAISER[0]);
    printf("%d --- test statistics %f \n", rank, 
            average_relative_error_0_2b);

    // max_average_relative_error_Q_0_2 = fabs((local_map_MAPPRAISER_output[0]-CMB_map_nest[index+1])/CMB_map_nest[index+1]);
    // max_average_relative_error_U_0_2 = fabs((local_map_MAPPRAISER_output[0+2]-CMB_map_nest[index+2])/CMB_map_nest[index+2]);
    for(index=1;index<Obj_pointing_matrix.lcount;index++){
        average_relative_error_0_2 += fabs((local_map_MAPPRAISER_output[index]-local_map_MAPPRAISER[index])/local_map_MAPPRAISER[index]);
        if (max_average_relative_error_0_2 < fabs((local_map_MAPPRAISER_output[index]-local_map_MAPPRAISER[index])/local_map_MAPPRAISER[index]));
            max_average_relative_error_0_2 = fabs((local_map_MAPPRAISER_output[index]-local_map_MAPPRAISER[index])/local_map_MAPPRAISER[index]);
        }
    
    printf("%d --- Average 2bis relative error on all pixels with %d %d : Total %f \n", rank, nstokes, npix, 
            average_relative_error_0_2/(Obj_pointing_matrix.lcount));
    printf("%d --- Maximum 2bis relative error on all pixels with %d %d : Total %f \n", rank, nstokes, npix, 
            max_average_relative_error_0_2); fflush(stdout);

    // free_PCG_var(PCG_variable);
    printf("%d --- Free step \n", rank); fflush(stdout);
    free_harmonic_superstruct(&Harmonic_struct, rank);
    free(Obj_pointing_matrix.lindices);
    printf("%d --- Free step 2 \n", rank); fflush(stdout);
    if (S2HAT_params->Local_param_s2hat.gangrank >= 0)
        free(local_alm_s2hat);

    free(local_map_MAPPRAISER);
    free(local_map_MAPPRAISER_output);
    free(CMB_map_nest);
    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}


// int main_multiplying_cov_matrix(int argc, char** argv){
int main(int argc, char** argv){

    // char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_band_limited_1024_0.fits";
    char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_band_limited_1024_0.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits";
    int rank, nprocs;
    int i, ell_value;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int nside = 512;
    // int nside = 4;
    int lmax = 2*nside+2;
    int domain_PCG_computation = 0;
    int bool_apply_filter = 0;
    int nstokes = 3;
    int number_correlations = 4;


    int npix = 12*nside*nside;
    int index = 0;
    int max_size_test = 10;

    double *CMB_map_ring = (double *) malloc( nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, nstokes);
    // printf("%d --- Map read !!! %f %f %f\n", rank, CMB_map_ring[0], CMB_map_ring[1570048], CMB_map_ring[npix+1570048]); fflush(stdout);
    // printf("%d ---- Initial ring CMB map : From pixel 0 to %d - %f -", rank, max_size_test, CMB_map_ring[index]);
    // for (index=0;index<max_size_test;index++){
    //     printf("- %f -", CMB_map_ring[index]);
    // }
    
    printf(" \n"); fflush(stdout);
    double *CMB_map_nest = (double *) malloc( nstokes*npix*sizeof(double));
    convert_full_map_ring2nest(CMB_map_ring, CMB_map_nest, nside, nstokes);
    printf("%d --- Map transformed !!! %f %f \n", rank, CMB_map_nest[0], CMB_map_nest[1]); fflush(stdout);
    free(CMB_map_ring);


    int overlap = 0; // 10;
    Mat Obj_pointing_matrix;
    Obj_pointing_matrix.trash_pix = 0;
    Obj_pointing_matrix.nnz = nstokes;
    Obj_pointing_matrix.lcount = (npix/nprocs + overlap)*nstokes;
    Obj_pointing_matrix.comm = worldcomm;

    if (rank == nprocs-1)
        Obj_pointing_matrix.lcount += (npix%nprocs - overlap)*nstokes;
    printf("%d --- Producing lcount !!! %d over %d\n", rank, Obj_pointing_matrix.lcount, nstokes*npix); fflush(stdout);

    int number_pixels_MAPP = Obj_pointing_matrix.lcount-(Obj_pointing_matrix.nnz)*(Obj_pointing_matrix.trash_pix);
    
    Obj_pointing_matrix.lindices = (int *)malloc(number_pixels_MAPP*sizeof(int));

    int first_pixel = rank*(npix/nprocs)*nstokes;
    for (i=0; i<Obj_pointing_matrix.lcount; i++){
        Obj_pointing_matrix.lindices[i] = first_pixel + i;
    }
    printf("%d --- Producing lindices !!! first_pixel %d -- 0: %d last: %d \n", rank, first_pixel, Obj_pointing_matrix.lindices[0], Obj_pointing_matrix.lindices[Obj_pointing_matrix.lcount-1]); 
    fflush(stdout);

    double *local_map_MAPPRAISER, *local_map_MAPPRAISER_output;
    s2hat_dcomplex *local_alm_s2hat;
    
    local_map_MAPPRAISER = (double *)malloc(Obj_pointing_matrix.lcount*sizeof(double));

    printf("%d -- Getting min, max of lindices \n", rank);
    int min_val = Obj_pointing_matrix.lindices[0];
    int max_val = Obj_pointing_matrix.lindices[0];
    for (i=1; i<Obj_pointing_matrix.lcount; i++){
        if (min_val > Obj_pointing_matrix.lindices[i])
            min_val = Obj_pointing_matrix.lindices[i];
        if (max_val < Obj_pointing_matrix.lindices[i])
            max_val = Obj_pointing_matrix.lindices[i];
    }
    printf("%d -- Min indice : %d ; Max indice : %d \n", rank, min_val, max_val); fflush(stdout);
    printf("%d -- Test last value CMB_map_nest %f \n", rank, CMB_map_nest[nstokes*npix - 1]); fflush(stdout);
    printf("%d --- Test before setting local_map_MAPPRAISER !!! \n", rank); fflush(stdout);
    for (i=0; i<Obj_pointing_matrix.lcount; i++){
        local_map_MAPPRAISER[i] = CMB_map_nest[Obj_pointing_matrix.lindices[i]];
    }
    printf("%d --- Producing local_map_MAPPRAISER !!! first %f last: %f \n", rank, local_map_MAPPRAISER[0], local_map_MAPPRAISER[Obj_pointing_matrix.lcount-1]); fflush(stdout);
    printf("%d --- From CMB_map_nest !!! first %f last: %f \n", rank, CMB_map_nest[Obj_pointing_matrix.lindices[0]], CMB_map_nest[Obj_pointing_matrix.lindices[Obj_pointing_matrix.lcount-1]]); 
    fflush(stdout);


    // PCG_var PCG_variable;
    // initialize_PCG_var_struct(&(PCG_variable), local_map_MAPPRAISER);

    int *mask_binary = NULL;
    Harmonic_superstruct Harmonic_struct;
    printf("%d --- initializing harmonic superstruct !!! \n", rank); fflush(stdout);
    init_harmonic_superstruct(&Obj_pointing_matrix, &(Harmonic_struct), mask_binary, nside, lmax, c_ell_path, number_correlations);
    printf("%d --- finishing producing harmonic superstruct !!!!!!! \n", rank); fflush(stdout);
    S2HAT_parameters *S2HAT_params = &(Harmonic_struct.S2HAT_params);

    // double **inverse_covariance_matrix = calloc(lmax, sizeof(double *));
    // for(ell_value=0; ell_value<lmax; ell_value++){
    //     inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double)); 
    //     // (nstokes*nstokes) for each element of the covariance matrix either [TT], [EE, EB, BE, BB], or [TT, TE, TB, ET, EE, EB, BT, BE, BB]
    //     }

    // printf("%d --- Setting inverse covariance matrix to identity \n", rank);
    // for(ell_value=0; ell_value<lmax; ell_value++){
    //     for (index=0; index<nstokes; index++)
    //         inverse_covariance_matrix[ell_value][index*nstokes + index] = 1; 
    //         // (nstokes*nstokes) for each element of the covariance matrix either [TT], [EE, EB, BE, BB], or [TT, TE, TB, ET, EE, EB, BT, BE, BB]
    //     }
    // // get_inverse_covariance_matrix_NxN(&(Harmonic_sup.S2HAT_params), inverse_covariance_matrix);

    if (S2HAT_params->Local_param_s2hat.gangrank >= 0)
        local_alm_s2hat = (s2hat_dcomplex *)malloc(S2HAT_params->nstokes*S2HAT_params->size_alm*sizeof(s2hat_dcomplex));

    printf("%d --- map2harmonic !!! \n", rank); fflush(stdout);
    global_map_2_harmonic(local_map_MAPPRAISER, local_alm_s2hat, &Obj_pointing_matrix, &(Harmonic_struct));
    printf("%d --- map2harmonic end !!! \n", rank); fflush(stdout);

    max_size_test = 50;
    printf("%d ---- ***alms*** from pixel 0 to %d - %f %f -", rank, max_size_test, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    for (index=0;index<max_size_test;index++){
        printf("- %f %f -", local_alm_s2hat[index].re, local_alm_s2hat[index].im);
        }
    printf(" \n");
    fflush(stdout);

    int number_of_nan = 0;
    for (index=0;index<S2HAT_params->nstokes*S2HAT_params->size_alm;index++){
        if (!(local_alm_s2hat[index].re == local_alm_s2hat[index].re)){
            // local_alm_s2hat[index].re = 0;
            number_of_nan++;
            }
        if (!(local_alm_s2hat[index].im == local_alm_s2hat[index].im)){
            // local_alm_s2hat[index].im = 0;
            number_of_nan++;
            }
    }
    printf(" \n"); fflush(stdout);

    printf("%d ----!!!! ***alms*** Number of nans %d \n", rank, number_of_nan); fflush(stdout);
    // double *local_map_pix = (double *)malloc(nstokes*S2HAT_params->Local_param_s2hat.map_size*sizeof(double));
    // apply_alm2pix(local_map_pix, local_alm_s2hat, S2HAT_params);

    //  max_size_test = 10;
    // printf("%d ---- INTERMEDIARY STEP alm2pix IN S2HAT : From pixel 0 to %d - %f -", rank, max_size_test, local_map_pix[index]);
    // for (index=0;index<max_size_test;index++){
    //     printf("- %f -", local_map_pix[index]);
    //     }
    // printf(" \n"); fflush(stdout);
    // free(local_map_pix);

    local_map_MAPPRAISER_output = (double *)malloc(number_pixels_MAPP*sizeof(double));
    printf("%d --- harmonic2map !!! \n", rank); fflush(stdout);
    global_harmonic_2_map(local_map_MAPPRAISER_output, local_alm_s2hat, &Obj_pointing_matrix, &(Harmonic_struct));

    printf("%d ---- //////////////////////// Comparison new map vs old \n", rank); fflush(stdout);

    max_size_test = 10;
    printf("%d ---- From pixel 0 to %d - %f %f -", rank, max_size_test, local_map_MAPPRAISER_output[0], local_map_MAPPRAISER[0]);
    for (index=0;index<max_size_test;index++){
        printf("- %f %f -", local_map_MAPPRAISER_output[index*nstokes], local_map_MAPPRAISER[index*nstokes]);
        }
    printf(" \n");
    fflush(stdout);

    // int new_pixel = npix;
    int new_pixel = 0;
    printf("%d ---2- From pixel %d to %d - %f %f -", rank, new_pixel, new_pixel+max_size_test, local_map_MAPPRAISER_output[new_pixel*nstokes+1], local_map_MAPPRAISER[new_pixel*nstokes+1]);
    for (index=new_pixel;index<new_pixel+max_size_test;index++){
        printf("- %f %f -", local_map_MAPPRAISER_output[index*nstokes+1], local_map_MAPPRAISER[index*nstokes+1]);
            }
    printf(" \n");
    fflush(stdout);

    // new_pixel = 2*npix;
    printf("%d ---3- From pixel %d to %d - %f %f -", rank, new_pixel, new_pixel+max_size_test, local_map_MAPPRAISER_output[new_pixel*nstokes+2], local_map_MAPPRAISER[new_pixel*nstokes+2]);
    for (index=new_pixel;index<new_pixel+max_size_test;index++){
        printf("- %f %f -", local_map_MAPPRAISER_output[index*nstokes+2], local_map_MAPPRAISER[index*nstokes+2]);
        }
    printf(" \n");

    double average_relative_error_T_0_1 = 0;
    double average_relative_error_Q_0_1 = 0;
    double average_relative_error_U_0_1 = 0;
    double max_average_relative_error_T_0_1 = 0;
    double max_average_relative_error_Q_0_1 = 0;
    double max_average_relative_error_U_0_1 = 0;
    double average_relative_error_0_2 = 0;
    double average_relative_error_0_2b = 0;
    double max_average_relative_error_0_2 = 0;

    int size_total_probed = Obj_pointing_matrix.lcount/nstokes;
    // size_total_probed = 10;
    int index_complete_map;
    int init_index = first_pixel;
    int index_nest = Obj_pointing_matrix.lindices[i];
    int number_0 = 0;
    int index_U_max = 0;

    index = init_index;
    printf("%d --- calculating statistics!!! \n", rank); fflush(stdout);
    max_average_relative_error_T_0_1 = fabs((local_map_MAPPRAISER_output[0]-CMB_map_nest[0])/CMB_map_nest[index]);
    max_average_relative_error_Q_0_1 = fabs((local_map_MAPPRAISER_output[0+1]-CMB_map_nest[index+1])/CMB_map_nest[index+1]);
    max_average_relative_error_U_0_1 = fabs((local_map_MAPPRAISER_output[0+2]-CMB_map_nest[index+2])/CMB_map_nest[index+2]);
    for(index=0;index<size_total_probed;index++){
        index_complete_map = first_pixel + index*nstokes;
        average_relative_error_T_0_1 += fabs((local_map_MAPPRAISER_output[index*nstokes]-CMB_map_nest[index_complete_map])/CMB_map_nest[index_complete_map]);
        if (max_average_relative_error_T_0_1 < fabs((local_map_MAPPRAISER_output[index*nstokes]-CMB_map_nest[index_complete_map])/CMB_map_nest[index_complete_map]));
            max_average_relative_error_T_0_1 = fabs((local_map_MAPPRAISER_output[index*nstokes]-CMB_map_nest[index_complete_map])/CMB_map_nest[index_complete_map]);

        if(CMB_map_nest[index] == 0.000000){
            number_0++;
        }
        if(nstokes>1){
            average_relative_error_Q_0_1 += fabs((local_map_MAPPRAISER_output[index*nstokes+1]-CMB_map_nest[index_complete_map+1])/CMB_map_nest[index_complete_map+1]);
            if (max_average_relative_error_Q_0_1 < fabs((local_map_MAPPRAISER_output[index*nstokes+1]-CMB_map_nest[index_complete_map+1])/CMB_map_nest[index_complete_map+1]));
                max_average_relative_error_Q_0_1 = fabs((local_map_MAPPRAISER_output[index*nstokes+1]-CMB_map_nest[index_complete_map+1])/CMB_map_nest[index_complete_map+1]);


            average_relative_error_U_0_1 += fabs((local_map_MAPPRAISER_output[index*nstokes+2]-CMB_map_nest[index_complete_map+2])/CMB_map_nest[index_complete_map+2]);
            if (max_average_relative_error_U_0_1 < fabs((local_map_MAPPRAISER_output[index*nstokes+2]-CMB_map_nest[index_complete_map+2])/CMB_map_nest[index_complete_map+2])){
                if (fabs(local_map_MAPPRAISER_output[index*nstokes+2]) > 0.01){
                    max_average_relative_error_U_0_1 = fabs((local_map_MAPPRAISER_output[index*nstokes+2]-CMB_map_nest[index_complete_map+2])/CMB_map_nest[index_complete_map+2]);
                    index_U_max = index;
                }
            }

        }
    }
    printf("\n");
    printf("%d --- Average 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", rank, nstokes, npix, 
            average_relative_error_T_0_1/(size_total_probed), average_relative_error_Q_0_1/(size_total_probed), average_relative_error_U_0_1/(size_total_probed));
    printf("%d --- Maximum 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", rank, nstokes, npix, 
            max_average_relative_error_T_0_1, max_average_relative_error_Q_0_1, max_average_relative_error_U_0_1);
    printf("%d --- Number of 0 : %d ; Pixel for max value U : index %d value %f \n", rank, number_0, index_U_max, local_map_MAPPRAISER_output[index_U_max*nstokes+2]); fflush(stdout);
    
    printf("%d --- calculating statistics!!! \n", rank); fflush(stdout);
    max_average_relative_error_0_2 = fabs((local_map_MAPPRAISER_output[0]-local_map_MAPPRAISER[0])/local_map_MAPPRAISER[0]);

    average_relative_error_0_2b = fabs((local_map_MAPPRAISER_output[0]-local_map_MAPPRAISER[0])/local_map_MAPPRAISER[0]);
    printf("%d --- test statistics %f \n", rank, 
            average_relative_error_0_2b);

    // max_average_relative_error_Q_0_2 = fabs((local_map_MAPPRAISER_output[0]-CMB_map_nest[index+1])/CMB_map_nest[index+1]);
    // max_average_relative_error_U_0_2 = fabs((local_map_MAPPRAISER_output[0+2]-CMB_map_nest[index+2])/CMB_map_nest[index+2]);
    for(index=1;index<Obj_pointing_matrix.lcount;index++){
        average_relative_error_0_2 += fabs((local_map_MAPPRAISER_output[index]-local_map_MAPPRAISER[index])/local_map_MAPPRAISER[index]);
        if (max_average_relative_error_0_2 < fabs((local_map_MAPPRAISER_output[index]-local_map_MAPPRAISER[index])/local_map_MAPPRAISER[index]));
            max_average_relative_error_0_2 = fabs((local_map_MAPPRAISER_output[index]-local_map_MAPPRAISER[index])/local_map_MAPPRAISER[index]);
        }
    
    printf("%d --- Average 2bis relative error on all pixels with %d %d : Total %f \n", rank, nstokes, npix, 
            average_relative_error_0_2/(Obj_pointing_matrix.lcount));
    printf("%d --- Maximum 2bis relative error on all pixels with %d %d : Total %f \n", rank, nstokes, npix, 
            max_average_relative_error_0_2); fflush(stdout);


    // free_covariance_matrix(inverse_covariance_matrix, lmax);
    // free_PCG_var(PCG_variable);
    printf("%d --- Free step \n", rank); fflush(stdout);
    free_harmonic_superstruct(&Harmonic_struct, rank);
    free(Obj_pointing_matrix.lindices);
    printf("%d --- Free step 2 \n", rank); fflush(stdout);
    if (S2HAT_params->Local_param_s2hat.gangrank >= 0)
        free(local_alm_s2hat);

    free(local_map_MAPPRAISER);
    free(local_map_MAPPRAISER_output);
    free(CMB_map_nest);
    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}
