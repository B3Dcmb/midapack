
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
#include <chealpix.h>

#include "mappraiser.h"
// #include "domain_generalization.h"

int multiply_square_matrices(int order_matrix, double *matrix_1, double *matrix_2, double *matrix_output);


int multiply_square_matrices(int order_matrix, double *matrix_1, double *matrix_2, double *matrix_output){
    int i, j, k;
    double sum=0;
    for (i=0; i<order_matrix; i++){
        for (k=0; k<order_matrix; k++){
            for(j=0; j<order_matrix; j++){
                sum += matrix_1[i*order_matrix + j] * matrix_2[j*order_matrix + k];
            }
            matrix_output[i*order_matrix + k] = sum;
            sum = 0;
        }
    }
}


int main_first_ver(int argc, char** argv){
// int main(int argc, char** argv){

    char *path_CMB_map = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/Map_test_band_limited.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024.fits";
    int rank, size;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int i, j, ell_value;

    int nside = 512;
    int lmax = 2*nside;

    int nstokes = 3; // T, Q, U
    int number_correlations = 4; // TT, EE, BB, TE


    int npix = 12*nside*nside;

    double *CMB_map_ring = (double *) malloc( nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, nstokes);

    double *CMB_map_nest = (double *) malloc( nstokes*npix*sizeof(double));
    convert_full_map_ring2nest(CMB_map_ring, CMB_map_nest, nside, nstokes);
    printf("%d --- Map transformed !!!  %f %f \n", rank, CMB_map_nest[0], CMB_map_nest[1]); fflush(stdout);
    free(CMB_map_ring);
    


    int overlap = 0; // 10;
    Mat Obj_pointing_matrix;
    Obj_pointing_matrix.trash_pix = 0;
    Obj_pointing_matrix.nnz = nstokes;
    Obj_pointing_matrix.lcount = (npix/size + overlap)*nstokes;
    Obj_pointing_matrix.comm = worldcomm;

    if (rank == size-1)
        Obj_pointing_matrix.lcount += (npix%size - overlap)*nstokes;
    printf("%d --- Producing lcount !!! %d over %d\n", rank, Obj_pointing_matrix.lcount, nstokes*npix); fflush(stdout);

    int number_pixels_MAPP = Obj_pointing_matrix.lcount-(Obj_pointing_matrix.nnz)*(Obj_pointing_matrix.trash_pix);
    
    Obj_pointing_matrix.lindices = (int *)malloc(number_pixels_MAPP*sizeof(int));

    int first_pixel = rank*(npix/size)*nstokes;
    for (i=0; i<Obj_pointing_matrix.lcount; i++){
        Obj_pointing_matrix.lindices[i] = first_pixel + i;
    }
    printf("%d --- Producing lindices !!! first_pixel %d -- 0: %d last: %d \n", rank, first_pixel, Obj_pointing_matrix.lindices[0], Obj_pointing_matrix.lindices[Obj_pointing_matrix.lcount-1]); 
    fflush(stdout);

    double *local_map_MAPPRAISER, *local_map_MAPPRAISER_output;

    local_map_MAPPRAISER = (double *)malloc(Obj_pointing_matrix.lcount*sizeof(double));

    // printf("%d -- Getting min, max of lindices \n", rank);
    // int min_val = Obj_pointing_matrix.lindices[0];
    // int max_val = Obj_pointing_matrix.lindices[0];
    // for (i=1; i<Obj_pointing_matrix.lcount; i++){
    //     if (min_val > Obj_pointing_matrix.lindices[i])
    //         min_val = Obj_pointing_matrix.lindices[i];
    //     if (max_val < Obj_pointing_matrix.lindices[i])
    //         max_val = Obj_pointing_matrix.lindices[i];
    // }
    // printf("%d -- Min indice : %d ; Max indice : %d \n", rank, min_val, max_val); fflush(stdout);
    // printf("%d -- Test last value CMB_map_nest %f \n", rank, CMB_map_nest[nstokes*npix - 1]); fflush(stdout);
    // printf("%d --- Test before setting local_map_MAPPRAISER !!! \n", rank); fflush(stdout);
    for (i=0; i<Obj_pointing_matrix.lcount; i++){
        local_map_MAPPRAISER[i] = CMB_map_nest[Obj_pointing_matrix.lindices[i]];
    }
    free(CMB_map_nest);
    // printf("%d --- Producing local_map_MAPPRAISER !!! first %f last: %f \n", rank, local_map_MAPPRAISER[0], local_map_MAPPRAISER[Obj_pointing_matrix.lcount-1]); fflush(stdout);
    // printf("%d --- From CMB_map_nest !!! first %f last: %f \n", rank, CMB_map_nest[Obj_pointing_matrix.lindices[0]], CMB_map_nest[Obj_pointing_matrix.lindices[Obj_pointing_matrix.lcount-1]]); 
    // fflush(stdout);


    // PCG_var PCG_variable;
    // initialize_PCG_var_struct(&(PCG_variable), local_map_MAPPRAISER);

    int *mask_binary = NULL;
    Harmonic_superstruct Harmonic_struct;
    printf("%d --- initializing harmonic superstruct !!! \n", rank); fflush(stdout);
    init_harmonic_superstruct(&Obj_pointing_matrix, &(Harmonic_struct), mask_binary, nside, lmax, c_ell_path, number_correlations);
    printf("%d --- finishing producing harmonic superstruct !!!!!!! \n", rank); fflush(stdout);
    S2HAT_parameters *S2HAT_params = &(Harmonic_struct.S2HAT_params);

    
    s2hat_dcomplex *local_alm_s2hat;
    if (S2HAT_params->Local_param_s2hat.gangrank >= 0)
        local_alm_s2hat = (s2hat_dcomplex *)malloc(S2HAT_params->nstokes*S2HAT_params->size_alm*sizeof(s2hat_dcomplex));

    printf("%d --- map2harmonic !!! \n", rank); fflush(stdout);
    global_map_2_harmonic(local_map_MAPPRAISER, local_alm_s2hat, &Obj_pointing_matrix, &(Harmonic_struct));
    printf("%d --- map2harmonic end !!! \n", rank); fflush(stdout);

    s2hat_dcomplex *local_alm_s2hat_inverted;
    if (S2HAT_params->Local_param_s2hat.gangrank >= 0)
        local_alm_s2hat_inverted = (s2hat_dcomplex *)malloc(S2HAT_params->nstokes*S2HAT_params->size_alm*sizeof(s2hat_dcomplex));

    double **inverse_covariance_matrix;
    inverse_covariance_matrix = calloc(lmax, sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }

    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    // get_inverse_covariance_matrix_NxN(S2HAT_params, inverse_covariance_matrix);
    // double power_inv_cov = 1/2.;
    double power_inv_cov = 1;
    printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);
    apply_inv_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, S2HAT_params);

    
    printf("%d --- harmonic2map !!! \n", rank); fflush(stdout);
    local_map_MAPPRAISER_output = (double *)malloc(Obj_pointing_matrix.lcount*sizeof(double));
    global_harmonic_2_map(local_map_MAPPRAISER_output, local_alm_s2hat_inverted, &Obj_pointing_matrix, &(Harmonic_struct));
    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);
    
    if (S2HAT_params->Local_param_s2hat.gangrank >= 0){
        free(local_alm_s2hat);
        free(local_alm_s2hat_inverted);
    }

    Obj_pointing_matrix.flag = ALLREDUCE;

    double *full_sky_map_Stokes, *full_map_to_send;
    
    if (rank==0)
        full_sky_map_Stokes = (double *)calloc(nstokes*npix,sizeof(double));
    full_map_to_send = (double *)calloc(nstokes*npix,sizeof(double));
    int *full_sky_map_indices = (int *)calloc(nstokes*npix,sizeof(int));

    for(i=0;i<nstokes*npix;i++)
        full_sky_map_indices[i] = i;

    m2m(local_map_MAPPRAISER_output, Obj_pointing_matrix.lindices, Obj_pointing_matrix.lcount, full_map_to_send, full_sky_map_indices, nstokes*npix);

    // memcpy(full_map_to_send, local_map_MAPPRAISER_output, Obj_pointing_matrix.lcount*sizeof(double));
    free(local_map_MAPPRAISER_output);

    printf("%d --- MPI reduce step ! \n", rank); fflush(stdout);
    MPI_Reduce(full_map_to_send, full_sky_map_Stokes, nstokes*npix, MPI_DOUBLE, MPI_SUM, 0, Obj_pointing_matrix.comm);    
    printf("%d --- End of MPI reduce step ! \n", rank); fflush(stdout);
    free(full_map_to_send);

    if (rank == 0){
        printf("--- Recording maps ! \n"); fflush(stdout);
        char filename_save[80];
        float full_sky_map_1_Stokes[npix];
        for(i=0;i<nstokes; i++)
            {
            // sprintf(filename_save, "/global/cscratch1/sd/mag/WF_work/WF_Tests/mapfile_filtered_bis_%d.fits", i);
            sprintf(filename_save, "/global/cscratch1/sd/mag/WF_work/WF_Tests/mapfile_filtered_wpow_%d.fits", i);
            printf("--- %d -- Saving Stokes components maps \n", i); fflush(stdout);
            for(j=0; j<npix; j++)
                full_sky_map_1_Stokes[j] = full_sky_map_Stokes[j + i*npix];
            write_healpix_map(full_sky_map_1_Stokes, nside, filename_save, 0, "C");
        }
        free(full_sky_map_Stokes);
    }

    printf("%d --- Free step 0 ! \n", rank); fflush(stdout);
    free(local_map_MAPPRAISER);
    free(Obj_pointing_matrix.lindices);
    printf("%d --- Free step 1 ! \n", rank); fflush(stdout);
    
    printf("%d --- Free step 2 ! \n", rank); fflush(stdout);
    printf("%d --- Free step 3 ! \n", rank); fflush(stdout);
    free_covariance_matrix(inverse_covariance_matrix, lmax);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    free_harmonic_superstruct(&Harmonic_struct, rank);
    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}

int main_sec_ver(int argc, char** argv){
// int main(int argc, char** argv){

    char *path_CMB_map = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/Map_test_band_limited.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTE.fits";
    int rank, size;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int i, j, ell_value;

    int nside = 512;
    int lmax = 2*nside;

    int nstokes = 3; // T, Q, U
    int number_correlations = 4; // TT, EE, BB, TE


    int npix = 12*nside*nside;

    double *CMB_map_ring = (double *) malloc( nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, nstokes);

    printf("%d --- Setting T to 0 \n", rank);
    for (i=0; i<npix; i++)
        CMB_map_ring[i] = 0;

    // double *CMB_map_nest = (double *) malloc( nstokes*npix*sizeof(double));
    // convert_full_map_ring2nest(CMB_map_ring, CMB_map_nest, nside, nstokes);
    // printf("%d --- Map transformed !!!  %f %f \n", rank, CMB_map_nest[0], CMB_map_nest[1]); fflush(stdout);
    // free(CMB_map_ring);
    printf("%d --- Initializing S2HAT_param \n", rank); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    int *mask_binary=NULL;
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);
    printf("%d --- End of initializing S2HAT_param \n", rank); fflush(stdout);

    double *local_map_pix = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("%d --- Distributing S2HAT_param \n", rank); fflush(stdout);
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring, local_map_pix, &S2HAT_params);
    printf("%d --- End of distributing S2HAT_param \n", rank); fflush(stdout);
    // int overlap = 0; // 10;
    // Mat Obj_pointing_matrix;
    // Obj_pointing_matrix.trash_pix = 0;
    // Obj_pointing_matrix.nnz = nstokes;
    // Obj_pointing_matrix.lcount = (npix/size + overlap)*nstokes;
    // Obj_pointing_matrix.comm = worldcomm;

    // if (rank == size-1)
    //     Obj_pointing_matrix.lcount += (npix%size - overlap)*nstokes;
    // printf("%d --- Producing lcount !!! %d over %d\n", rank, Obj_pointing_matrix.lcount, nstokes*npix); fflush(stdout);

    // int number_pixels_MAPP = Obj_pointing_matrix.lcount-(Obj_pointing_matrix.nnz)*(Obj_pointing_matrix.trash_pix);
    
    // Obj_pointing_matrix.lindices = (int *)malloc(number_pixels_MAPP*sizeof(int));

    // int first_pixel = rank*(npix/size)*nstokes;
    // for (i=0; i<Obj_pointing_matrix.lcount; i++){
    //     Obj_pointing_matrix.lindices[i] = first_pixel + i;
    // }
    // printf("%d --- Producing lindices !!! first_pixel %d -- 0: %d last: %d \n", rank, first_pixel, Obj_pointing_matrix.lindices[0], Obj_pointing_matrix.lindices[Obj_pointing_matrix.lcount-1]); 
    // fflush(stdout);

    // double *local_map_MAPPRAISER, *local_map_MAPPRAISER_output;

    // local_map_MAPPRAISER = (double *)malloc(Obj_pointing_matrix.lcount*sizeof(double));

    // printf("%d -- Getting min, max of lindices \n", rank);
    // int min_val = Obj_pointing_matrix.lindices[0];
    // int max_val = Obj_pointing_matrix.lindices[0];
    // for (i=1; i<Obj_pointing_matrix.lcount; i++){
    //     if (min_val > Obj_pointing_matrix.lindices[i])
    //         min_val = Obj_pointing_matrix.lindices[i];
    //     if (max_val < Obj_pointing_matrix.lindices[i])
    //         max_val = Obj_pointing_matrix.lindices[i];
    // }
    // printf("%d -- Min indice : %d ; Max indice : %d \n", rank, min_val, max_val); fflush(stdout);
    // printf("%d -- Test last value CMB_map_nest %f \n", rank, CMB_map_nest[nstokes*npix - 1]); fflush(stdout);
    // printf("%d --- Test before setting local_map_MAPPRAISER !!! \n", rank); fflush(stdout);
    // for (i=0; i<Obj_pointing_matrix.lcount; i++){
    //     local_map_MAPPRAISER[i] = CMB_map_nest[Obj_pointing_matrix.lindices[i]];
    // }
    // free(CMB_map_nest);
    // printf("%d --- Producing local_map_MAPPRAISER !!! first %f last: %f \n", rank, local_map_MAPPRAISER[0], local_map_MAPPRAISER[Obj_pointing_matrix.lcount-1]); fflush(stdout);
    // printf("%d --- From CMB_map_nest !!! first %f last: %f \n", rank, CMB_map_nest[Obj_pointing_matrix.lindices[0]], CMB_map_nest[Obj_pointing_matrix.lindices[Obj_pointing_matrix.lcount-1]]); 
    // fflush(stdout);


    // PCG_var PCG_variable;
    // initialize_PCG_var_struct(&(PCG_variable), local_map_MAPPRAISER);

    // int *mask_binary = NULL;
    // Harmonic_superstruct Harmonic_struct;
    // printf("%d --- initializing harmonic superstruct !!! \n", rank); fflush(stdout);
    // init_harmonic_superstruct(&Obj_pointing_matrix, &(Harmonic_struct), mask_binary, nside, lmax, c_ell_path, number_correlations);
    // printf("%d --- finishing producing harmonic superstruct !!!!!!! \n", rank); fflush(stdout);
    // S2HAT_parameters *S2HAT_params = &(Harmonic_struct.S2HAT_params);

    
    s2hat_dcomplex *local_alm_s2hat;
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0)
        local_alm_s2hat = (s2hat_dcomplex *)malloc(S2HAT_params.nstokes*S2HAT_params.size_alm*sizeof(s2hat_dcomplex));

    

    printf("%d --- map2harmonic !!! \n", rank); fflush(stdout);
    // global_map_2_harmonic(local_map_MAPPRAISER, local_alm_s2hat, &Obj_pointing_matrix, &(Harmonic_struct));
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);
    printf("%d --- map2harmonic end !!! \n", rank); fflush(stdout);


    s2hat_dcomplex *local_alm_s2hat_inverted;
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0)
        local_alm_s2hat_inverted = (s2hat_dcomplex *)malloc(S2HAT_params.nstokes*S2HAT_params.size_alm*sizeof(s2hat_dcomplex));

    double **inverse_covariance_matrix;
    inverse_covariance_matrix = calloc(lmax, sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    // get_inverse_covariance_matrix_NxN(&S2HAT_params, inverse_covariance_matrix);
    
    printf("%d --- printing inverse of covariance matrix !!! \n", rank); fflush(stdout);
    int order_matrix = nstokes;
    int size_ell_max = 10;
    int index, index_2;
    for (ell_value=2;ell_value<size_ell_max;ell_value++){
        printf("\n #### ell= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
            printf("\n");
        }
    }
    fflush(stdout);
    // double power_inv_cov = 1/2.;
    // double power_inv_cov = 1/2.;
    // double power_inv_cov = 0;
    printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);
    apply_inv_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, &S2HAT_params);

    
    printf("%d --- harmonic2map !!! \n", rank); fflush(stdout);
    // local_map_MAPPRAISER_output = (double *)malloc(Obj_pointing_matrix.lcount*sizeof(double));
    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    // global_harmonic_2_map(local_map_MAPPRAISER_output, local_alm_s2hat_inverted, &Obj_pointing_matrix, &(Harmonic_struct));
    apply_alm2pix(local_map_output, local_alm_s2hat_inverted, &S2HAT_params);
    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);
    
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0){
        free(local_alm_s2hat);
        free(local_alm_s2hat_inverted);
    }

    // Obj_pointing_matrix.flag = ALLREDUCE;

    // double *full_sky_map_Stokes, *full_map_to_send;
    
    // if (rank==0)
    //     full_sky_map_Stokes = (double *)calloc(nstokes*npix,sizeof(double));
    // full_map_to_send = (double *)calloc(nstokes*npix,sizeof(double));
    // int *full_sky_map_indices = (int *)calloc(nstokes*npix,sizeof(int));

    // for(i=0;i<nstokes*npix;i++)
    //     full_sky_map_indices[i] = i;

    // m2m(local_map_MAPPRAISER_output, Obj_pointing_matrix.lindices, Obj_pointing_matrix.lcount, full_map_to_send, full_sky_map_indices, nstokes*npix);

    // // memcpy(full_map_to_send, local_map_MAPPRAISER_output, Obj_pointing_matrix.lcount*sizeof(double));
    // free(local_map_MAPPRAISER_output);
    double *full_sky_map_2;
    full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
    gather_map(local_map_output, full_sky_map_2, nstokes, &S2HAT_params);

    free(local_map_pix);
    free(local_map_output);
    // printf("%d --- MPI reduce step ! \n", rank); fflush(stdout);
    // MPI_Reduce(full_map_to_send, full_sky_map_Stokes, nstokes*npix, MPI_DOUBLE, MPI_SUM, 0, Obj_pointing_matrix.comm);    
    // printf("%d --- End of MPI reduce step ! \n", rank); fflush(stdout);
    // free(full_map_to_send);

    if (rank == 0){
        printf("--- Recording maps ! \n"); fflush(stdout);
        char filename_save[80];
        float *full_sky_map_1_Stokes = (float *)malloc(npix*sizeof(float));
        for(i=0;i<nstokes; i++)
            {
            // sprintf(filename_save, "/global/cscratch1/sd/mag/WF_work/WF_Tests/mapfile_filtered_bis_%d.fits", i);
            sprintf(filename_save, "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/map_files/mapfile_filtered_wonest_woTT_woTE_neworder_%d.fits", i);
            printf("--- %d -- Saving Stokes components maps \n", i); fflush(stdout);
            for(j=0; j<npix; j++)
                full_sky_map_1_Stokes[j] = full_sky_map_2[j + i*npix];
            write_healpix_map(full_sky_map_1_Stokes, nside, filename_save, 0, "C");
        }
        
        free(full_sky_map_1_Stokes);
    }
    

    printf("%d --- Free step 0 ! \n", rank); fflush(stdout);
    // free(local_map_MAPPRAISER);
    // free(Obj_pointing_matrix.lindices);
    printf("%d --- Free step 1 ! \n", rank); fflush(stdout);
    
    printf("%d --- Free step 2 ! \n", rank); fflush(stdout);
    printf("%d --- Free step 3 ! \n", rank); fflush(stdout);
    free_covariance_matrix(inverse_covariance_matrix, lmax);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    free(full_sky_map_2);
    // free_harmonic_superstruct(&Harmonic_struct, rank);
    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}

// int main_test_read_cells(int argc, char** argv){
int main(int argc, char** argv){

    char *path_CMB_map = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/Map_test_band_limited.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTE.fits";
    int rank, size;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int i, j, ell_value;

    int nside = 512;
    int lmax = 2*nside;
    // int lmax = 20;

    int nstokes = 3; // T, Q, U
    int number_correlations = 4; // TT, EE, BB, TE


    int npix = 12*nside*nside;

    double *CMB_map_ring = (double *) malloc( nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, nstokes);

    // printf("%d --- Setting T to 0 \n", rank);
    // for (i=0; i<npix; i++)
    //     CMB_map_ring[i] = 0;

    printf("%d --- Initializing S2HAT_param \n", rank); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    int *mask_binary=NULL;
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);
    printf("%d --- End of initializing S2HAT_param \n", rank); fflush(stdout);

    printf("%d --- getting cells of covariance matrix !!! \n", rank); fflush(stdout);
    double *c_ell_array = (double *)calloc(number_correlations*(lmax),sizeof(double));
    int not_block_diagonal = 0; // We want block diagonal c_ells
    read_fits_cells(lmax, number_correlations, c_ell_array, c_ell_path, 1, not_block_diagonal); // Reading cell_fits_file

    int correl_index = 0;
    int max_size_ells_to_probe = 5;
    printf("%d --- cell_TT = %.9f -", rank, c_ell_array[correl_index*(lmax)]);
    for(ell_value=1; ell_value<max_size_ells_to_probe; ell_value++){
        printf("- %.9f -", rank, c_ell_array[correl_index*(lmax)+ell_value]);
    }
    printf("\n"); fflush(stdout);
    
    correl_index = 1;
    printf("%d --- cell_EE = %.9f -", rank, c_ell_array[correl_index*(lmax)]);
    for(ell_value=1; ell_value<max_size_ells_to_probe; ell_value++){
        printf("- %.9f -", rank, c_ell_array[correl_index*(lmax)+ell_value]);
    }
    printf("\n"); fflush(stdout);

    correl_index = 2;
    printf("%d --- cell_BB = %.9f -", rank, c_ell_array[correl_index*(lmax)]);
    for(ell_value=1; ell_value<max_size_ells_to_probe; ell_value++){
        printf("- %.9f -", rank, c_ell_array[correl_index*(lmax)+ell_value]);
    }
    printf("\n"); fflush(stdout);

    correl_index = 3;
    printf("%d --- cell_TE = %.9f -", rank, c_ell_array[correl_index*(lmax)]);
    for(ell_value=1; ell_value<max_size_ells_to_probe; ell_value++){
        printf("- %.9f -", rank, c_ell_array[correl_index*(lmax)+ell_value]);
    }
    printf("\n"); fflush(stdout);

    double **covariance_matrix_NxN = (double **)malloc(lmax*sizeof(double *));
    for (ell_value=0; ell_value<lmax; ell_value++){
        covariance_matrix_NxN[ell_value] = (double *)calloc(nstokes*nstokes, sizeof(double));
    }

    // double **inv_covariance_matrix_NxN = (double **)malloc(lmax*sizeof(double *));
    // for (ell_value=0; ell_value<lmax; ell_value++){
    //     inv_covariance_matrix_NxN[ell_value] = (double *)calloc(nstokes*nstokes, sizeof(double));
    // }
    double *inv_covariance_matrix_NxN = (double *)calloc(nstokes*nstokes,sizeof(double));
    // for (ell_value=0; ell_value<lmax; ell_value++){
    //     inv_covariance_matrix_NxN[ell_value] = (double *)calloc(nstokes*nstokes, sizeof(double));
    // }
    

    printf("--- Getting covariance matrix \n"); fflush(stdout);
    get_covariance_matrix_block_diagonal(c_ell_path, number_correlations, covariance_matrix_NxN, &S2HAT_params);
    printf("--- Covariance matrix obtained ! \n"); fflush(stdout);
    
    int index, index_2, ell_index;
    int order_matrix = nstokes;
    char cholesky_part = 'L';
    int info;
    double cholesky_factor[nstokes*nstokes];

    
    for (ell_value=1;ell_value<max_size_ells_to_probe;ell_value++){
        printf("\n ####0000 ell= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.8f \t", covariance_matrix_NxN[ell_value][index*order_matrix + index_2]);
            }
            printf("\n"); fflush(stdout);
        }
        memcpy(cholesky_factor, covariance_matrix_NxN[ell_value], nstokes*nstokes*sizeof(double));
        // printf("\n ####Cholesky ell= %d \n", ell_value);
        // for(ell_index=0; ell_index<nstokes*(nstokes+1)/2; ell_index++){
        //     cholesky_factor[ell_index] = covariance_matrix_NxN[ell_value][ell_index + ell_index/nstokes];
        //     printf("%.8f %d+%d \t", cholesky_factor[ell_index], ell_index, ell_index/nstokes);
        // }
        // if (nstokes == 3){
        //         cholesky_factor[nstokes*(nstokes+1)/2 - 1] = covariance_matrix_NxN[ell_value][nstokes*(nstokes+1)/2 + 2];
        //         printf("corr %.8f %d+%d \t", cholesky_factor[ell_index], nstokes*(nstokes+1)/2 - 1, nstokes*(nstokes+1)/2 + 2);
        // }
        // printf("\n");
        double cholesky_factor_2[nstokes*nstokes];
        memcpy(cholesky_factor_2, covariance_matrix_NxN[ell_value], nstokes*nstokes*sizeof(double));
        dpotrf_(&cholesky_part, &order_matrix, cholesky_factor_2, &order_matrix, &info);
        

        printf("\n ####Cholesky_true info %d ell= %d \n", info, ell_value); fflush(stdout);
        // for(ell_index=0; ell_index<nstokes*(nstokes+1)/2; ell_index++){
        for(ell_index=0; ell_index<nstokes*nstokes; ell_index++){
            printf("%.8f , \t", cholesky_factor_2[ell_index]); fflush(stdout);
        }
        printf("\n"); fflush(stdout);

        // memcpy(cholesky_factor_2, covariance_matrix_NxN[ell_value], nstokes*nstokes*sizeof(double));
        dpotri_(&cholesky_part, &order_matrix, cholesky_factor_2, &order_matrix, &info);

        printf("\n ####Cholesky_true_inv info %d ell= %d \n", info, ell_value); fflush(stdout);
        // for(ell_index=0; ell_index<nstokes*(nstokes+1)/2; ell_index++){
        for(ell_index=0; ell_index<nstokes*nstokes; ell_index++){
            printf("%.8f , \t", cholesky_factor_2[ell_index]); fflush(stdout);
        }
        printf("\n"); fflush(stdout);

        printf("--- Getting Cholesky decomposition and inversion \n", ell_value); fflush(stdout);
        get_cholesky_decomposition_inverted(order_matrix, cholesky_factor, cholesky_part);
        printf("--- Cholesky decomposition and inversion obtained ! \n", ell_value); fflush(stdout);

        printf("\n ####Cholesky_inv ell= %d \n", ell_value); fflush(stdout);
        // for(ell_index=0; ell_index<nstokes*(nstokes+1)/2; ell_index++){
        for(ell_index=0; ell_index<nstokes*nstokes; ell_index++){
            printf("%.8f , \t", cholesky_factor[ell_index]); fflush(stdout);
        }
        printf("\n"); fflush(stdout);

        double test_matrix[nstokes*nstokes];
        int corr=0;
        for (index=0; index < order_matrix; index++){
            // if (index*nstokes == nstokes*(nstokes+1)/2)
            //     corr = -1;

            // test_matrix[index*order_matrix + index] = cholesky_factor[index*nstokes - corr];
            // for (index_2=0; index_2 < order_matrix - index; index_2++){
                // test_matrix[index*order_matrix + index_2] = cholesky_factor[index*nstokes + index_2];
                // test_matrix[index_2*order_matrix + index] = cholesky_factor[index*nstokes + index_2];
            for (index_2=0; index_2 < order_matrix; index_2++){
                test_matrix[index*order_matrix + index_2] = cholesky_factor[index*nstokes + index_2];
            }
        }
        test_matrix[3] = test_matrix[1];
        test_matrix[6] = test_matrix[2];
        test_matrix[7] = test_matrix[5];

        printf("\n ####test_reconstruct_Cholesky_0 ell= %d \n", ell_value); fflush(stdout);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.8f \t", test_matrix[index*order_matrix + index_2]);
            }
            printf("\n"); fflush(stdout);
        }
        double output_test_matrix[nstokes*nstokes];
        multiply_square_matrices(order_matrix, test_matrix, covariance_matrix_NxN[ell_value], output_test_matrix);

        printf("\n ####test_reconstruct_unit_Cholesky ell= %d \n", ell_value); fflush(stdout);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.8f \t", output_test_matrix[index*order_matrix + index_2]);
            }
            printf("\n"); fflush(stdout);
        }

    
        // dpptri_(&cholesky_part, &order_matrix, cholesky_factor, &info);

        // double test_matrix_2[nstokes*nstokes];
        // corr=0;
        // for (index=0; index < order_matrix; index++){
        //     if (index*nstokes == nstokes*(nstokes+1)/2)
        //         corr = -1;

        //     test_matrix_2[index*order_matrix + index] = cholesky_factor[index*nstokes - corr];
        //     for (index_2=0; index_2 < order_matrix - index; index_2++){
        //         test_matrix_2[index*order_matrix + index_2] = cholesky_factor[index*nstokes + index_2];
        //     }
        // }
        // double output_test_matrix_2[nstokes*nstokes];
        // multiply_square_matrices(order_matrix, output_test_matrix, covariance_matrix_NxN[ell_value], output_test_matrix_2);

        // printf("\n ####test_reconstruct_unit_from_Cholesky ell= %d \n", ell_value); fflush(stdout);
        // for (index=0; index < order_matrix; index++){
        //     for (index_2=0; index_2<order_matrix; index_2++){
        //         printf("%.8f \t", output_test_matrix_2[index*order_matrix + index_2]);
        //     }
        //     printf("\n"); fflush(stdout);
        // }

        printf("\n ####Cholesky_inv ell= %d \n", ell_value); fflush(stdout);
        for(ell_index=0; ell_index<nstokes*(nstokes+1)/2; ell_index++){
            inv_covariance_matrix_NxN[ell_index + (ell_index/nstokes)*nstokes] = cholesky_factor[ell_index];
            printf("%.8f \t", cholesky_factor[ell_index]);
        }
        printf("\n"); fflush(stdout);
        inv_covariance_matrix_NxN[nstokes] = cholesky_factor[nstokes-2];
        inv_covariance_matrix_NxN[2*nstokes] = cholesky_factor[nstokes-1];
        inv_covariance_matrix_NxN[2*nstokes+1] = cholesky_factor[2*nstokes-1];


        printf("\n ####inv ell= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.8f \t", inv_covariance_matrix_NxN[index*order_matrix + index_2]);
            }
            printf("\n");
        }
    }
    fflush(stdout);
    free(inv_covariance_matrix_NxN);

    double **inverse_covariance_matrix, **cholesky_decomposition;
    inverse_covariance_matrix = calloc(lmax, sizeof(double *));
    cholesky_decomposition = calloc(lmax, sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
        cholesky_decomposition[ell_value] = calloc(nstokes*(nstokes+1)/2,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    get_inverse_covariance_matrix_NxN(&S2HAT_params, inverse_covariance_matrix);
    printf("%d --- printing inverse of covariance matrix !!! \n", rank); fflush(stdout);

    int size_ell_max = 10;
    
    for (ell_value=2;ell_value<size_ell_max;ell_value++){
        printf("\n #### ell= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
            printf("\n");
        }
    }
    fflush(stdout);
    // double power_inv_cov = 1/2.;
    // double power_inv_cov = 1/2.;
    double power_inv_cov = 0;
    // printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);
    // apply_inv_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, power_inv_cov, &S2HAT_params);

    
    printf("%d --- Free step 0 ! \n", rank); fflush(stdout);
    // free(local_map_MAPPRAISER);
    // free(Obj_pointing_matrix.lindices);
    free(c_ell_array);
    printf("%d --- Free step 1 ! \n", rank); fflush(stdout);
    
    printf("%d --- Free step 2 ! \n", rank); fflush(stdout);
    printf("%d --- Free step 3 ! \n", rank); fflush(stdout);
    free_covariance_matrix(inverse_covariance_matrix, lmax-1);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    // free(full_sky_map_2);
    // free_harmonic_superstruct(&Harmonic_struct, rank);
    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}

// With all Stokes params
int main_clean_ver_0(int argc, char** argv){
// int main(int argc, char** argv){

    // char *path_CMB_map = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/Map_test_band_limited.fits";
    char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_band_limited_1024_0.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTE.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTTTE.fits";
    int rank, size;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int i, j, ell_value;

    int nside = 512;
    int lmax = 2*nside;

    int nstokes = 3; // T, Q, U
    int number_correlations = 4; // TT, EE, BB, TE


    int npix = 12*nside*nside;
    int index, index_2;
    int order_matrix = nstokes;

    double *CMB_map_ring = (double *) malloc( nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, nstokes);

    printf("%d --- Setting T to 0 !!!! \n", rank);
    for (i=0; i<npix; i++)
        CMB_map_ring[i] = 0;

    int max_size_pixels = 10;
    printf("%d --- Reading first elems of map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f -", CMB_map_ring[i]);
    printf("\n");

    printf("%d --- Initializing S2HAT_param \n", rank); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    int *mask_binary=NULL;
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);
    printf("%d --- End of initializing S2HAT_param \n", rank); fflush(stdout);

    double *local_map_pix = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("%d --- Distributing S2HAT_param \n", rank); fflush(stdout);
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring, local_map_pix, &S2HAT_params);
    printf("%d --- End of distributing S2HAT_param \n", rank); fflush(stdout);
    
    printf("%d --- Reading first elems of local map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", local_map_pix[i], CMB_map_ring[i]);
    printf("\n");

    s2hat_dcomplex *local_alm_s2hat;
    // if (S2HAT_params.Local_param_s2hat.gangrank >= 0)
    local_alm_s2hat = (s2hat_dcomplex *)malloc(nstokes*S2HAT_params.size_alm*sizeof(s2hat_dcomplex));

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);
    // int nrings = Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1;
    // double local_w8ring[nrings*nstokes];
    // int i_ring;
    // for( i_ring=0; i_ring< nstokes*nrings; i_ring++)
    //         local_w8ring[i_ring]=1.;
    int spin=0;            
    // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, 1, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, S2HAT_params.lda, local_alm_s2hat,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    // s2hat_map2alm(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, 1, 1, 
    //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, local_w8ring, Local_param_s2hat->map_size, local_map_pix, S2HAT_params.lda, local_alm_s2hat, 
    //     Local_param_s2hat->nplm, NULL,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    
    int max_size_test = 50;
    
    // printf("%d --- Local_alm after 1st pix2alm - %f %f -", rank, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm_s2hat[index].re, local_alm_s2hat[index].im);
    // }
    // printf(" \n");
    spin=2;            
    // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, 1, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, S2HAT_params.lda, local_alm_s2hat+S2HAT_params.size_alm,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    // printf("%d --- map2harmonic end !!! \n", rank); fflush(stdout);

   


    double **covariance_matrix;
    covariance_matrix = calloc(lmax, sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    get_covariance_matrix_block_diagonal(c_ell_path, number_correlations, covariance_matrix, &S2HAT_params);
    double **inverse_covariance_matrix;
    inverse_covariance_matrix = calloc(lmax, sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    get_inverse_covariance_matrix_NxN(&S2HAT_params, inverse_covariance_matrix);

    // Getting unit
    // for (ell_value=0;ell_value<lmax;ell_value++){
    //     for (index=0; index < order_matrix; index++){
    //         for (index_2=0; index_2 < order_matrix; index_2++){
    //             if (index == index_2)
    //                 inverse_covariance_matrix[ell_value][index*order_matrix + index_2] = 1;
    //             else{
    //                 inverse_covariance_matrix[ell_value][index*order_matrix + index_2] = 0;
    //             }
    //         }
    //     }
    // }

    printf("%d --- printing inverse of covariance matrix !!! \n", rank); fflush(stdout);
    
    int size_ell_max = 10;
    
    for (ell_value=0;ell_value<size_ell_max;ell_value++){
        printf("\n #### ell= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
            printf("\n");
        }
    }
    fflush(stdout);

    char cholesky_part = 'L';
    int info;
    double *transpose_test = (double *)calloc(nstokes*nstokes,sizeof(double));
    double transpose_test_2[nstokes*nstokes];
    double transpose_test_3[nstokes*nstokes];
    for(ell_value=0; ell_value<lmax; ell_value++){
        dpotrf_(&cholesky_part, &nstokes, inverse_covariance_matrix[ell_value], &nstokes, &info);
        inverse_covariance_matrix[ell_value][nstokes] = 0;
        if (nstokes == 3){
            inverse_covariance_matrix[ell_value][6] = 0;
            inverse_covariance_matrix[ell_value][7] = 0;
        }
        // inverse_covariance_matrix[ell_value][nstokes] = inverse_covariance_matrix[ell_value][1];
        // if (nstokes == 3){
        //     inverse_covariance_matrix[ell_value][6] = inverse_covariance_matrix[ell_value][2];
        //     inverse_covariance_matrix[ell_value][7] = inverse_covariance_matrix[ell_value][5];
        // }

        // inverse_covariance_matrix[ell_value][1] = 0;
        // if (nstokes == 3){
        //     inverse_covariance_matrix[ell_value][2] = 0;
        //     inverse_covariance_matrix[ell_value][5] = 0;
        // }
        // if (ell_value<15){
        //     for(i=0; i<nstokes; i++)
        //         transpose_test[i*nstokes + i] = inverse_covariance_matrix[ell_value][i*nstokes + i];
            
        //     transpose_test[nstokes] = inverse_covariance_matrix[ell_value][nstokes];
        //     if (nstokes == 3){
        //         transpose_test[6] = inverse_covariance_matrix[ell_value][6];
        //         transpose_test[7] = inverse_covariance_matrix[ell_value][7];
        //     }
        //     transpose_test[1] = 0;
        //     if (nstokes == 3){
        //         transpose_test[2] = 0;
        //         transpose_test[5] = 0;
        //     }
        //     // transpose_test[1] = inverse_covariance_matrix[ell_value][nstokes];
        //     // if (nstokes == 3){
        //     //     transpose_test[2] = inverse_covariance_matrix[ell_value][6];
        //     //     transpose_test[5] = inverse_covariance_matrix[ell_value][7];
        //     // }

        //     multiply_square_matrices(nstokes, covariance_matrix[ell_value], transpose_test, transpose_test_2);
        //     // transpose_test[nstokes] = transpose_test[1];
        //     // if (nstokes == 3){
        //     //     transpose_test[6] = transpose_test[2];
        //     //     transpose_test[7] = transpose_test[5];
        //     // }
        //     // transpose_test[1] = 0;
        //     // if (nstokes == 3){
        //     //     transpose_test[2] = 0;
        //     //     transpose_test[5] = 0;
        //     // }
        //     transpose_test[1] = inverse_covariance_matrix[ell_value][nstokes];
        //     if (nstokes == 3){
        //         transpose_test[2] = inverse_covariance_matrix[ell_value][6];
        //         transpose_test[5] = inverse_covariance_matrix[ell_value][7];
        //     }
        //     transpose_test[nstokes] = 0;
        //     if (nstokes == 3){
        //         transpose_test[6] = 0;
        //         transpose_test[7] = 0;
        //     }
        //     multiply_square_matrices(nstokes, transpose_test, transpose_test_2, transpose_test_3);
        //     printf("\n ####inv_transp ell= %d \n", ell_value);
        //     for (index=0; index < order_matrix; index++){
        //         for (index_2=0; index_2<order_matrix; index_2++){
        //             printf("%.7f \t", transpose_test_3[index*order_matrix + index_2]);
        //         }
        //     printf("\n");
        // }
        // }
    }
    free(transpose_test);
    

    max_size_test = 50;
    
    // printf("%d --- Local_alm long after 1st pix2alm - %f %f -", rank, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm_s2hat[index].re, local_alm_s2hat[index].im);
    // }
    // printf(" \n");

    int min_alms = lmax, max_alms = min_alms + max_size_test;
    // printf("%d --- Local_alm long after 1st pix2alm --- alms %d to %d - %f %f -", rank, min_alms, max_alms, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    // for (index=min_alms;index<max_size_test+max_alms;index++){
    //     printf("- %f %f -", local_alm_s2hat[index].re, local_alm_s2hat[index].im);
    // }
    // printf(" \n");

    s2hat_dcomplex *local_alm_s2hat_inverted;
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0)
        local_alm_s2hat_inverted = (s2hat_dcomplex *)malloc(S2HAT_params.nstokes*S2HAT_params.size_alm*sizeof(s2hat_dcomplex));

    printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);

    apply_inv_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, &S2HAT_params);
    // int m_value, index_stokes, line_index, nmvals;
    // double res_real, res_imag;
    // for(ell_value=0; ell_value < lmax; ell_value++){
    //     for(m_value=0; m_value < S2HAT_params.Local_param_s2hat.nmvals; m_value++){
    //         for (index_stokes=0; index_stokes<nstokes; index_stokes++){
    //             res_real = 0;
    //             res_imag = 0;
    //             for (line_index=0; line_index < nstokes; line_index++){
    //                 res_real += inverse_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * local_alm_s2hat[line_index*nmvals*(lmax) + m_value*(lmax) + ell_value].re;
    //                 res_imag += inverse_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * local_alm_s2hat[line_index*nmvals*(lmax) + m_value*(lmax) + ell_value].im;
    //             }
    //             local_alm_s2hat_inverted[index_stokes*nmvals*(lmax) + m_value*(lmax) + ell_value].re = res_real;
    //             local_alm_s2hat_inverted[index_stokes*nmvals*(lmax) + m_value*(lmax) + ell_value].im = res_imag;
    //         }
    //     }
    // }

    // printf("%d --- Local_alm after apply inv cov matrix - %f %f -", rank, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm_s2hat[index].re, local_alm_s2hat[index].im);
    // }
    // printf(" \n");
    // printf("%d --- Local_alm long after apply inv cov matrix --- alms %d to %d - %f %f -", rank, min_alms, max_alms, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    // for (index=min_alms;index<max_size_test+max_alms;index++){
    //     printf("- %f %f -", local_alm_s2hat[index].re, local_alm_s2hat[index].im);
    // }
    // printf(" \n");
    // printf("%d --- Local_alm after apply inv cov matrix - %f %f -", rank, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm_s2hat_inverted[index].re, local_alm_s2hat_inverted[index].im);
    // }
    // printf(" \n");
    // printf("%d --- Local_alm long after apply inv cov matrix --- alms %d to %d - %f %f -", rank, min_alms, max_alms, local_alm_s2hat[0].re, local_alm_s2hat[0].im);
    // for (index=min_alms;index<max_size_test+max_alms;index++){
    //     printf("- %f %f -", local_alm_s2hat_inverted[index].re, local_alm_s2hat_inverted[index].im);
    // }
    // printf(" \n");
    // // printf("%d --- Local_alm after apply_inv_cov_matrix - %f %f -", rank, local_alm_s2hat_inverted[0].re, local_alm_s2hat_inverted[0].im);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_alm_s2hat_inverted[index].re, local_alm_s2hat_inverted[index].im);
    // // }
    // // printf(" \n");
    
    printf("%d --- harmonic2map !!! \n", rank); fflush(stdout);
    // local_map_MAPPRAISER_output = (double *)malloc(Obj_pointing_matrix.lcount*sizeof(double));
    double *local_map_output = (double *)calloc(nstokes*S2HAT_params.Local_param_s2hat.map_size,sizeof(double));
    // global_harmonic_2_map(local_map_MAPPRAISER_output, local_alm_s2hat_inverted, &Obj_pointing_matrix, &(Harmonic_struct));
    // apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);
    apply_alm2pix(local_map_output, local_alm_s2hat_inverted, &S2HAT_params);

    // s2hat_alm2map(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, 1, 1, 
    //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_output, S2HAT_params.lda, 
    //     local_alm_s2hat_inverted, Local_param_s2hat->nplm, NULL, 
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    
    // printf("%d --- Reading0 first elems of local map output \n", rank);
    // for (i=0; i<max_size_pixels; i++)
    //     printf("- %f %f -", local_map_output[i], CMB_map_ring[i]);
    // printf("\n");

    // spin=2;
    // s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, 1,
    //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_output+Local_param_s2hat->map_size, S2HAT_params.lda, local_alm_s2hat_inverted+S2HAT_params.size_alm,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);

    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);

    printf("%d --- Reading first elems of local map output \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", local_map_output[i], CMB_map_ring[i]);
    printf("\n");

    double *full_sky_map_2;
    full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
    gather_map(local_map_output, full_sky_map_2, nstokes, &S2HAT_params);

    printf("%d --- Reading first elems of final maps retrieved \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", full_sky_map_2[i], CMB_map_ring[i]);
    printf("\n");

    // printf("%d --- MPI reduce step ! \n", rank); fflush(stdout);
    // MPI_Reduce(full_map_to_send, full_sky_map_Stokes, nstokes*npix, MPI_DOUBLE, MPI_SUM, 0, Obj_pointing_matrix.comm);    
    // printf("%d --- End of MPI reduce step ! \n", rank); fflush(stdout);
    // free(full_map_to_send);

    if (rank == 0){
        printf("--- Recording maps ! \n"); fflush(stdout);
        char filename_save[80];
        float *full_sky_map_1_Stokes = (float *)malloc(npix*sizeof(float));
        for(i=0;i<nstokes; i++)
            {
            // sprintf(filename_save, "/global/cscratch1/sd/mag/WF_work/WF_Tests/mapfile_filtered_bis_%d.fits", i);
            // sprintf(filename_save, "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/map_files/mapfile_filtered_newver0_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_unit6_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver6_%d.fits", i);
            sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_woTTTE_0_%d.fits", i);
            printf("--- %d -- Saving Stokes components maps \n", i); fflush(stdout);
            for(j=0; j<npix; j++)
                full_sky_map_1_Stokes[j] = full_sky_map_2[j + i*npix];
            write_healpix_map(full_sky_map_1_Stokes, nside, filename_save, 0, "C");
        }
        free(full_sky_map_1_Stokes);
    }


    printf("%d --- Free step 0 ! \n", rank); fflush(stdout);
    free(CMB_map_ring);
    printf("%d --- Free step 1 ! \n", rank); fflush(stdout);
    free(full_sky_map_2);
    printf("%d --- Free step 2 ! \n", rank); fflush(stdout);
    // free(local_map_pix);
    printf("%d --- Free step 2b ! \n", rank); fflush(stdout);
    free(local_map_output);
    printf("%d --- Free step 3 ! \n", rank); fflush(stdout);
    free_covariance_matrix(inverse_covariance_matrix, lmax-1);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    // free_harmonic_superstruct(&Harmonic_struct, rank);
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0){
        printf("%d --- Free step 1b ! \n", rank); fflush(stdout);
        free(local_alm_s2hat_inverted);
        printf("%d --- Free step 1a ! \n", rank); fflush(stdout);
        free(local_alm_s2hat);
    }
    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}
