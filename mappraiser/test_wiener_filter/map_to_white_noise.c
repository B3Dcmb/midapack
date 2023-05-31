
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
#include <chealpix.h>

#include "mappraiser.h"
// #include "domain_generalization.h"

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
    inverse_covariance_matrix = calloc(lmax+1, sizeof(double *));
    for(ell_value=0; ell_value<lmax+1; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }

    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    get_inverse_covariance_matrix_NxN(S2HAT_params, inverse_covariance_matrix);
    // double power_inv_cov = 1/2.;
    double power_inv_cov = 1;
    printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);
    apply_inv_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, power_inv_cov, S2HAT_params);

    
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

// int main_sec_ver(int argc, char** argv){
int main(int argc, char** argv){

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
    inverse_covariance_matrix = calloc(lmax+1, sizeof(double *));
    for(ell_value=0; ell_value<lmax+1; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    get_inverse_covariance_matrix_NxN(&S2HAT_params, inverse_covariance_matrix);
    
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
    double power_inv_cov = 0;
    printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);
    apply_inv_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, power_inv_cov, &S2HAT_params);

    
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
