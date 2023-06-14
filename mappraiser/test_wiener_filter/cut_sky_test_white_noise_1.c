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

void write_fits_vect( int nele, double *vect, char *outfile, int hdu);

void write_fits_vect( int nele, double *vect, char *outfile, int hdu)
{
  fitsfile *fptr;

  int irow,status = 0;
  char newfilename[1000];
  if( hdu == 1) sprintf(newfilename,"!%s",outfile);
  else sprintf(newfilename,"%s",outfile);

  char * coltype[1] ={"VECT"};
  char * colform[1] ={"1D"};
  char *tunit[1] = { "toto"};
  char extname[] = "ARRAY";

  if( hdu == 1) ffinit( &fptr, newfilename, &status);
  else ffopen( &fptr, newfilename, 1, &status);

  fits_create_tbl(fptr, BINARY_TBL, 0, 1, coltype, colform, tunit, extname, &status);

  fits_movabs_hdu( fptr, hdu+1, NULL, &status);
  for(irow=1; irow<=nele; irow++)
  {
    // printf(" // %e %e %e %e %e %e %e %e %e %e %e %e %e %e \n", mll[10], mll[20], mll[30], mll[40], mll[50], mll[60], mll[70], mll[80], mll[90], mll[100], mll[110], mll[120], mll[130], mll[140]);
    fits_write_col( fptr, TDOUBLE, 1, irow, 1, 1, &vect[irow-1], &status);
  }

  ffclos(fptr,&status);

}

// int main_clean_ver_nside_64_cut_sky(int argc, char** argv){
int main(int argc, char** argv){

    printf("############################################################ START !!! ############################################################ \n");
    // char *path_CMB_map = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/Map_test_band_limited.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE_128.fits";
    char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";
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
    int rank, size;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int i, j, ell_value;

    int nside = 64;
    int lmax = 2*nside;

    int nstokes = 3; // T, Q, U
    int number_correlations = 4; // TT, EE, BB, TE
    // int nstokes = 2; //Q, U
    // int number_correlations = 2; // EE, BB


    int npix = 12*nside*nside;
    int index, index_2;
    int order_matrix = nstokes;

    double *CMB_map_ring = (double *) malloc( 3*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, 3);

    int *mask_binary=NULL;
    // printf("%d --- Setting T to 0 !!!! \n", rank);
    // for (i=0; i<npix; i++)
    //     CMB_map_ring[i] = 0;
    int max_size_pixels = 40;
    printf("%d --- Reading first elems of map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f %f -", CMB_map_ring[i], CMB_map_ring[i+npix], CMB_map_ring[i+2*npix]);
    printf("\n");

    // double *CMB_map_polar = CMB_map_ring+npix;

    

    printf("%d --- Initializing S2HAT_param \n", rank); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);
    printf("%d --- End of initializing S2HAT_param \n", rank); fflush(stdout);

    double *local_map_pix = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("%d --- Distributing S2HAT_param \n", rank); fflush(stdout);
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring, local_map_pix, &S2HAT_params);
    // distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_polar, local_map_pix, &S2HAT_params);
    printf("%d --- End of distributing S2HAT_param \n", rank); fflush(stdout);

    // printf("%d --- Reading first elems of local map \n", rank);
    // for (i=0; i<max_size_pixels; i++)
    //     printf("- %f %f -", local_map_pix[i], CMB_map_polar[i]);
    // printf("\n");
    // printf("%d --- Reading first+npix elems of local map \n", rank);
    // for (i=0; i<max_size_pixels; i++)
    //     printf("- %f %f %f -", local_map_pix[i+npix], CMB_map_polar[i+npix], CMB_map_ring[i+2*npix]);
    // printf("\n");

    s2hat_dcomplex *local_alm_s2hat;
    // if (S2HAT_params.Local_param_s2hat.gangrank >= 0)
    // local_alm_s2hat = (s2hat_dcomplex *)malloc(nstokes*S2HAT_params.size_alm*sizeof(s2hat_dcomplex));
    local_alm_s2hat = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);


    // double *c_ell_array_output = (double *)malloc(nstokes*lmax*sizeof(double));
    // int nspec = 3;
    // alm2cls(local_alm_s2hat, c_ell_array_output, nspec, &S2HAT_params);
    // // char save_filename[100] = "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cell_prewitened_newver_woTTTE_vtuned_1b.fits";
    // // char save_filename[100] = "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cell_prewitened_newver_woTE_vtuned_1b.fits";
    // char save_filename[100] = "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cell_prewitened_newver_only_BB_vtuned_1b.fits";
    // // sprintf( save_filename, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cell_prewitened_newver_woTTTE_vtuned_0.fits", imap1, imap2);
    // // save_filename 
	// write_fits_vect( nstokes*lmax, c_ell_array_output, save_filename, 1); 
    // free(c_ell_array_output);
    // printf("%d --- Free step 2 ! \n", rank); fflush(stdout);
    // free(local_map_pix);


    int max_size_test = 50;
    
    // double **covariance_matrix;
    // covariance_matrix = calloc(lmax, sizeof(double *));
    // for(ell_value=0; ell_value<lmax; ell_value++){
    //     covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    // }
    // get_covariance_matrix_NxN(c_ell_path, number_correlations, covariance_matrix, &S2HAT_params);

    double **inverse_covariance_matrix;
    inverse_covariance_matrix = malloc(lmax*sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    get_inverse_covariance_matrix_diagonal(&S2HAT_params, inverse_covariance_matrix);

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

    printf("%d --- printing inverse of covariance matrix !!! lmax %d \n", rank, lmax); fflush(stdout);
    
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

    double normalization_factor = 0; //0;
    for (ell_value=0; ell_value<lmax; ell_value++){
        normalization_factor += (2*ell_value + 1)/(4*M_PI);
    }

    printf("Value of normalization factor %f \n", normalization_factor);
    normalization_factor = normalization_factor;
    char cholesky_part = 'L';
    int info;

    for(ell_value=0; ell_value<lmax; ell_value++){
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                inverse_covariance_matrix[ell_value][index*order_matrix + index_2] *= 1/normalization_factor;
            }
        }
        dpotrf_(&cholesky_part, &nstokes, inverse_covariance_matrix[ell_value], &nstokes, &info);

        if (ell_value<size_ell_max){
            printf("\n #### ell_Cholesky= %d \n", ell_value);
            for (index=0; index < order_matrix; index++){
                for (index_2=0; index_2<order_matrix; index_2++){
                    printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
                }
                printf("\n");
            }
        }
        inverse_covariance_matrix[ell_value][nstokes] = 0;
        if (nstokes == 3){
            inverse_covariance_matrix[ell_value][6] = 0;
            inverse_covariance_matrix[ell_value][7] = 0;
        }
        if (ell_value<size_ell_max){
            printf("\n #### ell_final= %d \n", ell_value);
            for (index=0; index < order_matrix; index++){
                for (index_2=0; index_2<order_matrix; index_2++){
                    printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
                }
                printf("\n");
            }
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
    // // free(transpose_test);

    // // Changing the value of inv_cov_matrix to match exactly the Cholesky value
    // normalization_factor = 1;
    // for(ell_value=1; ell_value<lmax; ell_value++){
    //     // if ((ell_value == 76)||(ell_value == 88)||(ell_value == 112)){
    //     //     printf("Value pre-Cholesky nstokes %d ell %d --- inv_cov_matrix %f \n", nstokes, ell_value,inverse_covariance_matrix[ell_value][0]);
    //     //     printf("Rest : [%f, %f] [%f, %f] \n", inverse_covariance_matrix[ell_value][0], inverse_covariance_matrix[ell_value][1], inverse_covariance_matrix[ell_value][2], inverse_covariance_matrix[ell_value][3]);
    //     // }
    //     // dpotrf_(&cholesky_part, &nstokes, inverse_covariance_matrix[ell_value], &nstokes, &info);
    //     // inverse_covariance_matrix[ell_value][0] = 1/sqrt(covariance_matrix[ell_value][0]); //*(2*ell_value + 1));
    //     // inverse_covariance_matrix[ell_value][1] = 0;
    //     // inverse_covariance_matrix[ell_value][2] = 0;
    //     // inverse_covariance_matrix[ell_value][3] = 1/sqrt(covariance_matrix[ell_value][3]); //*(2*ell_value + 1));
    //     inverse_covariance_matrix[ell_value][0] = 0; //1/sqrt(covariance_matrix[ell_value][0])/sqrt(normalization_factor); //*(2*ell_value + 1));
    //     inverse_covariance_matrix[ell_value][1] = 0;
    //     inverse_covariance_matrix[ell_value][2] = 0;
    //     inverse_covariance_matrix[ell_value][3] = 0; //1/sqrt(covariance_matrix[ell_value][3])/sqrt(normalization_factor); //*(2*ell_value + 1));
    //     inverse_covariance_matrix[ell_value][4] = 1/sqrt(covariance_matrix[ell_value][4])/sqrt(normalization_factor); //*(2*ell_value + 1));
    //     inverse_covariance_matrix[ell_value][8] = 1/sqrt(covariance_matrix[ell_value][8])/sqrt(normalization_factor); //*(2*ell_value + 1));
    //     if ((covariance_matrix[ell_value][4] <= 0)||(covariance_matrix[ell_value][8]<=0)){
    //         printf("NEGATIVE ELEMENT IN COV_MATRIX %d - %f %f !!!", ell_value, covariance_matrix[ell_value][4], covariance_matrix[ell_value][8]);
    //     }
        
    //     // if (nstokes == 3){
    //     //     inverse_covariance_matrix[ell_value][6] = 0;
    //     //     inverse_covariance_matrix[ell_value][7] = 0;
    //     // }
    // }
    // printf("\n");

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
    printf("TEST 9 -- lmax %d \n", lmax);
    printf("%d ----!!!! ***local_map_output*** nstokes %d ; map_size %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("TEST 9c -- lmax %d \n", lmax);
    s2hat_dcomplex *local_alm_s2hat_inverted;
    // if (S2HAT_params.Local_param_s2hat.gangrank >= 0)
    // local_alm_s2hat_inverted = (s2hat_dcomplex *)malloc(S2HAT_params.nstokes*S2HAT_params.size_alm*sizeof(s2hat_dcomplex));
    local_alm_s2hat_inverted = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);

    int number_of_nan = 0;
    for (index=0;index<S2HAT_params.nstokes*S2HAT_params.size_alm;index++){
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

    int nmvals = Local_param_s2hat->nmvals, index_stokes, m_value;
    printf("%d ----!!!! ***local_alm_s2hat pre inv cov*** size_alm %d ; Number of nans %d \n", rank, S2HAT_params.size_alm, number_of_nan); fflush(stdout);
    printf("%d ----!!!! Values nans : ", rank); fflush(stdout);
    ell_value = 1; m_value = 0; index_stokes = 1; 
    printf("- ell %d m_ %d stokes %d alm.re %f alm.im %f -", ell_value, m_value, nstokes, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].re, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].im);
    ell_value = 1; m_value = 1; index_stokes = 0; 
    printf("- ell %d m_ %d stokes %d alm.re %f alm.im %f -", ell_value, m_value, nstokes, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].re, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].im);
    ell_value = 1; m_value = 1; index_stokes = 1; 
    printf("- ell %d m_ %d stokes %d alm.re %f alm.im %f -", ell_value, m_value, nstokes, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].re, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].im);
    ell_value = 1; m_value = 1; index_stokes = 2; 
    printf("- ell %d m_ %d stokes %d alm.re %f alm.im %f -", ell_value, m_value, nstokes, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].re, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].im);
    ell_value = 1; m_value = 3; index_stokes = 2; 
    printf("- ell %d m_ %d stokes %d alm.re %f alm.im %f -", ell_value, m_value, nstokes, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].re, local_alm_s2hat[nstokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].im);


    printf("\n");

    number_of_nan = 0;
    printf("%d ----!!!! ***inverse_covariance_matrix*** lmax %d number_of_nan reset to %d  \n", rank, lmax, number_of_nan); fflush(stdout);
    for (ell_value=0;ell_value<lmax;ell_value++){
        
        // printf("\n #### ell= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                if (!(inverse_covariance_matrix[ell_value][index*nstokes+index_2] == inverse_covariance_matrix[ell_value][index*nstokes+index_2])){
                    // printf(" --- NAN HERE : ell %d index %d index_2 %d value  EE %f BB %.10f -- ", ell_value, index, index_2, covariance_matrix[ell_value][0], covariance_matrix[ell_value][3]);
                    number_of_nan++; 
                }
            }
        }
        // printf(" --- number_of_nan inv cov : %d \n", number_of_nan);
    }
    printf(" --- number_of_nan of inv cov : %d \n", number_of_nan);
    
    // printf("TEST ?! \n");
    // printf(" \n"); fflush(stdout);
    printf("%d ----!!!! ***local_alm_s2hat post inv cov*** size_alm %d ; Number of nans %d \n", rank, S2HAT_params.size_alm, number_of_nan); fflush(stdout);

    // printf("%d ----!!!! ***local_alm_s2hattest*** size_alm %d : \n", rank, S2HAT_params.size_alm); fflush(stdout);
    // for (ell_value=0; ell_value<2; ell_value++){
    //     for(m_value=0; m_value<2*ell_value+1)
    // }
    // input_local_alm[line_index*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value].re
    apply_inv_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, &S2HAT_params);

    number_of_nan = 0;
    for (index=0;index<S2HAT_params.nstokes*S2HAT_params.size_alm;index++){
        if (!(local_alm_s2hat_inverted[index].re == local_alm_s2hat_inverted[index].re)){
            // local_alm_s2hat[index].re = 0;
            number_of_nan++;
            }
        if (!(local_alm_s2hat_inverted[index].im == local_alm_s2hat_inverted[index].im)){
            // local_alm_s2hat[index].im = 0;
            number_of_nan++;
            }
    }
    printf(" \n"); fflush(stdout);
    printf("%d ----!!!! ***local_alm_s2hat_inverted post inv cov*** size_alm %d ; Number of nans %d \n", rank, S2HAT_params.size_alm, number_of_nan); fflush(stdout);
    
    // double *c_ell_array_output_2 = (double *)malloc(nstokes*lmax*sizeof(double));
    // nspec = 3;
    // alm2cls(local_alm_s2hat_inverted, c_ell_array_output_2, nspec, &S2HAT_params);
    // // char save_filename_2[100] = "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cell_witened_newver_woTTTE_vtuned_1b.fits";
    // // char save_filename_2[100] = "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cell_witened_newver_vtuned_1b.fits";
    // // char save_filename_2[100] = "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cell_witened_newver_WOTE_vtuned_1b.fits";
    // char save_filename_2[100] = "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cell_witened_newver_only_BB_vtuned_1b.fits";

	// write_fits_vect( nstokes*lmax, c_ell_array_output_2, save_filename_2, 1); 
    // free(c_ell_array_output_2);

    number_of_nan = 0;
    for (index=0;index<S2HAT_params.nstokes*S2HAT_params.size_alm;index++){
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
    printf("%d ----!!!! ***local_alm_s2hat post inv cov*** size_alm %d ; Number of nans %d \n", rank, S2HAT_params.size_alm, number_of_nan); fflush(stdout);

    number_of_nan = 0;
    for (index=0;index<S2HAT_params.nstokes*S2HAT_params.size_alm;index++){
        if (!(local_alm_s2hat_inverted[index].re == local_alm_s2hat_inverted[index].re)){
            // local_alm_s2hat[index].re = 0;
            number_of_nan++;
            }
        if (!(local_alm_s2hat_inverted[index].im == local_alm_s2hat_inverted[index].im)){
            // local_alm_s2hat[index].im = 0;
            number_of_nan++;
            }
    }
    printf(" \n"); fflush(stdout);
    printf("%d ----!!!! ***local_alm_s2hat_inverted post alm2cls*** size_alm %d ; Number of nans %d \n", rank, S2HAT_params.size_alm, number_of_nan); fflush(stdout);
    
    // printf("%d --- Free step 1a ! \n", rank); fflush(stdout);
    // free(local_alm_s2hat);
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

    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    // local_map_MAPPRAISER_output = (double *)malloc(Obj_pointing_matrix.lcount*sizeof(double));
    // free(CMB_map_ring);
    // printf("%d --- harmonic2map Test 0 !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    // double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    // double *local_map_output = (double *)malloc(2*49152*sizeof(double));
    // printf("%d --- harmonic2map Test !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    // global_harmonic_2_map(local_map_MAPPRAISER_output, local_alm_s2hat_inverted, &Obj_pointing_matrix, &(Harmonic_struct));
    // apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);
    number_of_nan = 0;
    for (index=0;index<S2HAT_params.nstokes*S2HAT_params.size_alm;index++){
        if (!(local_alm_s2hat_inverted[index].re == local_alm_s2hat_inverted[index].re)){
            // local_alm_s2hat[index].re = 0;
            number_of_nan++;
            }
        if (!(local_alm_s2hat_inverted[index].im == local_alm_s2hat_inverted[index].im)){
            // local_alm_s2hat[index].im = 0;
            number_of_nan++;
            }
    }
    printf(" \n"); fflush(stdout);

    printf("%d ----!!!! ***local_alm_s2hat_inverted*** size_alm %d ; Number of nans %d \n", rank, S2HAT_params.size_alm, number_of_nan); fflush(stdout);
    
    printf("%d ----TEST 9b \n", rank); fflush(stdout);
    // double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("TEST 10 -- lmax %d \n", lmax);
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

    // printf("%d --- Reading first elems of local map output \n", rank);
    // for (i=0; i<max_size_pixels; i++)
    //     printf("- %f %f -", local_map_output[i], CMB_map_polar[i]);
    // printf("\n");

    double *full_sky_map_2;
    full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
    gather_map(local_map_output, full_sky_map_2, nstokes, &S2HAT_params);

    // printf("%d --- Reading first elems of final maps retrieved \n", rank);
    // for (i=0; i<max_size_pixels; i++)
    //     printf("- %f %f -", full_sky_map_2[i], CMB_map_polar[i]);
    // printf("\n");

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
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_woTTTE_v2_1_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_woTTTE_vtuned_5b_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_vtuned_1b_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_vtuned_woTE_1d_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_vtuned_only_BB_1a_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_vtuned_only_EEBB_1e_%d.fits", i);
            sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_vtuned_2B_%d.fits", i);
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
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0){
        printf("%d --- Free step 1b ! \n", rank); fflush(stdout);
        free(local_alm_s2hat_inverted);
        // printf("%d --- Free step 1a ! \n", rank); fflush(stdout);
        // free(local_alm_s2hat);
    }
    printf("%d --- Free step 2b ! \n", rank); fflush(stdout);
    free(local_map_output);
    // free(c_ell_array_output_2);
    // free(c_ell_array_output);
    printf("%d --- Free step 3 ! \n", rank); fflush(stdout);
    free_covariance_matrix(inverse_covariance_matrix, lmax-1);
    free_covariance_matrix(covariance_matrix, lmax-1);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    // free_harmonic_superstruct(&Harmonic_struct, rank);
    free_s2hat_parameters_struct(&S2HAT_params);

    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}

