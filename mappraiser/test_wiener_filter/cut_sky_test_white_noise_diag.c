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
// #include "domain_generalization.h"
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

int main_clean_ver_nside_64_all(int argc, char** argv){
// int main(int argc, char** argv){

    printf("############################################################ START !!! ############################################################ \n");
    // char *path_CMB_map = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/Map_test_band_limited.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_SO.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_cut_latS.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_cut_equator.fits";
    char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_pt_src.fits";
    

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

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int i, j, ell_value;

    // int nside = 512;
    int nside = 64;
    int lmax = 2*nside;
    int npix = 12*nside*nside;
    
    // double *mask_binary=NULL;
    double *mask_binary= (double *)malloc(npix*sizeof(double));
    read_fits_mask(nside, mask_binary, path_mask, 1);
    
    int nstokes = 3; // T, Q, U
    int number_correlations = 4; // TT, EE, BB, TE
    // int nstokes = 2; //Q, U
    // int number_correlations = 2; // EE, BB


    
    int index, index_2;
    int order_matrix = nstokes;

    double *CMB_map_ring = (double *) malloc( 3*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, 3);


    int max_size_pixels = 40;
    printf("%d --- Reading first elems of map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f %f -", CMB_map_ring[i], CMB_map_ring[i+npix], CMB_map_ring[i+2*npix]);
    printf("\n");

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
    

    s2hat_dcomplex *local_alm_s2hat;

    local_alm_s2hat = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);

   

    // double **covariance_matrix;
    // covariance_matrix = calloc(lmax, sizeof(double *));
    // for(ell_value=0; ell_value<lmax; ell_value++){
    //     covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    // }
    // get_covariance_matrix_block_diagonal(c_ell_path, number_correlations, covariance_matrix, &S2HAT_params);

    double **inverse_covariance_matrix;
    inverse_covariance_matrix = malloc(lmax*sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    get_inverse_covariance_matrix_diagonal(&S2HAT_params, inverse_covariance_matrix);

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
    }

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

    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("TEST 9c -- lmax %d \n", lmax);
    s2hat_dcomplex *local_alm_s2hat_inverted;

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
    
    printf("%d ----!!!! ***local_alm_s2hat post inv cov*** size_alm %d ; Number of nans %d \n", rank, S2HAT_params.size_alm, number_of_nan); fflush(stdout);


    apply_inv_block_diag_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, &S2HAT_params);

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

    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);

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

    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);

    double *full_sky_map_2;
    full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
    gather_map(local_map_output, full_sky_map_2, nstokes, &S2HAT_params);

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
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_vtuned_2B_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_0b_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_SO_1c_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_cut_latS_1c_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_cut_equator_1c_%d.fits", i);
            sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_pt_src_1c_%d.fits", i);
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
    // free_covariance_matrix(covariance_matrix, lmax-1);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    // free_harmonic_superstruct(&Harmonic_struct, rank);
    free_s2hat_parameters_struct(&S2HAT_params);

    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}


// int main_test_no_inv(int argc, char** argv){
int main(int argc, char** argv){

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

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int i, j, ell_value;

    // int nside = 512;
    int nside = 64;
    int lmax = 2*nside;
    int npix = 12*nside*nside;
    
    double *mask_binary=NULL;
    // double *mask_binary= (double *)malloc(npix*sizeof(double));
    // read_fits_mask(nside, mask_binary, path_mask, 1);
    
    int nstokes = 3; // T, Q, U
    int number_correlations = 4; // TT, EE, BB, TE
    // int nstokes = 2; //Q, U
    // int number_correlations = 2; // EE, BB


    
    int index, index_2;
    int order_matrix = nstokes;

    double *CMB_map_ring = (double *) malloc( 3*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring, path_CMB_map, 3);


    int max_size_pixels = 40;
    printf("%d --- Reading first elems of map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f %f -", CMB_map_ring[i], CMB_map_ring[i+npix], CMB_map_ring[i+2*npix]);
    printf("\n");

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
    

    s2hat_dcomplex *local_alm_s2hat;

    local_alm_s2hat = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);

   

    // double **covariance_matrix;
    // covariance_matrix = calloc(lmax, sizeof(double *));
    // for(ell_value=0; ell_value<lmax; ell_value++){
    //     covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    // }
    // get_covariance_matrix_block_diagonal(c_ell_path, number_correlations, covariance_matrix, &S2HAT_params);

    double **inverse_covariance_matrix;
    inverse_covariance_matrix = malloc(lmax*sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);
    get_inverse_covariance_matrix_diagonal(&S2HAT_params, inverse_covariance_matrix);

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
    }

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

    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("TEST 9c -- lmax %d \n", lmax);
    s2hat_dcomplex *local_alm_s2hat_inverted;

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
    
    printf("%d ----!!!! ***local_alm_s2hat post inv cov*** size_alm %d ; Number of nans %d \n", rank, S2HAT_params.size_alm, number_of_nan); fflush(stdout);


    apply_inv_block_diag_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, &S2HAT_params);

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

    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);

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
    apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);
    // apply_alm2pix(local_map_output, local_alm_s2hat_inverted, &S2HAT_params);

    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);

    double *full_sky_map_2;
    full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
    gather_map(local_map_output, full_sky_map_2, nstokes, &S2HAT_params);

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
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/mapfile_filtered_newver_vtuned_2B_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_0b_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_SO_1c_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_cut_latS_1c_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_cut_equator_1c_%d.fits", i);
            // sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/cutsky_diag_mapfile_whiten_pt_src_1c_%d.fits", i);
            sprintf(filename_save, "/pscratch/sd/m/mag/WF_work/map_files_test_WF/full_sky_test_0a_%d.fits", i);
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
    // free_covariance_matrix(covariance_matrix, lmax-1);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    // free_harmonic_superstruct(&Harmonic_struct, rank);
    free_s2hat_parameters_struct(&S2HAT_params);

    printf("%d --- Done !!! \n", rank); fflush(stdout);

    return 0;
}

