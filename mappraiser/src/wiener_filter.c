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

void apply_Wiener_filter_pixel(int nside, int lmax, int nstokes, double *CMB_map_ring, double *CMB_map_output, double *c_ells, int number_correlations, double *mask_binary, MPI_comm worldcomm){
    /* Apply Wiener filter in pixel space to CMB map */

    printf("############################################################ START !!! ############################################################ \n");
    int rank, size;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm worldcomm = MPI_COMM_WORLD;

    int i, j, ell_value;

    int npix = 12*nside*nside;
    int index, index_2;
    int order_matrix = nstokes;

    // printf("%d --- Initializing S2HAT_param \n", rank); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);

    double *local_map_pix = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));

    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring, local_map_pix, &S2HAT_params);

    s2hat_dcomplex *local_alm_s2hat;

    local_alm_s2hat = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    // printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);

    double **inverse_covariance_matrix;
    inverse_covariance_matrix = malloc(lmax*sizeof(double *));
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    // printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);

    get_covariance_matrix_block_diagonal(c_ells, inverse_covariance_matrix, &S2HAT_params);

    for(ell_value=0; ell_value<lmax+1; ell_value++){
        get_cholesky_decomposition_inverted(nstokes, inverse_covariance_matrix[ell_value], 'L');
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

    // printf("%d ----!!!! ***local_map_output*** nstokes %d ; map_size %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));

    s2hat_dcomplex *local_alm_s2hat_inverted;

    local_alm_s2hat_inverted = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    // printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);
    apply_inv_block_diag_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, &S2HAT_params);

    // printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    
    apply_alm2pix(local_map_output, local_alm_s2hat_inverted, &S2HAT_params);

    // printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);


    gather_map(local_map_output, CMB_map_output, nstokes, &S2HAT_params);

    // printf("%d --- Free step 2 ! \n", rank); fflush(stdout);
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0){
        // printf("%d --- Free step 1b ! \n", rank); fflush(stdout);
        free(local_alm_s2hat_inverted);
    }
    // printf("%d --- Free step 2b ! \n", rank); fflush(stdout);
    free(local_map_output);
    // printf("%d --- Free step 3 ! \n", rank); fflush(stdout);
    free_covariance_matrix(inverse_covariance_matrix, lmax-1);
    // printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    free_s2hat_parameters_struct(&S2HAT_params);

    // printf("%d --- Done !!! \n", rank); fflush(stdout);
}
