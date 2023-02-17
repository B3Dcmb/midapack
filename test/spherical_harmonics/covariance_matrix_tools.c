#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>

// choose header based on compilation option
#ifdef W_MKL
#include <mkl.h>
#else
#include <lapacke.h>
#endif

// #include "fitsio.h"
#include <unistd.h>
#include "s2hat.h"
#include "midapack.h"
#include "s2hat_tools.h"



int alm2cls(s2hat_dcomplex* local_alm, double *c_ell_array, int nspec, int nstokes, S2HAT_parameters *S2HAT_params){
    /* Transform alm to c_ell coefficients
     local_alm is a 4-dimensional array in the form :
        (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
        (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
    Here the HEALpix convention has been chosen

    Output :  c_ell_array in the ordering [0:nlmax,1:nspec] with nspec corresponding to TT, EE, BB, TE, TB, EB
              If nspec = 6; // All 6 spectras TT, EE, BB, TE, TB, EB computed
    */

    S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters Local_param_s2hat = *(S2HAT_params->Local_param_s2hat);

    int lmax = Global_param_s2hat.nlmax;
    int ell= 0;
    
    // int nstokes = 3;
    int nmaps = 1;
    int mapnum = 0;
    int ncomp = nstokes; // Number for alm components (T, E, B)
    
    int lda = ncomp; // Healpix convention chosen

    // int nspec = 6; // All 6 spectras TT, EE, BB, TE, TB, EB computed

    // c_ell_array = (double *) calloc( nstokes*(lmax+1), sizeof(double));
    // printf("Test - gangrank %d", Local_param_s2hat.gangrank);
    collect_cls(nmaps, mapnum, ncomp, lmax, Local_param_s2hat.nmvals, Local_param_s2hat.mvals, lda, 
                local_alm, nspec, c_ell_array, Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    // printf("Test2 - gangrank %d %f %f", Local_param_s2hat.gangrank, *(c_ell_array), *(c_ell_array+lmax-1));

    return 0;
}

int get_inverse_matrix(int order_matrix, double* matrix_to_be_inverted){
    int errorHandler;
    int pivotArray[order_matrix];

    int lda = order_matrix;
    int lwork = order_matrix*order_matrix;
    double work[lwork];

    dgetrf_(&order_matrix, &order_matrix, matrix_to_be_inverted, &lda, pivotArray, &errorHandler);
    // LU decomposition of matrix_to_be_inverted; give result in matrix_to_be_inverted
    // printf("LU decomposition with dgetrf : %d should be zero\n", errorHandler);

    double result[order_matrix*order_matrix];
    dgetri_(&order_matrix, matrix_to_be_inverted, &lda, pivotArray, work, &lwork, &errorHandler);
    // Inversion of system matrix_to_be_inverted
    // printf("Inversion of matrix with dgetri : %d should be zero\n", errorHandler);
}




int get_covariance_matrix_NxN(char *c_ell_path, int number_correl, int nstokes, double **covariance_matrix_NxN, S2HAT_GLOBAL_parameters Global_param_s2hat){
    /* Read c_ell file to compute covariance matrix

    Number_correl is expected to be :
    - 1 : only TT
    - 2 : EE and BB given in this order
    - 3 : EE, BB and EB
    - 4 : TT, EE, BB and TE are given in this order
    - 6 : TT, EE, BB, TE, TB and EB are given in this order
    
    Output : covariance matrix will be a 2 dim arary with lmax as it first dimension, and maximum 9 for its second dimension to contain :
        TT TE TB
        ET EE EB
        BT BE BB
        in this order, ravelled in 1D
        so covariance_matrix_NxN[lmax][9] with 9 being [TT, TE, TB, ET, EE, EB, BT, BE, BB] (with TE=ET, TB=BT and BE=EB)
        
        If the polarization only is focused (so if number_correl = 2 or 3), then the covariance matrix will be in the form :
        EE EB
        BE BB
        in this order, ravelled in 1D

        Finally, if number_correl is 1, then only the information about the intensity will be contained in the covariance matrix
        
        */
    int lmax = Global_param_s2hat.nlmax;
    int correl_index, ell_value;
    double *c_ell_array;

    if ((number_correl != 5) && (number_correl <= 6)){
        printf("Error : number_correl must be either 1: TT ; 2: EE, BB ; 3: EE, BB, BE  ; 4: TT, EE, BB and TE ; or 6: TT, EE, BB, TE, TB and EB \n");
        fflush(stdout);
    }

    c_ell_array = malloc(number_correl*(lmax+1)*sizeof(double));
    read_fits_cells(lmax+1, number_correl, c_ell_array, c_ell_path, 1); // Reading cell_fits_file


    for (ell_value=0; ell_value<lmax+1; ell_value++){
        for (correl_index=0; correl_index<nstokes; correl_index++){
            covariance_matrix_NxN[ell_value][correl_index*nstokes + correl_index] = c_ell_array[ (lmax+1)*correl_index + ell_value ]; // Diagonal part : TT (0), EE (4), BB (8)
            // printf("Test2 : %d - %d - %f %f \n", ell_value, correl_index, covariance_matrix_NxN[ell_value][correl_index*3 + correl_index], c_ell_array[ (lmax+1)*correl_index + ell_value ]);
        }

        if (nstokes>1){
            covariance_matrix_NxN[ell_value][1] = c_ell_array[ nstokes*(lmax+1) + ell_value ]; // Cross-correlation TE (if case with intensity+polarization) or EB (if only polarization) (up-right block)
            covariance_matrix_NxN[ell_value][nstokes] = c_ell_array[ nstokes*(lmax+1) + ell_value ]; // Cross-correlation TE (if case with intensity+polarization) or EB (if only polarization) (middle-left block)
            if(number_correl == 6){
                covariance_matrix_NxN[ell_value][2] = c_ell_array[ 4*(lmax+1) + ell_value ]; // Cross-correlation TB (up-right block)
                covariance_matrix_NxN[ell_value][6] = c_ell_array[ 4*(lmax+1) + ell_value ]; // Cross-correlation TB (bottom-left block)
                
                covariance_matrix_NxN[ell_value][5] = c_ell_array[ 5*(lmax+1) + ell_value ]; // Cross-correlation EB (bottom-middle block)
                covariance_matrix_NxN[ell_value][7] = c_ell_array[ 5*(lmax+1) + ell_value ]; // Cross-correlation EB (middle-right block)
            }
        }
    }

    free(c_ell_array);
    return 0;
}


int get_inverse_covariance_matrix_NxN(int nstokes, S2HAT_parameters *S2HAT_params, double **inverse_covariance_matrix){
    /* Function to obtain inverse of covariance matrix in harmonic domain, from given c_ells
    */

    Files_path_WIENER_FILTER Files_path_WF_struct = *(S2HAT_params->Files_WF_struct);
    S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);

    double **covariance_matrix;
    int ell_value, index_1;
    int lmax = Global_param_s2hat.nlmax;

    // covariance_matrix = calloc(lmax+1, sizeof(double *));
    // for(ell_value=0; ell_value<lmax+1; ell_value++){
    //     covariance_matrix[ell_value] = calloc(9,sizeof(double));
    // }

    char *c_ell_path = Files_path_WF_struct.c_ell_path;
    int number_correlations = Files_path_WF_struct.number_correlations;
    get_covariance_matrix_NxN(c_ell_path, number_correlations, nstokes, inverse_covariance_matrix, Global_param_s2hat);

    for(ell_value=0; ell_value<lmax+1; ell_value++){
        get_inverse_matrix(nstokes*nstokes, inverse_covariance_matrix[ell_value]);
    }

    // It's possible covariance_matrix will be returned as [3][3], which is not what we want
    // To maybe modify/verify later

    // for (index_1=0; index_1<lmax+1; index_1++){
    //         free(covariance_matrix[index_1]);
    // }
    // free(covariance_matrix);
    return 0;
}


