#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
#include <unistd.h>

// #include "s2hat.h"
#include "midapack.h"
// #include "s2hat_tools.h"


int apply_alm2pix(double *local_map_pix, s2hat_dcomplex *local_alm, S2HAT_parameters *S2HAT_params){
    /* Transform alm coefficients local_alm into a pixel map local_map_pix, 
    all details here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Calm2map.html 

    local_alm is a 4-dimensional array in the form :
        (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
        (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
    in the form I, E, B
    */
    
    int nstokes = S2HAT_params->nstokes; // Getting number of Stokes parameters

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params->Local_param_s2hat);
    // Getting global and local S2HAT parameters needed for the computation
    
    if (Local_param_s2hat->gangrank < 0)
        return 0;
    // Only apply the transformation to the processes used with S2HAT

    int nmaps = 1; // We only provide 1 input set of alm coefficient
    // int nstokes = 3; // We want all T, Q and U maps
    // int lda = nstokes; // We choose the HEALPIX convention with local_alm in the form (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps)
    // int lda = Global_param_s2hat->nlmax; // We choose the S2HAT convention with local_alm in the form (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps)
    int lda = S2HAT_params->lda;
    int spin;

    switch(nstokes)
    {
        case 1: // Case only intensity
            printf("<<<< Test 39 -- \n"); fflush(stdout);
            spin=0;
            s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            break;
        case 3: // Case with both intensity and polarization
            printf("<<<< Test 41 -- lmax %d map_size %d lda %d \n", Global_param_s2hat->nlmax, Local_param_s2hat->map_size, lda); fflush(stdout);
            // s2hat_alm2map(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
            //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, 
            //     local_alm, Local_param_s2hat->nplm, NULL, 
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0

            // s2hat_alm2map(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, 1, 
            //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, 
            //     local_alm, Local_param_s2hat->nplm, NULL, 
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);

            spin=0;
            s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
            printf("<<<< Test 41b -- \n"); fflush(stdout);
            spin=2;
            s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm+S2HAT_params->size_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
            break;

        case 2: // Case with only polarization
            printf("<<<< Test 40 -- \n"); fflush(stdout);
            spin=nstokes;
            s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
            printf("<<<< Test 42 -- \n"); fflush(stdout);
            break;
    }
    return 0;
}



int apply_pix2alm(double *local_map_pix, s2hat_dcomplex *local_alm, S2HAT_parameters *S2HAT_params){
    /* Transform pixel map into alm coefficients, all details here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Cmap2alm.html 
     local_alm is a 4-dimensional array in the form :
        (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
        (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
    Here the HEALpix convention has been chosen
    Output will be in the form I, E, B
    */
    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params->Local_param_s2hat);
    // Getting global and local S2HAT parameters needed for the computation
    
    if (Local_param_s2hat->gangrank < 0)
        return 0;
    // Only apply the transformation to the processes used with S2HAT

    int nstokes = S2HAT_params->nstokes; // Getting number of Stokes parameters

    int nmaps = 1; // We only provide 1 input set of alm coefficient
    // int lda = nstokes; // We choose the HEALPIX convention with local_alm in the form (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps)
    // int lda = Global_param_s2hat->nlmax; // We choose the S2HAT convention with local_alm in the form (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps)
    int lda = S2HAT_params->lda;

    int nrings = Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1;
    // double *local_w8ring = (double *)malloc(nrings*nstokes*sizeof(double));
    double local_w8ring[nrings*nstokes];
    int i_ring;

    int spin;

    // if (Local_param_s2hat->gangrank != -1){
    //   local_w8ring = (double *) malloc( nstokes*(Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1)*sizeof(double));
    //   for( i_ring=0; i_ring< nstokes*(Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1); i_ring++)
    //           local_w8ring[i_ring]=1.;
    // }
      for( i_ring=0; i_ring< nstokes*nrings; i_ring++)
              local_w8ring[i_ring]=1.;

    // int max_size_test = 50, number_of_nans = 0, index;
    // if (Local_param_s2hat->gangrank == 0){
    //     printf("%d ²²²²²² local_map_pix from pixel 0 to %d - %f -", Local_param_s2hat->gangrank, max_size_test, local_map_pix[0]);
    //     for (index=0+Local_param_s2hat->map_size;index<max_size_test+Local_param_s2hat->map_size;index++){
    //         printf("- %f -", local_map_pix[index]);
    //         if (!(local_map_pix[index]==local_map_pix[index]))
    //             number_of_nans++;

    //         }
    //     printf(" \n"); fflush(stdout);
    //     printf("%d ²²²²²² number_of_nans %d \n", Local_param_s2hat->gangrank, number_of_nans); fflush(stdout);
    // }

    switch(nstokes)
    {
        case 1: // Case only intensity
            spin=0;            
            s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		        local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            break;
        case 3: // Case with both intensity and polarization
            // s2hat_map2alm(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
            //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm, 
            //     Local_param_s2hat->nplm, NULL,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0
            spin=0;            
            s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		        local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            spin=2;            
            s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		        local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm+S2HAT_params->size_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            break;

        case 2: // Case with only polarization
            spin=nstokes;            
            s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		        local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
            break;
    }

    // free(local_w8ring);
    return 0;
}


int apply_inv_covariance_matrix_to_alm(s2hat_dcomplex *input_local_alm, s2hat_dcomplex *out_local_alm, double **inv_covariance_matrix, S2HAT_parameters *S2HAT_params){
    /* Apply inverse of covariance matrix to input_local_alm */

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params->Local_param_s2hat);


    int ell_value, m_value, index_stokes, line_index, nmvals;
    int nstokes = S2HAT_params->nstokes;
    // int nstokes = 3;
    int lmax = Global_param_s2hat->nlmax +1;

    int index, max_size_test=50;
    printf("%d --- Local_alm during apply inv cov matrix - %f %f -", Local_param_s2hat->gangrank, input_local_alm[0].re, input_local_alm[0].im);
    for (index=1;index<max_size_test;index++){
        printf("- %f %f -", input_local_alm[index].re, input_local_alm[index].im);
    }
    printf(" \n");

    if (Local_param_s2hat->gangrank != -1){
        nmvals = Local_param_s2hat->nmvals; // Total number of m values
        // int *mvals = Local_param_s2hat->mvals; // Values of m the considered processor contain

        double res_real, res_imag;

        if(S2HAT_params->lda == Global_param_s2hat->nlmax){
            printf("~~~~ S2HAT convention !! %d \n", S2HAT_params->lda); fflush(stdout);
            for(ell_value=0; ell_value < lmax; ell_value++){
                for(m_value=0; m_value < nmvals; m_value++){
                    for (index_stokes=0; index_stokes<nstokes; index_stokes++){
                        res_real = 0;
                        res_imag = 0;
                        for (line_index=0; line_index < nstokes; line_index++){
                            res_real += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[line_index*nmvals*(lmax) + m_value*(lmax) + ell_value].re;
                            res_imag += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[line_index*nmvals*(lmax) + m_value*(lmax) + ell_value].im;
                        }
                        out_local_alm[index_stokes*nmvals*(lmax) + m_value*(lmax) + ell_value].re = res_real;
                        out_local_alm[index_stokes*nmvals*(lmax) + m_value*(lmax) + ell_value].im = res_imag;
                    }
                }
            }
        }
        else{
            printf("~~~~ HEALPIX convention !! %d \n", S2HAT_params->lda); fflush(stdout);
            for(ell_value=0; ell_value < lmax; ell_value++){
                for(m_value=0; m_value < nmvals; m_value++){
                    for (index_stokes=0; index_stokes<nstokes; index_stokes++){
                        res_real = 0;
                        res_imag = 0;
                        for (line_index=0; line_index < nstokes; line_index++){
                            // res_real += input_local_alm[line_index*(lmax)*nmvals + ell_value*nmvals + m_value].re * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                            // res_imag += input_local_alm[line_index*(lmax)*nmvals + ell_value*nmvals + m_value].im * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                            
                            // res_real += input_local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + line_index].re * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                            // res_imag += input_local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + line_index].im * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];

                            res_real += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + line_index].re;
                            res_imag += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + line_index].im;
                            // Maybe with [line_index*nstokes + index_stokes] ?
                        }
                        // input_local_alm[index_stokes*(lmax)*nmvals + ell_value*nmvals + m_value].re = res_real;
                        // input_local_alm[index_stokes*(lmax)*nmvals + ell_value*nmvals + m_value].im = res_imag;
                        out_local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + index_stokes].re = res_real;
                        out_local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + index_stokes].im = res_imag;
                    }
                    // Verify it is not applied to part where a_lm not defined !!!
                }
            }
        }
    
    printf("%d --- Local_alm just after apply inv cov matrix - %f %f -", Local_param_s2hat->gangrank, input_local_alm[0].re, input_local_alm[0].im);
    for (index=1;index<max_size_test;index++){
        printf("- %f %f -", input_local_alm[index].re, input_local_alm[index].im);
    }
    printf(" \n");
    // long int index;
    // for(index=0; index<nstokes*(lmax)*m_value; index++){
    //     local_alm[index].re = out_local_alm[index].re;
    //     local_alm[index].im = out_local_alm[index].im;
    // }
    // free(out_local_alm);
    }
    return 0;
}


int alm2cls(s2hat_dcomplex* local_alm, double *c_ell_array, int nspec, S2HAT_parameters *S2HAT_params){
    /* Transform alm to c_ell coefficients
     local_alm is a 4-dimensional array in the form :
        (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
        (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
    Here the HEALpix convention has been chosen

    Output :  c_ell_array in the ordering [0:nlmax,1:nspec] with nspec corresponding to TT, EE, BB, TE, TB, EB
              If nspec = 6; // All 6 spectras TT, EE, BB, TE, TB, EB computed

    Note : 
        nspec == ncomp -- only auto-spectra are computed;
        nspec == ncomp+1 ( ncomp>1) -- all ncomp auto spectra and one cross spectrum (e.g., TE if ncomp == 3) are computed;
        nspec == [ncomp (ncomp+1)]/2 (ncomp > 1) -- all auto and all cross spectra are computed.

    */

    S2HAT_GLOBAL_parameters Global_param_s2hat = S2HAT_params->Global_param_s2hat;
    S2HAT_LOCAL_parameters Local_param_s2hat = S2HAT_params->Local_param_s2hat;
    int nstokes = S2HAT_params->nstokes;

    int lmax = Global_param_s2hat.nlmax;
    int ell= 0;
    
    // int nstokes = 3;
    int nmaps = 1;
    int mapnum = 0;
    int ncomp = nstokes; // Number for alm components (T, E, B)
    
    // int lda = ncomp; // Healpix convention chosen
    int lda = S2HAT_params->lda;

    // int nspec = 6; // All 6 spectras TT, EE, BB, TE, TB, EB computed

    // c_ell_array = (double *) calloc( nstokes*(lmax+1), sizeof(double));
    if (Local_param_s2hat.gangrank == 0)
        collect_cls(nmaps, mapnum, ncomp, lmax, Local_param_s2hat.nmvals, Local_param_s2hat.mvals, lda, 
                local_alm, nspec, c_ell_array, Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);

    return 0;
}
