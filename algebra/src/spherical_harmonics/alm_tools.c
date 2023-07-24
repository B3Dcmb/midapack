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

    s2hat_dcomplex *local_alm_T;
    int elem_alm;
    switch(nstokes)
    {
        case 1: // Case only intensity
            // printf("<<<< Test 39 -- \n"); fflush(stdout);
            spin=0;
            s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            break;
        case 3: // Case with both intensity and polarization
            // printf("<<<< Test 41 -- lmax %d map_size %d lda %d \n", Global_param_s2hat->nlmax, Local_param_s2hat->map_size, lda); fflush(stdout);
            // s2hat_alm2map(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
            //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, 
            //     local_alm, Local_param_s2hat->nplm, NULL, 
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0

            // if (S2HAT_params->iter_alm){
            //     s2hat_alm2map(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
            //         Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, 1, 
            //         Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, 
            //         local_alm, Local_param_s2hat->nplm, NULL, 
            //         Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // }
            // else {
            local_alm_T = (s2hat_dcomplex *)malloc(2*S2HAT_params->size_alm*sizeof(s2hat_dcomplex));
            for (elem_alm=0;elem_alm<S2HAT_params->size_alm; elem_alm++){
                local_alm_T[elem_alm].re = -local_alm[elem_alm].re;
                local_alm_T[elem_alm].im = -local_alm[elem_alm].im;
                local_alm_T[elem_alm + S2HAT_params->size_alm].re = 0;
                local_alm_T[elem_alm + S2HAT_params->size_alm].im = 0;
            }
            spin=0;
            s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm_T,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            free(local_alm_T);
            // }
            // Computing alm2map in the specific case where only Q,U are provided
            // printf("<<<< Test 41b -- \n"); fflush(stdout);
            spin=2;
            s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, 
                local_alm+S2HAT_params->size_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
            break;

        case 2: // Case with only polarization
            // printf("<<<< Test 40 -- \n"); fflush(stdout);
            spin=nstokes;
            s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
            // printf("<<<< Test 42 -- \n"); fflush(stdout);
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
    // printf("~~~~ Intermediate step \n"); fflush(stdout);
    int spin=0;

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
    // printf("~~~~ Starting map2alm -- nstokes %d \n", nstokes); fflush(stdout);
    switch(nstokes)
    {
        case 1: // Case only intensity
            // printf("~~~~ Starting map2alm case 1 \n"); fflush(stdout);
            spin=0;            
            s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		        local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // spin=2;
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm+S2HAT_params->size_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            break;
        case 3: // Case with both intensity and polarization
            // s2hat_map2alm(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
            //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm, 
            //     Local_param_s2hat->nplm, NULL,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0
            
            // printf("~~~~ Starting map2alm case 2 \n"); fflush(stdout);
            
            
            // double *local_map_pix_T = (double *)calloc(3*Local_param_s2hat->map_size,sizeof(double));
            // memcpy(local_map_pix_T, local_map_pix, Local_param_s2hat->map_size*sizeof(double));
            // double *local_map_pix_QU = (double *)calloc(3*Local_param_s2hat->map_size,sizeof(double));
            // memcpy(local_map_pix_QU, local_map_pix+Local_param_s2hat->map_size, 2*Local_param_s2hat->map_size*sizeof(double));
            
            // s2hat_dcomplex *local_alm_copy = (s2hat_dcomplex *) calloc( ((2*Global_param_s2hat->nlmax*Global_param_s2hat->nmmax)),sizeof(s2hat_dcomplex));
            spin=0;
            // printf("~~~~ Starting map2alm 1/2 -- spin %d \n", spin); fflush(stdout);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix_T, lda, local_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // free(local_map_pix_T);
            // s2hat_dcomplex *local_alm_T = (s2hat_dcomplex *) calloc( (2*S2HAT_params->size_alm),sizeof(s2hat_dcomplex));
            // s2hat_dcomplex *local_alm_P = (s2hat_dcomplex *) calloc( (2*S2HAT_params->size_alm),sizeof(s2hat_dcomplex));
            s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		        local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm_T,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            spin=2;
            // printf("~~~~ Starting map2alm 2/2 -- spin %d \n", spin); fflush(stdout);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix_QU, lda, local_alm+S2HAT_params->size_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		        local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm+S2HAT_params->size_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm_P,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // memcpy(local_alm, local_alm_T, S2HAT_params->size_alm*sizeof(s2hat_dcomplex));
            // memcpy(local_alm+S2HAT_params->size_alm, local_alm_P, 2*S2HAT_params->size_alm*sizeof(s2hat_dcomplex));
            // printf("~~~~ Finishing map2alm 2/2 -- spin %d \n", spin); fflush(stdout);
            // free(local_map_pix_QU);
            // free(local_alm_T);
            // free(local_alm_P);s

            // printf("~~~~ Starting map2alm case 2 \n"); fflush(stdout);
            // spin=0;            
            // printf("~~~~ Starting map2alm 1/2 -- spin %d \n", spin); fflush(stdout);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // spin=2;            
            // printf("~~~~ Starting map2alm 2/2 -- spin %d \n", spin); fflush(stdout);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm+S2HAT_params->size_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // printf("~~~~ Finishing map2alm 2/2 -- spin %d \n", spin); fflush(stdout);
            break;

        case 2: // Case with only polarization
            // printf("~~~~ Starting map2alm case 3 \n"); fflush(stdout);
            spin=nstokes;            
            s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		        local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
            break;
    }
    // printf("~~~~ Finish map2alm \n"); fflush(stdout);
    // free(local_w8ring);
    return 0;
}


int apply_pix2alm_iter(double *local_map_pix, s2hat_dcomplex *local_alm, S2HAT_parameters *S2HAT_params){
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

    // int niter_max = 1000;
    // double epsilon = 1.0e-10;
    int niter_max = S2HAT_params->iter_alm;
    double epsilon = S2HAT_params->error_alm;
    
    if (Local_param_s2hat->gangrank < 0)
        return 0;
    // Only apply the transformation to the processes used with S2HAT

    int nstokes = S2HAT_params->nstokes; // Getting number of Stokes parameters

    int nmaps = 1; // We only provide 1 input set of alm coefficient
    int lda = S2HAT_params->lda;

    // s2hat_int4 *iter_out = (s2hat_int4 *)calloc( nmaps, sizeof(s2hat_int4));
    // s2hat_flt8 *eps_out  = (s2hat_flt8 *)calloc( nmaps, sizeof(s2hat_flt8)); 
    // int *iter_out = (int *)calloc( nmaps, sizeof(int));
    // float *eps_out  = (float *)calloc( nmaps, sizeof(float)); 
    // int *iter_out_2 = (int *)calloc( nmaps, sizeof(int));
    // float *eps_out_2  = (float *)calloc( nmaps, sizeof(float)); 
    int iter_out[nmaps];
    float eps_out[nmaps];
    int iter_out_2[nmaps];
    float eps_out_2[nmaps];
    // int iter_out = 0;
    // float eps_out = 0;

    int nrings = Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1;
    // double *local_w8ring = (double *)malloc(nrings*nstokes*sizeof(double));
    double local_w8ring[nrings*nstokes];
    int i_ring;
    
    int spin=0;

    for( i_ring=0; i_ring< nstokes*nrings; i_ring++)
            local_w8ring[i_ring]=1.;
    
    int pixel_T;
    double *local_map_pix_T, *local_map_pix_P;
    s2hat_dcomplex *local_alm_T, *local_alm_P;
    // printf("~~~~ Starting map2alm -- nstokes %d \n", nstokes); fflush(stdout);
    switch(nstokes)
    {
        case 1: // Case only intensity
            // printf("~~~~ Starting map2alm case 1 \n"); fflush(stdout);
            spin=0;            
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);

            // s2hat_map2alm_cg( Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //         Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
            //         local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda,
            //         local_alm, nplm, NULL, Local_param_s2hat->gangsize, Local_param_s2hat->gangrank,
            //         niter_max, epsilon, &iter_out, &eps_out, /* define the convergence */
            //         Local_param_s2hat->gangcomm);
            s2hat_map2alm_spin_cg( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                    Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
	                local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm, 
                    Local_param_s2hat->gangsize, Local_param_s2hat->gangrank,
                    niter_max, epsilon, iter_out, eps_out, /* define the convergence */
                    Local_param_s2hat->gangcomm);
            break;
        case 3: // Case with both intensity and polarization
            // s2hat_map2alm(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
            //     Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm, 
            //     Local_param_s2hat->nplm, NULL,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0
            
            // printf("~~~~ Starting map2alm case 2 \n"); fflush(stdout);
            
            
            // double *local_map_pix_T = (double *)calloc(3*Local_param_s2hat->map_size,sizeof(double));
            // memcpy(local_map_pix_T, local_map_pix, Local_param_s2hat->map_size*sizeof(double));
            // double *local_map_pix_QU = (double *)calloc(3*Local_param_s2hat->map_size,sizeof(double));
            // memcpy(local_map_pix_QU, local_map_pix+Local_param_s2hat->map_size, 2*Local_param_s2hat->map_size*sizeof(double));
            
            // s2hat_dcomplex *local_alm_copy = (s2hat_dcomplex *) calloc( ((2*Global_param_s2hat->nlmax*Global_param_s2hat->nmmax)),sizeof(s2hat_dcomplex));
            spin=0;
            // printf("~~~~ Starting map2alm 1/2 -- spin %d \n", spin); fflush(stdout);
            // // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            // //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    // //     local_w8ring, Local_param_s2hat->map_size, local_map_pix_T, lda, local_alm,
            // //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // // free(local_map_pix_T);
            // // s2hat_dcomplex *local_alm_T = (s2hat_dcomplex *) calloc( (2*S2HAT_params->size_alm),sizeof(s2hat_dcomplex));
            // // s2hat_dcomplex *local_alm_P = (s2hat_dcomplex *) calloc( (2*S2HAT_params->size_alm),sizeof(s2hat_dcomplex));
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            // //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    // //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm_T,
            // //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // s2hat_map2alm_cg( Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
            //         Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, 1, 
            //         Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda,
	        //         local_alm, Local_param_s2hat->nplm, NULL, 
            //         Local_param_s2hat->gangsize, Local_param_s2hat->gangrank,
		    //         niter_max, epsilon, iter_out, eps_out,                             /* define the convergence */
		    //         Local_param_s2hat->gangcomm);
            local_map_pix_T = (double *) calloc(2*Local_param_s2hat->map_size, sizeof(double));
            for (pixel_T=0; pixel_T<Local_param_s2hat->map_size; pixel_T++){
                local_map_pix_T[pixel_T] = -local_map_pix[pixel_T];
                // local_map_pix_T[pixel_T+Local_param_s2hat->map_size] = 0;
                // local_map_pix_T[pixel_T+Local_param_s2hat->map_size] = local_map_pix[pixel_T];
            }
            local_alm_T = (s2hat_dcomplex *)calloc(2*S2HAT_params->size_alm,sizeof(s2hat_dcomplex));
            // for (elem_alm=0;elem_alm<2*S2HAT_params->size_alm; elem_alm++){
            //     local_alm_T[elem_alm].re = 0;
            //     local_alm_T[elem_alm].im = 0;
            // }
            s2hat_map2alm_spin_cg( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                    Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
	                local_w8ring, Local_param_s2hat->map_size, local_map_pix_T, lda, local_alm_T, 
                    Local_param_s2hat->gangsize, Local_param_s2hat->gangrank,
                    niter_max, epsilon, iter_out, eps_out, /* define the convergence */
                    Local_param_s2hat->gangcomm);
            free(local_map_pix_T);
            memcpy(local_alm, local_alm_T, S2HAT_params->size_alm*sizeof(s2hat_dcomplex));
            free(local_alm_T);

            spin=2;
            local_map_pix_P = (double *) calloc(2*Local_param_s2hat->map_size, sizeof(double));
            for (pixel_T=0; pixel_T<Local_param_s2hat->map_size; pixel_T++){
                local_map_pix_P[pixel_T] = local_map_pix[pixel_T+Local_param_s2hat->map_size];
                local_map_pix_P[pixel_T+Local_param_s2hat->map_size] = local_map_pix[pixel_T+2*Local_param_s2hat->map_size];
                // local_map_pix_T[pixel_T+Local_param_s2hat->map_size] = 0;
                // local_map_pix_T[pixel_T+Local_param_s2hat->map_size] = -local_map_pix[pixel_T];
            }
            // printf("~~~~ Starting map2alm 2/2 -- spin %d \n", spin); fflush(stdout);
            // // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            // //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    // //     local_w8ring, Local_param_s2hat->map_size, local_map_pix_QU, lda, local_alm+S2HAT_params->size_alm,
            // //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm+S2HAT_params->size_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            // //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    // //     local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm_P,
            // //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // s2hat_map2alm_spin_cg( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //         Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
	        //         local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm+S2HAT_params->size_alm, 
            //         Local_param_s2hat->gangsize, Local_param_s2hat->gangrank,
            //         niter_max, epsilon, iter_out_2, eps_out_2, /* define the convergence */
            //         Local_param_s2hat->gangcomm);
            s2hat_map2alm_spin_cg( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                    Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
	                local_w8ring, Local_param_s2hat->map_size, local_map_pix_P, lda, local_alm+S2HAT_params->size_alm, 
                    Local_param_s2hat->gangsize, Local_param_s2hat->gangrank,
                    niter_max, epsilon, iter_out_2, eps_out_2, /* define the convergence */
                    Local_param_s2hat->gangcomm);
            free(local_map_pix_P);
            // memcpy(local_alm, local_alm_T, S2HAT_params->size_alm*sizeof(s2hat_dcomplex));
            // memcpy(local_alm+S2HAT_params->size_alm, local_alm_P, 2*S2HAT_params->size_alm*sizeof(s2hat_dcomplex));
            // printf("~~~~ Finishing map2alm 2/2 -- spin %d \n", spin); fflush(stdout);
            // free(local_map_pix_QU);
            // free(local_alm_T);
            // free(local_alm_P);s

            // printf("~~~~ Starting map2alm case 2 \n"); fflush(stdout);
            // spin=0;            
            // printf("~~~~ Starting map2alm 1/2 -- spin %d \n", spin); fflush(stdout);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // spin=2;            
            // printf("~~~~ Starting map2alm 2/2 -- spin %d \n", spin); fflush(stdout);
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix+Local_param_s2hat->map_size, lda, local_alm+S2HAT_params->size_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // printf("~~~~ Finishing map2alm 2/2 -- spin %d \n", spin); fflush(stdout);
            break;

        case 2: // Case with only polarization
            // printf("~~~~ Starting map2alm case 3 \n"); fflush(stdout);
            spin=nstokes;            
            // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
            //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
		    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
            //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            s2hat_map2alm_spin_cg( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
                    Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
	                local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm, 
                    Local_param_s2hat->gangsize, Local_param_s2hat->gangrank,
                    niter_max, epsilon, &iter_out, &eps_out, /* define the convergence */
                    Local_param_s2hat->gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
            break;
    }
    // printf("~~~~ Finish map2alm \n"); fflush(stdout);
    // free(local_w8ring);
    return 0;
}


int apply_inv_block_diag_covariance_matrix_to_alm(s2hat_dcomplex *input_local_alm, s2hat_dcomplex *out_local_alm, double **inv_covariance_matrix, S2HAT_parameters *S2HAT_params){
    /* Apply inverse of covariance matrix to input_local_alm */

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params->Local_param_s2hat);


    int ell_value, m_value, index_stokes, line_index, nmvals;
    int nstokes = S2HAT_params->nstokes;
    int lmax = Global_param_s2hat->nlmax;// +1;

    int number_of_nan = 0, number_of_nan_re = 0, number_of_nan_im = 0;
    int number_of_nan_inv_cov = 0;
    int number_of_nan_alm = 0;
    if (Local_param_s2hat->gangrank != -1){
        nmvals = Local_param_s2hat->nmvals; // Total number of m values
        // int *mvals = Local_param_s2hat->mvals; // Values of m the considered processor contain

        double res_real, res_imag;

        if(S2HAT_params->lda == Global_param_s2hat->nlmax){
            // printf("~~~~ S2HAT convention !! lda %d lmax %d nmvals %d nmax %d \n", S2HAT_params->lda, lmax, nmvals, Global_param_s2hat->nmmax); fflush(stdout);
            int element_alm;
            for (element_alm=0; element_alm<nstokes*S2HAT_params->size_alm; element_alm++){
                out_local_alm[element_alm].re = 0;
                out_local_alm[element_alm].im = 0;

            }

            int max_value_ell = 3;
            for(ell_value=0; ell_value < lmax+1; ell_value++){
                for(m_value=0; m_value < nmvals; m_value++){
                    // if ((ell_value < max_value_ell) && (m_value < 2*ell_value + 2)){
                    //     printf(" Multiply ell %d m_value %d :\n", ell_value, m_value); 
                    // }
                // for(m_value=0; m_value < min(2*(ell_value+1)+1, Global_param_s2hat->nmmax); m_value++){
                    // for (index_stokes=0; index_stokes<nstokes; index_stokes++){
                    //     res_real = 0;
                    //     res_imag = 0;
                    //     for (line_index=0; line_index<nstokes; line_index++){
                    //         res_real += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].re;
                    //         res_imag += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].im;
                    //         // res_real += inv_covariance_matrix[ell_value][line_index*nstokes + index_stokes] * input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].re;
                    //         // res_imag += inv_covariance_matrix[ell_value][line_index*nstokes + index_stokes] * input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].im;
                    //     }
                    //     out_local_alm[(index_stokes*nmvals + m_value)*(lmax+1) + ell_value].re = res_real;
                    //     out_local_alm[(index_stokes*nmvals + m_value)*(lmax+1) + ell_value].im = res_imag;
                    // }
                    for (index_stokes=0; index_stokes<nstokes; index_stokes++){
                        // res_real = 0;
                        // res_imag = 0;
                        // if ((ell_value < max_value_ell) && (m_value < 2*ell_value + 2)){
                        //         printf("# nstokes %d :", index_stokes); 
                        //     }
                        for (line_index=0; line_index<nstokes; line_index++){
                            out_local_alm[(index_stokes*nmvals + m_value)*(lmax+1) + ell_value].re += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].re;
                            out_local_alm[(index_stokes*nmvals + m_value)*(lmax+1) + ell_value].im += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].im;
                            // res_real += inv_covariance_matrix[ell_value][line_index*nstokes + index_stokes] * input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].re;
                            // res_imag += inv_covariance_matrix[ell_value][line_index*nstokes + index_stokes] * input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].im;
                            // if ((ell_value < max_value_ell) &&   (m_value < 2*ell_value + 2)){
                            //     printf("- %f  re %f im %f -", inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index], input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].re, input_local_alm[(line_index*nmvals + m_value)*(lmax+1) + ell_value].im); 
                            // }
                            
                        }
                        // if ((ell_value < max_value_ell) &&   (m_value < 2*ell_value + 2)){
                        //         printf("\n"); 
                        // }
                            
                        // out_local_alm[(index_stokes*nmvals + m_value)*(lmax+1) + ell_value].re = res_real;
                        // out_local_alm[(index_stokes*nmvals + m_value)*(lmax+1) + ell_value].im = res_imag;
                    }
                }
            }
        }
        else{
            printf("~~~~ HEALPIX convention !! %d \n", S2HAT_params->lda); fflush(stdout);
            for(ell_value=0; ell_value < lmax+1; ell_value++){
                for(m_value=0; m_value < nmvals; m_value++){
                    for (index_stokes=0; index_stokes<nstokes; index_stokes++){
                        res_real = 0;
                        res_imag = 0;
                        for (line_index=0; line_index < nstokes; line_index++){
                            // res_real += input_local_alm[line_index*(lmax+1)*nmvals + ell_value*nmvals + m_value].re * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                            // res_imag += input_local_alm[line_index*(lmax+1)*nmvals + ell_value*nmvals + m_value].im * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                            
                            // res_real += input_local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + line_index].re * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                            // res_imag += input_local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + line_index].im * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];

                            res_real += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + line_index].re;
                            res_imag += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * input_local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + line_index].im;
                            // Maybe with [line_index*nstokes + index_stokes] ?
                        }
                        // input_local_alm[index_stokes*(lmax+1)*nmvals + ell_value*nmvals + m_value].re = res_real;
                        // input_local_alm[index_stokes*(lmax+1)*nmvals + ell_value*nmvals + m_value].im = res_imag;
                        out_local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + index_stokes].re = res_real;
                        out_local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + index_stokes].im = res_imag;
                    }
                    // Verify it is not applied to part where a_lm not defined !!!
                }
            }
        }
    }
    return 0;
}


// int apply_inv_full_covariance_matrix_to_alm(s2hat_dcomplex *input_local_alm, s2hat_dcomplex *out_local_alm, double **inv_covariance_matrix, S2HAT_parameters *S2HAT_params)
// {
//     /* Apply inverse of covariance matrix to input_local_alm */

//     S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params->Global_param_s2hat);
//     S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params->Local_param_s2hat);


//     int ell_value_1, ell_value_2, m_value, index_stokes, line_index, nmvals;
//     int nstokes = S2HAT_params->nstokes;
//     // int nstokes = 3;
//     int lmax = Global_param_s2hat->nlmax;// +1;

//     int index, max_size_test=50;
//     printf("%d --- Local_alm during apply inv cov matrix - %f %f -", Local_param_s2hat->gangrank, input_local_alm[0].re, input_local_alm[0].im);
//     for (index=1;index<max_size_test;index++){
//         printf("- %f %f -", input_local_alm[index].re, input_local_alm[index].im);
//     }
//     printf(" \n");

//     int number_of_nan = 0, number_of_nan_re = 0, number_of_nan_im = 0;
//     int number_of_nan_inv_cov = 0;
//     int number_of_nan_alm = 0;
//     if (Local_param_s2hat->gangrank != -1){
//         nmvals = Local_param_s2hat->nmvals; // Total number of m values
//         // int *mvals = Local_param_s2hat->mvals; // Values of m the considered processor contain

//         double res_real, res_imag;
//         if(S2HAT_params->lda == Global_param_s2hat->nlmax){
//             printf("~~~~ S2HAT convention !! lda %d lmax %d nmvals %d nmax %d \n", S2HAT_params->lda, lmax, nmvals, Global_param_s2hat->nmmax); fflush(stdout);
//             for(ell_value_1=0; ell_value_1<lmax; ell_value_1++){
//                 for(m_value=0; m_value < min(2*(ell_value_1+1)+1, Global_param_s2hat->nmmax); m_value++){
//                     for (index_stokes=0; index_stokes<nstokes; index_stokes++){
//                         res_real = 0;
//                         res_imag = 0;
//                         for (line_index=0; line_index < nstokes; line_index++){
//                             for(ell_value_2=0; ell_value_2<lmax; ell_value_2++){
//                                 res_real += inv_covariance_matrix[ell_value_1*(lmax+1) + ell_value_2][index_stokes*nstokes + line_index] * input_local_alm[line_index*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value_2].re;
//                                 res_imag += inv_covariance_matrix[ell_value_1*(lmax+1) + ell_value_2][index_stokes*nstokes + line_index] * input_local_alm[line_index*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value_2].im;
//                             }
//                             out_local_alm[index_stokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value_1].re = res_real;
//                             out_local_alm[index_stokes*nmvals*(lmax+1) + m_value*(lmax+1) + ell_value_1].im = res_imag;
//                         }
//                     }
//                 }
//             }
//         }
//         else{ /// NOT TESTED BELOW !!!
//             printf("~~~~ HEALPIX convention !! %d \n", S2HAT_params->lda); fflush(stdout);
//             for(ell_value_1=0; ell_value_1 < lmax-1; ell_value_1++){
//                 for(m_value=0; m_value < nmvals; m_value++){
//                     for (index_stokes=0; index_stokes<nstokes; index_stokes++){
//                         res_real = 0;
//                         res_imag = 0;
//                         for (line_index=0; line_index < nstokes; line_index++){
//                             for(ell_value_2=0; ell_value_2<lmax; ell_value_2++){
//                             res_real += inv_covariance_matrix[ell_value_1*(lmax+1) + ell_value_2][index_stokes*nstokes + line_index] * input_local_alm[m_value*(lmax+1)*nstokes + ell_value_2*nstokes + line_index].re;
//                             res_imag += inv_covariance_matrix[ell_value_1*(lmax+1) + ell_value_2][index_stokes*nstokes + line_index] * input_local_alm[m_value*(lmax+1)*nstokes + ell_value_2*nstokes + line_index].im;
//                             // Maybe with [line_index*nstokes + index_stokes] ?
//                             }
//                         }
//                         // input_local_alm[index_stokes*(lmax)*nmvals + ell_value*nmvals + m_value].re = res_real;
//                         // input_local_alm[index_stokes*(lmax)*nmvals + ell_value*nmvals + m_value].im = res_imag;
//                         out_local_alm[m_value*(lmax+1)*nstokes + ell_value_1*nstokes + index_stokes].re = res_real;
//                         out_local_alm[m_value*(lmax+1)*nstokes + ell_value_1*nstokes + index_stokes].im = res_imag;
//                     }
//                     // Verify it is not applied to part where a_lm not defined !!!
//                 }
//             }
//         }
//     // printf("----!!!! Calculation Number of nans *** size_alm %d ; Number of nans re %d ; Number of nans im %d ; Number of nans inv_cov %d ; Number of nans alms %d \n", S2HAT_params->size_alm, number_of_nan_im, number_of_nan_inv_cov, number_of_nan_alm); fflush(stdout);
//     // printf("%d --- Local_alm just after apply inv cov matrix - %f %f -", Local_param_s2hat->gangrank, input_local_alm[0].re, input_local_alm[0].im);
//     // for (index=1;index<max_size_test;index++){
//     //     printf("- %f %f -", input_local_alm[index].re, input_local_alm[index].im);
//     // }
//     // printf(" \n");
//     // number_of_nan = 0;
//     // for (index=0;index<S2HAT_params->nstokes*S2HAT_params->size_alm;index++){
//     //     if (!(input_local_alm[index].re == input_local_alm[index].re)){
//     //         // local_alm_s2hat[index].re = 0;
//     //         number_of_nan++;
//     //         }
//     //     if (!(input_local_alm[index].im == input_local_alm[index].im)){
//     //         // local_alm_s2hat[index].im = 0;
//     //         number_of_nan++;
//     //         }
//     // }
//     // printf(" \n"); fflush(stdout);
//     // printf("----!!!! ***input_local_alm post inv cov*** size_alm %d ; Number of nans %d \n", S2HAT_params->size_alm, number_of_nan); fflush(stdout);
//     // number_of_nan = 0;
//     // for (index=0;index<S2HAT_params->nstokes*S2HAT_params->size_alm;index++){
//     //     if (!(out_local_alm[index].re == out_local_alm[index].re)){
//     //         // local_alm_s2hat[index].re = 0;
//     //         number_of_nan++;
//     //         }
//     //     if (!(out_local_alm[index].im == out_local_alm[index].im)){
//     //         // local_alm_s2hat[index].im = 0;
//     //         number_of_nan++;
//     //         }
//     // }
//     // printf(" \n"); fflush(stdout);
//     // printf("----!!!! ***out_local_alm post inv cov*** size_alm %d ; Number of nans %d \n", S2HAT_params->size_alm, number_of_nan); fflush(stdout);
//     // int ell_value, index_2;
//     // number_of_nan = 0;
//     // for (ell_value=0;ell_value<lmax;ell_value++){
        
//     //     // printf("\n #### ell= %d \n", ell_value);
//     //     for (index=0; index < nstokes; index++){
//     //         for (index_2=0; index_2<nstokes; index_2++){
//     //             if (!(inv_covariance_matrix[ell_value][index*nstokes+index_2] == inv_covariance_matrix[ell_value][index*nstokes+index_2])){
//     //                 // printf(" --- NAN HERE : ell %d index %d index_2 %d value  EE %f BB %.10f -- ", ell_value, index, index_2, covariance_matrix[ell_value][0], covariance_matrix[ell_value][3]);
//     //                 number_of_nan++; 
//     //             }
//     //         }
//     //     }
//     //     // printf(" --- number_of_nan inv cov : %d \n", number_of_nan);
//     // }
//     // printf(" --- number_of_nan of inv cov post inv cov : %d \n", number_of_nan);
//     // long int index;
//     // for(index=0; index<nstokes*(lmax)*m_value; index++){
//     //     local_alm[index].re = out_local_alm[index].re;
//     //     local_alm[index].im = out_local_alm[index].im;
//     // }
//     // free(out_local_alm);
//     }
//     return 0;
// }


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
