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

#include "s2hat.h"
#include "midapack.h"
#include "s2hat_tools.h"





int apply_alm2pix(s2hat_dcomplex *local_alm, double *local_map_pix, int nstokes, S2HAT_parameters *S2HAT_params){
    /* Transform alm coefficients local_alm into a pixel map local_map_pix, 
    all details here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Calm2map.html 

    local_alm is a 4-dimensional array in the form :
        (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
        (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
    in the form I, E, B
    */
    
    S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters Local_param_s2hat = *(S2HAT_params->Local_param_s2hat);
    // Getting global and local S2HAT parameters needed for the computation

    int nmaps = 1; // We only provide 1 input set of alm coefficient
    // int nstokes = 3; // We want all T, Q and U maps
    int lda = nstokes; // We choose the HEALPIX convention with local_alm in the form (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps)
    int spin;

    switch(nstokes)
    {
        case 1: // Case only intensity
        case 3: // Case with both intensity and polarization
            s2hat_alm2map(Local_param_s2hat.plms, Global_param_s2hat.pixelization_scheme, Global_param_s2hat.scan_sky_structure_pixel, Global_param_s2hat.nlmax, Global_param_s2hat.nmmax, 
                Local_param_s2hat.nmvals, Local_param_s2hat.mvals, nmaps, nstokes, 
                Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size, local_map_pix, lda, 
                local_alm, Local_param_s2hat.nplm, NULL, 
                Local_param_s2hat.gangsize, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
            // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0

        case 2: // Case with only polarization
            spin=2;
            s2hat_alm2map_spin( Global_param_s2hat.pixelization_scheme, Global_param_s2hat.scan_sky_structure_pixel, spin, Global_param_s2hat.nlmax, Global_param_s2hat.nmmax,
                Local_param_s2hat.nmvals, Local_param_s2hat.mvals, nmaps,
                Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat.gangsize, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
    }
    return 0;
}



int apply_pix2alm(double *local_map_pix, s2hat_dcomplex *local_alm, int nstokes, S2HAT_parameters *S2HAT_params){
    /* Transform pixel map into alm coefficients, all details here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Cmap2alm.html 
     local_alm is a 4-dimensional array in the form :
        (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
        (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
    Here the HEALpix convention has been chosen
    Output will be in the form I, E, B
    */
    S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);
    S2HAT_LOCAL_parameters Local_param_s2hat = *(S2HAT_params->Local_param_s2hat);
    // Getting global and local S2HAT parameters needed for the computation

    int nmaps = 1; // We only provide 1 input set of alm coefficient
    // int nstokes = 3; // We provide all 3 T, Q and U maps
    int lda = nstokes; // We choose the HEALPIX convention with local_alm in the form (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps)

    double *local_w8ring;
    int i_ring;

    int spin;

    // printf("Test pix2alm entry !\n");
    // fflush(stdout);
    // local_w8ring=(double *)calloc( nstokes*(Local_param_s2hat.last_ring-Local_param_s2hat.first_ring+1), sizeof(double));
    local_w8ring = (double *) malloc( nstokes*(Local_param_s2hat.last_ring-Local_param_s2hat.first_ring+1)*sizeof(double));
    for( i_ring=0; i_ring< nstokes*(Local_param_s2hat.last_ring-Local_param_s2hat.first_ring+1); i_ring++) {
            local_w8ring[i_ring]=1.;
        }



    
    switch(nstokes)
    {
        case 1: // Case only intensity
        case 3: // Case with both intensity and polarization
            s2hat_map2alm(Local_param_s2hat.plms, Global_param_s2hat.pixelization_scheme, Global_param_s2hat.scan_sky_structure_pixel, Global_param_s2hat.nlmax, Global_param_s2hat.nmmax, 
                Local_param_s2hat.nmvals, Local_param_s2hat.mvals, nmaps, nstokes, 
                Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, local_w8ring, Local_param_s2hat.map_size, local_map_pix, lda, local_alm, 
                Local_param_s2hat.nplm, NULL,
                Local_param_s2hat.gangsize, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
            // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0

        case 2: // Case with only polarization
            spin=2;            
            s2hat_map2alm_spin( Global_param_s2hat.pixelization_scheme, Global_param_s2hat.scan_sky_structure_pixel, spin, Global_param_s2hat.nlmax, Global_param_s2hat.nmmax,
                Local_param_s2hat.nmvals, Local_param_s2hat.mvals, nmaps, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring,
		        local_w8ring, Local_param_s2hat.map_size, local_map_pix, lda, local_alm,
                Local_param_s2hat.gangsize, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
            // Computing alm2map in the specific case where only Q,U are provided
    }

    free(local_w8ring);
    return 0;
}

// int gather_map(double *local_map_pix, double *full_sky_map, int nstokes, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){

//     // int nstokes = 3;
//     collect_map(Global_param_s2hat.pixelization_scheme, 1, 0, nstokes, full_sky_map, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size,
//             local_map_pix, Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
// }


int apply_inv_covariance_matrix_to_alm(s2hat_dcomplex *input_local_alm, s2hat_dcomplex *out_local_alm, double **inv_covariance_matrix, int nstokes, S2HAT_parameters *S2HAT_params){
    /* Apply inverse of covariance matrix to input_local_alm */

    // S2HAT_GLOBAL_parameters Global_param_s2hat = *(S2HAT_params->Global_param_s2hat);
    // S2HAT_LOCAL_parameters Local_param_s2hat = *(S2HAT_params->Local_param_s2hat);


    int ell_value, m_value, index_stokes, line_index;
    // int nstokes = 3;
    int lmax = S2HAT_params->Global_param_s2hat->nlmax;

    int nmvals = S2HAT_params->Local_param_s2hat->nmvals; // Total number of m values
    // int *mvals = Local_param_s2hat->mvals; // Values of m the considered processor contain

    double res_real, res_imag;


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
    // long int index;
    // for(index=0; index<nstokes*(lmax+1)*m_value; index++){
    //     local_alm[index].re = out_local_alm[index].re;
    //     local_alm[index].im = out_local_alm[index].im;
    // }
    // free(out_local_alm);
}
