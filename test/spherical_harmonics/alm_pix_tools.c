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





int apply_alm2pix(s2hat_dcomplex *local_alm, double *local_map_pix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
    /* Transform alm coefficients local_alm into a pixel map local_map_pix, 
    all details here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Calm2map.html 

    local_alm is a 4-dimensional array in the form :
        (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
        (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
    in the form I, E, B
    */

    int nmaps = 1; // We only provide 1 input set of alm coefficient
    int nstokes = 3; // We want all T, Q and U maps
    int lda = nstokes; // We choose the HEALPIX convention with local_alm in the form (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps)

    // printf("Test 4 !\n");
    // fflush(stdout);
    s2hat_alm2map(Local_param_s2hat.plms, Global_param_s2hat.pixelization_scheme, Global_param_s2hat.scan_sky_structure_pixel, Global_param_s2hat.nlmax, Global_param_s2hat.nmmax, 
        Local_param_s2hat.nmvals, Local_param_s2hat.mvals, nmaps, nstokes, 
	    Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size, local_map_pix, lda, 
        local_alm, Local_param_s2hat.nplm, NULL, 
		Local_param_s2hat.gangsize, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
    // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0
    // printf("Test 5 !\n");
    // fflush(stdout);
    return 0;
}

// int apply_alm2pix_v2(s2hat_dcomplex *local_alm, double *local_map_pix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
//     /* Transform alm coefficients local_alm into a pixel map local_map_pix, 
//     all details here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Calm2map.html 

//     local_alm is a 4-dimensional array in the form :
//         (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
//         (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
//     in the form I, E, B
//     */

//     int nmaps = 1; // We only provide 1 input set of alm coefficient
//     int nstokes = 3; // We want all T, Q and U maps
//     int lda = nstokes; // We choose the HEALPIX convention with local_alm in the form (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps)

//     // printf("Test 4 !\n");
//     // fflush(stdout);
//     double *local_map_pix_T;
//     int index=0;
//     for(index=0;index<)
//     local_map_pix_T
//     s2hat_alm2map(Local_param_s2hat.plms, Global_param_s2hat.pixelization_scheme, Global_param_s2hat.scan_sky_structure_pixel, Global_param_s2hat.nlmax, Global_param_s2hat.nmmax, 
//         Local_param_s2hat.nmvals, Local_param_s2hat.mvals, nmaps, nstokes, 
// 	    Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size, local_map_pix, lda, 
//         local_alm, Local_param_s2hat.nplm, NULL, 
// 		Local_param_s2hat.gangsize, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
//     // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0
//     // printf("Test 5 !\n");
//     // fflush(stdout);
//     return 0;
// }




int apply_pix2alm(double *local_map_pix, s2hat_dcomplex *local_alm, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
    /* Transform pixel map into alm coefficients, all details here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/Cmap2alm.html 
     local_alm is a 4-dimensional array in the form :
        (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps), if lda == nstokes;      (HEALpix convention)
        (0:nlmax,0:nmvals-1,1:nstokes,1:nmaps), if lda == nlmax;      (S2HAT convention)
    Here the HEALpix convention has been chosen
    Output will be in the form I, E, B
    */

    int nmaps = 1; // We only provide 1 input set of alm coefficient
    int nstokes = 3; // We provide all 3 T, Q and U maps
    int lda = nstokes; // We choose the HEALPIX convention with local_alm in the form (1:nstokes,0:nlmax,0:nmvals-1,1:nmaps)

    double *local_w8ring;
    int i_ring;

    // printf("Test pix2alm entry !\n");
    // fflush(stdout);
    // local_w8ring=(double *)calloc( nstokes*(Local_param_s2hat.last_ring-Local_param_s2hat.first_ring+1), sizeof(double));
    local_w8ring = (double *) malloc( nstokes*(Local_param_s2hat.last_ring-Local_param_s2hat.first_ring+1)*sizeof(double));
    for( i_ring=0; i_ring< nstokes*(Local_param_s2hat.last_ring-Local_param_s2hat.first_ring+1); i_ring++) {
            local_w8ring[i_ring]=1.;
        }
    
    // printf("%d ## Global parameters ! --- lmax %d, mmax %d, nside %d \n", Local_param_s2hat.gangrank, Global_param_s2hat.nlmax, Global_param_s2hat.nmmax, Global_param_s2hat.nside);
    // s2hat_pixeltype pixelization_scheme_2 = Global_param_s2hat.pixelization_scheme;
    // s2hat_scandef scan_sky_structure_pixel_2 = Global_param_s2hat.scan_sky_structure_pixel;
    // s2hat_pixparameters pixpar_2 = Global_param_s2hat.pixpar;
    // printf("###### Test init_s2hat_global_parameters - %d \n", Local_param_s2hat.gangrank);
    // printf(" %d - pixelization_scheme :  type %d, npixsall %d, nphmx %d, nringsall %d \n", Local_param_s2hat.gangrank, pixelization_scheme_2.type, pixelization_scheme_2.npixsall, pixelization_scheme_2.nphmx, pixelization_scheme_2.nringsall);
    // printf(" %d - scan_sky_structure_pixel :  npixsobs %d, nringsobs %d \n", Local_param_s2hat.gangrank, scan_sky_structure_pixel_2.npixsobs, scan_sky_structure_pixel_2.nringsobs);
    // printf(" %d - pixpar :  par1 %d   par2 %d \n", Local_param_s2hat.gangrank, pixpar_2.par1, pixpar_2.par2);

    // printf("%d ## Test map2alm s2hat entry ! --- plms %d, nmvals %d, first_ring %d, last_ring %d, map_size %d, nplm %ld, mvals[0] %d \n", Local_param_s2hat.gangrank,
    //         Local_param_s2hat.plms, Local_param_s2hat.nmvals, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size, Local_param_s2hat.nplm,
    //         Local_param_s2hat.mvals[0]);
    // fflush(stdout);


    s2hat_map2alm(Local_param_s2hat.plms, Global_param_s2hat.pixelization_scheme, Global_param_s2hat.scan_sky_structure_pixel, Global_param_s2hat.nlmax, Global_param_s2hat.nmmax, 
            Local_param_s2hat.nmvals, Local_param_s2hat.mvals, nmaps, nstokes, 
            Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, local_w8ring, Local_param_s2hat.map_size, local_map_pix, lda, local_alm, 
            Local_param_s2hat.nplm, NULL,
            Local_param_s2hat.gangsize, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
        // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0
    // printf("Test pix2alm out ! - %d \n", Local_param_s2hat.gangrank);
    // fflush(stdout);
    free(local_w8ring);

    return 0;
}

int gather_map(double *local_map_pix, double *full_sky_map, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){

    int nstokes = 3;
    collect_map(Global_param_s2hat.pixelization_scheme, 1, 0, nstokes, full_sky_map, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size,
            local_map_pix, Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);

}


int apply_inv_covariance_matrix_to_alm(s2hat_dcomplex *local_alm, double **inv_covariance_matrix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
    /* Apply inverse of covariance matrix to local_alm */

    // distribute_alms to apply to covariance matrix ?
    
    

    int ell_value, m_value, index_stokes, line_index;
    int nstokes = 3;
    int lmax = Global_param_s2hat.nlmax;

    int nmvals = Local_param_s2hat.nmvals; // Total number of m values
    // int *mvals = Local_param_s2hat->mvals; // Values of m the considered processor contain

    double res_real, res_imag;

    s2hat_dcomplex *new_local_alm;
    new_local_alm = (s2hat_dcomplex *) calloc(nstokes*(lmax+1)*m_value,sizeof(s2hat_dcomplex));


    for(ell_value=0; ell_value < lmax+1; ell_value++){
        for(m_value=0; m_value < nmvals; m_value++){
            for (index_stokes=0; index_stokes<nstokes; index_stokes++){
                res_real = 0;
                res_imag = 0;
                for (line_index=0; line_index < nstokes; line_index++){
                    // res_real += local_alm[line_index*(lmax+1)*nmvals + ell_value*nmvals + m_value].re * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                    // res_imag += local_alm[line_index*(lmax+1)*nmvals + ell_value*nmvals + m_value].im * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                    
                    // res_real += local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + line_index].re * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];
                    // res_imag += local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + line_index].im * inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index];

                    res_real += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + line_index].re;
                    res_imag += inv_covariance_matrix[ell_value][index_stokes*nstokes + line_index] * local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + line_index].im;
                    // Maybe with [line_index*nstokes + index_stokes] ?
                }
                // local_alm[index_stokes*(lmax+1)*nmvals + ell_value*nmvals + m_value].re = res_real;
                // local_alm[index_stokes*(lmax+1)*nmvals + ell_value*nmvals + m_value].im = res_imag;
                new_local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + index_stokes].re = res_real;
                new_local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + index_stokes].im = res_imag;
            }
            // Verify it is not applied to part where a_lm not defined !!!
        }
    }
    long int index;
    for(index=0; index<nstokes*(lmax+1)*m_value; index++){
        local_alm[index].re = new_local_alm[index].re;
        local_alm[index].im = new_local_alm[index].im;
    }
    free(new_local_alm);
}
