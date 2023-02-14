/** @file domain_generalization.h
    @brief <b> Declaration of the routines in the generalization of the PCG variable domain.</b>
    @author Magdy Morshed
    @date November 2022 */

// #include "s2hat.h"
// For later : add if/if_not for including or not s2hat

typedef struct PCG_var{
    /* Local parameters of S2HAT, dependent on each processor */

    double *local_map_pix; // Local map in pixel domain, with MAPPRAISER convention
    s2hat_dcomplex *local_alm; // Local alm in harmonic domain, with S2HAT convention

    // Specifications of the computation
    int domain_PCG_computation; // Domain chosen for PCG computation : 0 pixel, 1 harmonic
    int bool_apply_filter; // 0 no filter applied, 1 Wiener filter applied
    int nstokes; // Number of Stokes parameters for the computation, 3 usually means T, Q, U


    // S2HAT structures necessary for harmonic domain
    S2HAT_parameters *S2HAT_parameters; 
    // contains : S2HAT_GLOBAL_parameters *Global_param_s2hat and S2HAT_LOCAL_parameters *Local_param_s2hat;
    // Note that all data_var will point to the same S2HAT_GLOBAL_parameters and S2HAT_LOCAL_parameters for each processor


    // Update flags to check if an update needs to be done, by retransforming pixel2alm or alm2pixel
    int does_map_pixel_need_update; // 0 no update needed, 1 need update from local_alm
    int does_local_alm_need_update; /// 0 no update needed, 1 need update from local_map_pix
    // Beware, those updates can be communication costly
    // Also note that in case we want to update from other domains of computation, the flags can take other values 2, 3, ...

} PCG_var;

/* Initalize PCG_var structure */
int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix, s2hat_dcomplex *local_alm, int domain_PCG_computation, int bool_apply_filter, int nstokes, S2HAT_parameters *S2HAT_params);

/* Transforms local map pixels into local alm in harmonic domain and vice-versa*/
int global_map_2_harmonic(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);
int global_harmonic_2_map(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);

/* Update PCG_var structure */
int update_PCG_var(PCG_var *PCG_variable, Mat *A);