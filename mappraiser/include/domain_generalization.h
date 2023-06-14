/** @file domain_generalization.h
    @brief <b> Declaration of the routines in the generalization of the PCG variable domain.</b>
    @author Magdy Morshed
    @date November 2022 */

// #include "s2hat.h"
// For later : add if/if_not for including or not s2hat
// #include "spherical_harmonics/s2hat_tools.h"
#include "midapack.h"

typedef struct Harmonic_superstruct Harmonic_superstruct;

typedef struct PCG_var{
    /* Local parameters of S2HAT, dependent on each processor */

    double *local_map_pix; // Local map in pixel domain, with MAPPRAISER convention
    // s2hat_dcomplex *local_alm; // Local alm in harmonic domain, with S2HAT convention

    // Specifications of the computation
    int domain_PCG_computation; // Domain chosen for PCG computation : 0 pixel, 1 harmonic
    int bool_apply_filter; // 0 no filter applied, 1 Wiener filter applied
    // int nstokes; // Number of Stokes parameters for the computation, where we assume : 1 for intensity only, 2 for polarization only, 3 for everything (T, Q, U)

    // Update flags to check if an update needs to be done, by retransforming pixel2alm or alm2pixel
    // int does_map_pixel_need_update; // 0 no update needed, 1 need update from local_alm
    // int does_local_alm_need_update; /// 0 no update needed, 1 need update from local_map_pix
    // Beware, those updates can be communication costly
    // Also note that in case we want to update from other domains of computation, the flags can take other values 2, 3, ...
} PCG_var;

struct Harmonic_superstruct{
    /* Harmonic superstructure for harmonic operations, defined here using S2HAT */

    // S2HAT structures necessary for harmonic domain
    S2HAT_parameters S2HAT_params;
    // contains : iles_path_WIENER_FILTER Files_WF_struct, S2HAT_GLOBAL_parameters Global_param_s2hat and S2HAT_LOCAL_parameters Local_param_s2hat;

    Butterfly_superstruct MAPPRAISER_to_S2HAT_butterfly;
    Butterfly_superstruct S2HAT_to_MAPPRAISER_butterfly;
    // Structures to use butterfly scheme, to transition from the MAPPRAISER distribution to the S2HAT one, and vice-versa
};


/* Initalize PCG_var structure */
int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix);

/* Initialize harmonic superstructure */
int init_harmonic_superstruct(Mat *A, Harmonic_superstruct *Harm_struct, int *mask_binary, int nside, int lmax, char *c_ell_path, int number_correlations);

/* Transforms local map pixels into local alm in harmonic domain and vice-versa*/
int global_map_2_harmonic(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, Harmonic_superstruct *Harmonic_sup);
int global_harmonic_2_map(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A,  Harmonic_superstruct *Harmonic_sup);

int get_mask_from_indices(Mat *A, int *mask_binary, int nside, int root);

int free_harmonic_superstruct(Harmonic_superstruct *Harmonic_sup, int rank);
// int free_PCG_var(PCG_var *PCG_var_obj);