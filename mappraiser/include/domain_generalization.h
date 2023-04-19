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
    int nstokes; // Number of Stokes parameters for the computation, where we assume : 1 for intensity only, 2 for polarization only, 3 for everything (T, Q, U)

    // Update flags to check if an update needs to be done, by retransforming pixel2alm or alm2pixel
    // int does_map_pixel_need_update; // 0 no update needed, 1 need update from local_alm
    // int does_local_alm_need_update; /// 0 no update needed, 1 need update from local_map_pix
    // Beware, those updates can be communication costly
    // Also note that in case we want to update from other domains of computation, the flags can take other values 2, 3, ...
} PCG_var;


typedef struct Butterfly_struct{

    int		*com_indices, com_count; // communicated indices, and size
    int		steps;			         // number of steps in the butterfly scheme
    int		*nS, *nR;		         // number of indices (to send and to receive); size = steps
    int		**R, **S;		         // sending or receiving indices

    int classic_or_reshuffle_butterfly; 
    // Flag to indicate if classic or reshuffled butterfly is used, e.g. if we expect the pixel distributions prior and after communication to be the same or different
    // 0 for classic butterfly scheme with same pixel distributions prior/after ; 1 for after

    // int do_we_need_to_project_into_different_scheme;
    // Flag to indicate if we need to project the input pixel distribution into a different scheme
    
    int *projector_values;
    // Projector for local values of map : allows to project the values from one pixel distribution to another
    
    int *ordered_indices;
    // Stored ordered indices, only used for ring2nest and nest2ring transitions in which case it corresponds to ordered ring indices
} Butterfly_struct;

struct Harmonic_superstruct{
    /* Harmonic superstructure for harmonic operations, defined here using S2HAT */

    // S2HAT structures necessary for harmonic domain
    S2HAT_parameters *S2HAT_parameters;
    // contains : S2HAT_GLOBAL_parameters *Global_param_s2hat and S2HAT_LOCAL_parameters *Local_param_s2hat;
    // Note that all data_var will point to the same S2HAT_GLOBAL_parameters and S2HAT_LOCAL_parameters for each processor

    Butterfly_struct *S2HAT_to_MAPPRAISER;
    Butterfly_struct *MAPPRAISER_to_S2HAT;
    // Structures to use butterfly scheme, here to store the necessary structures to transition ring2nest and nest2ring
};


/* Initalize PCG_var structure */
int initialize_PCG_var_struct(PCG_var *PCG_variable, double *local_map_pix, int domain_PCG_computation, int bool_apply_filter, int nstokes);

/* Initialize harmonic superstructure */
int init_harmonic_superstruct(int is_pixel_scheme_MAPPRAISER_ring, Mat *A, Harmonic_superstruct *Harm_struct, int *mask_binary);

/* Transforms local map pixels into local alm in harmonic domain and vice-versa*/
int global_map_2_harmonic(double* local_pixel_map_MAPPRAISER, s2hat_dcomplex *local_alm_s2hat, Mat *A, Harmonic_superstruct *Harmonic_sup);
int global_harmonic_2_map(s2hat_dcomplex *local_alm_s2hat, double* local_pixel_map_MAPPRAISER, Mat *A,  Harmonic_superstruct *Harmonic_sup);

/* Update PCG_var structure */
// int update_PCG_var(PCG_var *PCG_variable, Mat *A);

int get_mask_from_indices(Mat *A, int *mask_binary, int nside, int root);