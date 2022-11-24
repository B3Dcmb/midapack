/** @file domain_generalization.h
    @brief <b> Declaration of the routines in the generalization of the PCG variable domain.</b>
    @author Magdy Morshed
    @date November 2022 */


typedef struct data_var{
    /* Local parameters of S2HAT, dependent on each processor */

    double *local_map_pix; // Local map in pixel domain, with MAPPRAISER convention
    s2hat_dcomplex *local_alm; // Local alm in harmonic domain, with S2HAT convention

    // Specifications of the computation
    int domain_PCG_computation; // Domain chosen for PCG computation : 0 pixel, 1 harmonic
    int bool_apply_filter; // 0 no filter applied, 1 Wiener filter applied

} data_var;


