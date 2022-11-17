 /** @file s2hat_tools.h
    @brief <b> Declaration of the backbone routines of the Wiener-filtering extension of MAPPRAISER.</b>
    @author Magdy Morshed
    @date November 2022 */
 
 /* Full documentation for S2HAT here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    
#include <mpi.h>
#include "s2hat.h"


#ifndef DBL_MAX
#define DBL_MAX            1.79769313486231470e+308
#endif

// typedef struct S2HAT_GLOBAL_parameters S2HAT_GLOBAL_parameters;
typedef struct S2HAT_GLOBAL_parameters{
    /* Global parameters of S2HAT, to give to all processors */
    s2hat_pixeltype pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel;
    s2hat_pixparameters pixpar;

    int nside;
    int nlmax;
    int nmmax;
} S2HAT_GLOBAL_parameters;

// typedef struct S2HAT_LOCAL_parameters S2HAT_LOCAL_parameters;
typedef struct S2HAT_LOCAL_parameters{
    /* Local parameters of S2HAT, dependent on each processor */

    int gangrank;
    int gangsize;
    int gangroot;
    MPI_Comm gangcomm;


    int nmvals;
    int first_ring;
    int last_ring;
    int map_size;
    int* mvals; // size given by nmvals


    // Tools to precompute Legendre functions, but not taken into account by s2hat if plms=0, which is the default behaviour we choose
    int plms;
    long int nplm;
} S2HAT_LOCAL_parameters;

/* Get global s2hat structures which must be distributed to all processors*/
int get_main_s2hat_global_parameters(int nside, char *maskfile_path, s2hat_pixeltype pixelization_scheme, s2hat_scandef scan_sky_structure_pixel, s2hat_pixparameters pixpar);

/* Create wrapper structure s2hat of local parameters of s2hat, which will differ for all processors */
int init_s2hat_global_parameters(char *maskfile_path, int nside, int lmax, S2HAT_GLOBAL_parameters *Global_param_s2hat);

/* Create wrapper structure of s2hat of local parameters of s2hat, which will differ for all processors */
int init_s2hat_local_parameters_struct(S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters *Local_param_s2hat, int *mvals, int gangrank, int gangsize, int gangroot, MPI_Comm gangcomm);


/* Use s2hat routines to broadcast s2hat global structures */
void mpi_broadcast_s2hat_global_struc(S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat, int gangroot);

/* Free covariance matrix */
void free_covariance_matrix(double ** covariance_matrix_3x3, int lmax);

/* Free wrapper structures of s2hat */
void free_s2hat_parameters_struct(S2HAT_GLOBAL_parameters *Global_param_s2hat, S2HAT_LOCAL_parameters *Local_param_s2hat);






/* Function to read file corresponding to the mask */
void read_fits_mask(int nside, double *mask, char *path_mask_file, int col);

/* Function to transform the mask into binary (composed of 0 and 1 on pixel sky)*/
void make_mask_binary(double* mask, int* mask_binary, int f_sky, int npix);

/* Obtain c_ell array from c_ell path */
void read_fits_cells(int lmax, int number_correl, double *c_ell_array, char *path_list_file, int col);

/* Transform alm coefficients local_alm into a pixel map local_map_pix */
int apply_alm2pix(s2hat_dcomplex *local_alm, double *local_map_pix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);

/* Transform local pixel map into local alm coefficients */
int apply_pix2alm(double *local_map_pix, s2hat_dcomplex *local_alm, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);

/* Apply inverse of covariance matrix to local_alm */
int apply_inv_covariance_matrix_to_alm(s2hat_dcomplex *local_alm, double **inv_covariance_matrix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);





/* Transform alm to c_ell coefficients */
int alm2cls(s2hat_dcomplex *local_alm, double *c_ell_array, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);

/* General function to inverse matrix using LAPACK */
int get_inverse_matrix(int order_matrix, double* matrix_to_be_inverted);

/* Read c_ell to generate covariance matrix which will be in the form : covariance_matrix_3x3[lmax][9] with 9 being [TT, TE, TB, ET, EE, EB, BT, BE, BB] (with TE=ET, TB=BT and BE=EB) */
int get_covariance_matrix_3x3(char *c_ell_path, int number_correl, double **covariance_matrix_3x3, S2HAT_GLOBAL_parameters Global_param_s2hat);

/* Function to obtain inverse of covariance matrix in harmonic domain, from given c_ells */
int get_inverse_covariance_matrix_3x3(char *c_ell_path, int number_correl, double **inverse_covariance_matrix, S2HAT_GLOBAL_parameters Global_param_s2hat);


