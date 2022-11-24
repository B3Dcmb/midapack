
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




int get_main_s2hat_global_parameters(int nside, char *maskfile_path, s2hat_pixeltype *pixelization_scheme, s2hat_scandef *scan_sky_structure_pixel, s2hat_pixparameters *pixpar, bool use_mask_file){
    /* Get s2hat structure which must be distributed by all processors
    All processors must have those same s2hat structure 
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    
    double *mask;
    int *mask_binary;
    int f_sky=0;
    long npix = nside*nside*12;

    int pixchoice = PIXCHOICE_HEALPIX;   /* Choice  of convention to use for the pixelization scheme, here HEALPIX is chosen, 
    see https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/CparsAndStruct.html for more details*/

    pixpar->par1 = nside; /* NSIDE if HEALPIX convention is used for the transforms, 
    see https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/CsetPixelization.html */

    // printf("Pixchoice : %d \n", pixpar->par1);

    // printf( " //// Define s2hat structures\n"); 

    if (use_mask_file){
        // printf( " //// Reading mask\n");
        mask = (double *) malloc( npix*sizeof(double));
        read_fits_mask( nside, mask, maskfile_path, 1); /* Retrieve mask */

        mask_binary = (int *) malloc( npix*sizeof(int));
        make_mask_binary(mask, mask_binary, &f_sky, npix); /* Make mask binary (only with 0 and 1) */

        free(mask);
    }
    else{
        // printf( " //// No use of mask - fill fake mask with 1 \n");
        mask_binary = (int *) malloc( npix*sizeof(int));
        int index;
        for (index=0;index<npix;index++){
            mask_binary[index] = 1;
        }
        f_sky = npix;
    }
    // printf( " //// Setting s2hat structure \n"); fflush(stdout);

    set_pixelization(pixchoice, *pixpar, pixelization_scheme);   /* Set C pixelization structure */
    mask2scan(mask_binary, *pixelization_scheme, scan_sky_structure_pixel); /* Set scan pixelization : s2hat structure containing all the info about sky coverage needed by the transforms */
    /* Scan pixelization will be used to avoid doing unecessary calculations */
        
    // printf( "F_sky = %f%%\n", (double)f_sky/(double)npix*100.); fflush(stdout);
    return 0;
}


int init_s2hat_global_parameters(Files_path_WIENER_FILTER Files_WF_struct, int nside, int lmax, S2HAT_GLOBAL_parameters *Global_param_s2hat){
    /* Create s2hat structure of global parameters of s2hat, which must be distributed to all processors
    All processors must have those same s2hat structure 
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    
    s2hat_pixeltype pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel;
    s2hat_pixparameters pixpar;
    char *maskfile_path = Files_WF_struct.maskfile_path;
    bool use_mask_file = Files_WF_struct.use_mask_file;
    
    get_main_s2hat_global_parameters(nside, maskfile_path, &pixelization_scheme, &scan_sky_structure_pixel, &pixpar, use_mask_file);

    Global_param_s2hat->pixelization_scheme = pixelization_scheme;
    Global_param_s2hat->scan_sky_structure_pixel = scan_sky_structure_pixel;
    Global_param_s2hat->pixpar = pixpar;
    
    Global_param_s2hat->nside = nside;
    Global_param_s2hat->nlmax = lmax-1; // S2HAT will generate alms between 0 and nlmax included, so we have to give lamx decreased by 1 to get the lmax requested
    Global_param_s2hat->nmmax = lmax-1;
}

int init_MPI_struct_s2hat_local_parameters(S2HAT_LOCAL_parameters *Local_param_s2hat, int gangrank, int gangsize, int gangroot, MPI_Comm gangcomm){
    /* Initialize MPI parameters of s2hat structure of local parameters, which will differ for all processors
    Those local parameters are used to improve the computation efficiency
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 

    Local_param_s2hat->gangrank = gangrank;
    Local_param_s2hat->gangsize = gangsize;
    Local_param_s2hat->gangroot = gangroot;
    Local_param_s2hat->gangcomm = gangcomm;

    return 0;
}


int init_s2hat_local_parameters_struct(S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters *Local_param_s2hat){
    /* Create s2hat structure of local parameters of s2hat, which will differ for all processors
    
    !!! BEWARE : MPI structure is assumed to already have been assigned using init_MPI_struct_s2hat_local_parameters !!!

    Those local parameters are used to improve the computation efficiency
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    s2hat_pixeltype pixelization_scheme = Global_param_s2hat.pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel = Global_param_s2hat.scan_sky_structure_pixel;
    int nlmax = Global_param_s2hat.nlmax;
    int nmmax = Global_param_s2hat.nmmax;


    int plms=0, nmvals, first_ring, last_ring, map_size;

    long int nplm;
    int *mvals;

    // printf("%d ----- Test get_local_data_sizes : plms %d, nlmax %d, nmmax %d, gangrank %d, gangsize %d, gangroot %d \n",  Local_param_s2hat->gangrank, plms, nlmax, nmmax, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot);
    // fflush(stdout);
    
    get_local_data_sizes( plms, pixelization_scheme, scan_sky_structure_pixel, nlmax, nmmax, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize,
			&nmvals, &first_ring, &last_ring, &map_size, &nplm, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);
    
    // int nmvals_2;
    // nmvals_2 = nummvalues(Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, nmmax);
    // printf("nmvals - %d, %d for rank %d \n", nmvals, nmvals_2, Local_param_s2hat->gangrank );
    // fflush(stdout);
    /* Estimates size of local data object, and set first_ring, last_ring, map_size
    see more info on https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/CgetLocalDataSizes.html */

    // printf("--------Local_sizes obtained ! - %d : %d %d \n", Local_param_s2hat->gangrank, first_ring, last_ring);
    // fflush(stdout);

    mvals = (int *) calloc( nmvals, sizeof( int));    /// TO FREE LATER !!!!!!
    find_mvalues( Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, nmmax, nmvals, mvals);

    // printf("--------Mvalues obtained ! - %d : %d \n", Local_param_s2hat->gangrank, mvals[0]);
    // fflush(stdout);

    Local_param_s2hat->plms = plms;

    Local_param_s2hat->nmvals = nmvals;
    Local_param_s2hat->first_ring = first_ring;
    Local_param_s2hat->last_ring = last_ring;
    Local_param_s2hat->map_size = map_size;
    Local_param_s2hat->nplm = nplm;

    Local_param_s2hat->mvals = mvals;
    return 0;
}



void mpi_broadcast_s2hat_global_struc(S2HAT_GLOBAL_parameters *Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
    /* Use s2hat routines to broadcast s2hat structures */
    MPI_pixelizationBcast( &(Global_param_s2hat->pixelization_scheme), Local_param_s2hat.gangroot, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
    MPI_scanBcast(Global_param_s2hat->pixelization_scheme, &(Global_param_s2hat->scan_sky_structure_pixel), Local_param_s2hat.gangroot, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
    MPI_Bcast( &(Global_param_s2hat->pixpar.par1), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    MPI_Bcast( &(Global_param_s2hat->pixpar.par2), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    MPI_Bcast( &(Global_param_s2hat->nlmax), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    MPI_Bcast( &(Global_param_s2hat->nmmax), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    MPI_Bcast( &(Global_param_s2hat->nside), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
}


void free_covariance_matrix(double ** covariance_matrix_3x3, int lmax){
    int ell_value;
    for (ell_value=0; ell_value<lmax+1; ell_value++){
        free(covariance_matrix_3x3[ell_value]);
    }
    free(covariance_matrix_3x3);    
}

void free_s2hat_GLOBAL_parameters_struct(S2HAT_GLOBAL_parameters *Global_param_s2hat){
    destroy_pixelization(Global_param_s2hat->pixelization_scheme);
    destroy_scan(Global_param_s2hat->scan_sky_structure_pixel);
    free(Global_param_s2hat);
}

void free_s2hat_LOCAL_parameters_struct(S2HAT_LOCAL_parameters *Local_param_s2hat){
    free(Local_param_s2hat->mvals);
    free(Local_param_s2hat);
}
