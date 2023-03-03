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




int get_main_s2hat_global_parameters(int nside, int *mask_binary, s2hat_pixeltype *pixelization_scheme, s2hat_scandef *scan_sky_structure_pixel, s2hat_pixparameters *pixpar){
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

    // if (use_mask_file){
    //     // printf( " //// Reading mask\n");
    //     mask = (double *) malloc( npix*sizeof(double));
    //     read_fits_mask( nside, mask, maskfile_path, 1); /* Retrieve mask */

    //     mask_binary = (int *) malloc( npix*sizeof(int));
    //     make_mask_binary(mask, mask_binary, &f_sky, npix); /* Make mask binary (only with 0 and 1) */

    //     free(mask);
    // }
    // else{
    //     // printf( " //// No use of mask - fill fake mask with 1 \n");
    //     mask_binary = (int *) malloc( npix*sizeof(int));
    //     int index;
    //     for (index=0;index<npix;index++){
    //         mask_binary[index] = 1;
    //     }
    //     f_sky = npix;
    // }
    if (mask_binary == NULL){
        // printf( " //// Reading mask\n");
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
    // free(mask_binary); // ????? Problem ?
    return 0;
}


int init_s2hat_global_parameters(Files_path_WIENER_FILTER Files_WF_struct, int *mask_binary, int lmax, S2HAT_GLOBAL_parameters *Global_param_s2hat){
    /* Create s2hat structure of global parameters of s2hat, which must be distributed to all processors
    All processors must have those same s2hat structure 
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    
    s2hat_pixeltype pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel;
    s2hat_pixparameters pixpar;
    // char *maskfile_path = Files_WF_struct.maskfile_path;
    // bool use_mask_file = Files_WF_struct.use_mask_file;
    int nside = Files_WF_struct.nside;
    
    get_main_s2hat_global_parameters(nside, mask_binary, &pixelization_scheme, &scan_sky_structure_pixel, &pixpar);

    Global_param_s2hat->pixelization_scheme = pixelization_scheme;
    Global_param_s2hat->scan_sky_structure_pixel = scan_sky_structure_pixel;
    Global_param_s2hat->pixpar = pixpar;
    
    Global_param_s2hat->nside = nside;
    Global_param_s2hat->nlmax = lmax-1; // S2HAT will generate alms between 0 and nlmax included, so we have to give lamx decreased by 1 to get the lmax requested
    Global_param_s2hat->nmmax = lmax-1;
}


int init_MPI_struct_s2hat_local_parameters(S2HAT_LOCAL_parameters *Local_param_s2hat, int number_ranks_s2hat, MPI_Comm initcomm){
    /* Initialize MPI parameters of s2hat structure of local parameters, which will differ for all processors
    Create a MPI group with the processors which will be used for S2HAT operations
    
    Those local parameters are used to improve the computation efficiency
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */

    int s2hat_rank, s2hat_size, s2hat_root;

    int initrank, initsize;
    int initroot = 0;

    MPI_Comm_rank(initcomm, &initrank);
    MPI_Comm_size(initcomm, &initsize);

    MPI_Comm s2hat_comm;

   
    // Test case if we don't have to divide into two substructures of communicators
    if (number_ranks_s2hat >= initsize){
        mpi_create_subset(number_ranks_s2hat, initcomm, &s2hat_comm);

        MPI_Comm_rank ( s2hat_comm, &s2hat_rank ); 
        MPI_Comm_size ( s2hat_comm, &s2hat_size );
        if (s2hat_comm == MPI_COMM_NULL){
            s2hat_rank = -1;
            s2hat_size = 0;
            initroot = 0;
        }
    }
    else{
        s2hat_rank = initrank;
        s2hat_size = initsize;
        s2hat_comm = initcomm;
    }

    Local_param_s2hat->gangrank = s2hat_rank;
    Local_param_s2hat->gangsize = s2hat_size;
    Local_param_s2hat->gangroot = initroot;
    Local_param_s2hat->gangcomm = s2hat_comm;

    return 0;
}


int init_s2hat_local_parameters_struct(S2HAT_GLOBAL_parameters Global_param_s2hat, int nstokes, S2HAT_LOCAL_parameters *Local_param_s2hat){
    /* Create s2hat structure of local parameters of s2hat, which will differ for all processors
    
    !!! BEWARE : MPI structure is assumed to already have been assigned using init_MPI_struct_s2hat_local_parameters !!!

    Those local parameters are used to improve the computation efficiency
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    long int *pixel_numbered_ring;
    s2hat_pixeltype pixelization_scheme = Global_param_s2hat.pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel = Global_param_s2hat.scan_sky_structure_pixel;
    int nlmax = Global_param_s2hat.nlmax;
    int nmmax = Global_param_s2hat.nmmax;
    long int number_pixel_total = 12*((Global_param_s2hat.nside)**2);

    int plms=0, nmvals, first_ring, last_ring, map_size;

    long int nplm;
    int *mvals;

    // printf("%d ----- Test get_local_data_sizes : plms %d, nlmax %d, nmmax %d, gangrank %d, gangsize %d, gangroot %d \n",  Local_param_s2hat->gangrank, plms, nlmax, nmmax, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot);
    // fflush(stdout);
    
    if (Local_param_s2hat->gangrank != -1){
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

        Local_param_s2hat->first_pixel_number = Global_param_s2hat.pixelization_scheme.fpix[first_ring];
        Local_param_s2hat->last_pixel_number = Global_param_s2hat.pixelization_scheme.fpix[last_ring]; // -1 ???

        first_pixel_number_ring = Global_param_s2hat.pixelization_scheme.fpix[first_ring];
        // Getting first pixel number of ring probed by local proc, in S2HAT convention
        last_pixel_number_ring = Global_param_s2hat.pixelization_scheme.fpix[last_ring]; 
        // Getting last pixel number of ring probed by local proc, in S2HAT convention
        number_pixels_local = last_pixel_number - first_pixel_number;
        
        
        pixel_numbered_ring = (long int *)malloc(number_pixels_local*nstokes * sizeof(long int));
        for(i=0; i<number_pixels_local; i++)
        {
            for (j=0; j<nstokes; j++)
                pixel_numbered_ring[i + j*number_pixels_local] = i + first_pixel_number + j*number_pixel_total;
        }
        
        Local_param_s2hat->pixel_numbered_ring = pixel_numbered_ring;

        Local_param_s2hat->nmvals = nmvals;
        Local_param_s2hat->first_ring = first_ring;
        Local_param_s2hat->last_ring = last_ring;
        Local_param_s2hat->map_size = map_size;
        
        Local_param_s2hat->plms = plms; // Structures for storing Legendre polynomials, by default put to 0 and not used
        Local_param_s2hat->nplm = nplm; // Structures for storing Legendre polynomials, by default put to 0 and not used

        Local_param_s2hat->mvals = mvals;
    }
    else{
        Local_param_s2hat->map_size = 0;
        Local_param_s2hat->pixel_numbered_ring = NULL;
    }
    return 0;
}


int init_s2hat_parameters_superstruct(Files_path_WIENER_FILTER *Files_WF_struct, int *mask_binary, int nstokes, S2HAT_parameters *S2HAT_params, MPI_Comm world_comm)
{   // Initalize both S2HAT_GLOBAL_parameters and S2HAT_LOCAL_parameters for superstructure of S2HAT
    
    S2HAT_GLOBAL_parameters *Global_param_s2hat = (S2HAT_GLOBAL_parameters *) malloc( 1 * sizeof(S2HAT_GLOBAL_parameters));
    init_s2hat_global_parameters(*Files_WF_struct, mask_binary, Files_WF_struct->lmax_Wiener_Filter, Global_param_s2hat); 
    // Initialization of Global_param_s2hat structure, for sky pixelization scheme and lmax_WF choice


    S2HAT_LOCAL_parameters *Local_param_s2hat = (S2HAT_LOCAL_parameters *) malloc( 1 * sizeof(S2HAT_LOCAL_parameters));
    // init_MPI_struct_s2hat_local_parameters(Local_param_s2hat, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm); 
    init_MPI_struct_s2hat_local_parameters(Local_param_s2hat, Global_param_s2hat->scan_sky_structure_pixel.nringsobs, world_comm); 
    
    init_s2hat_local_parameters_struct(*Global_param_s2hat, nstokes, Local_param_s2hat);
    // Initialization of Local_param_s2hat structure, including MPI parameters, first/last rings studied, size of pixels cut sky per rank, etc. -- see Wiener filter extension directory for more details

    S2HAT_params->Global_param_s2hat = Global_param_s2hat;
    S2HAT_params->Local_param_s2hat = Local_param_s2hat;
    // Initialization of final superstructure S2HAT_params

    S2HAT_params->size_alm = (Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax;
    S2HAT_params->nstokes = nstokes;
}


void mpi_broadcast_s2hat_global_struc(S2HAT_GLOBAL_parameters *Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat){
    /* Use s2hat routines to broadcast s2hat structures */
    if (Local_param_s2hat.gangrank != -1){
        MPI_pixelizationBcast( &(Global_param_s2hat->pixelization_scheme), Local_param_s2hat.gangroot, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
        MPI_scanBcast(Global_param_s2hat->pixelization_scheme, &(Global_param_s2hat->scan_sky_structure_pixel), Local_param_s2hat.gangroot, Local_param_s2hat.gangrank, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->pixpar.par1), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->pixpar.par2), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->nlmax), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->nmmax), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
        MPI_Bcast( &(Global_param_s2hat->nside), 1, MPI_INT, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    }
}

int distribute_full_sky_map_into_local_maps_S2HAT(double* full_sky_map, double *local_map_s2hat, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat, int nstokes){
    /* Distribute full sky map in ring ordering, with convention [npix, nstokes] in column-wise order among procs, into local maps */
    distribute_map(Global_param_s2hat.pixelization_scheme, 1, 0, nstokes, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, Local_param_s2hat.map_size, 
        local_map_s2hat, full_sky_map, Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    // 1 for the number of maps, 0 for the index of the current map

    return 0;
}


int collect_partial_map_from_pixels(double* local_map_s2hat, double *output_submap, int first_pix, int last_pix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat, int nstokes){
    // Collect specific pixels from all local_map_s2hat to form a submap given first and last pixels 
    int submap_size = last_pix - first_pix; // Submapsize given by pixel numbers

    collect_partialmap(Global_param_s2hat.pixelization_scheme, 1, 0, nstokes, first_pix, last_pix, 
        output_submap, Local_param_s2hat.first_ring, Local_param_s2hat.last_ring, submap_size, local_map_s2hat, 
        Local_param_s2hat.gangrank, Local_param_s2hat.gangsize, Local_param_s2hat.gangroot, Local_param_s2hat.gangcomm);
    // 1 map given, 0 is the number of the map
}




void free_covariance_matrix(double ** covariance_matrix_NxN, int lmax){
    int ell_value;
    for (ell_value=0; ell_value<lmax+1; ell_value++){
        free(covariance_matrix_NxN[ell_value]);
    }
    free(covariance_matrix_NxN);    
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

void free_s2hat_parameters_struct(S2HAT_parameters *S2HAT_params){

    free_s2hat_GLOBAL_parameters_struct(S2HAT_params->Global_param_s2hat);
    if (S2HAT_params->Local_param_s2hat->gangrank != -1)
        free_s2hat_LOCAL_parameters_struct(S2HAT_params->Local_param_s2hat);

    free(S2HAT_params);
}

