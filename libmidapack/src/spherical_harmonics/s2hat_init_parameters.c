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




int get_main_s2hat_global_parameters(int nside, int *mask_binary, s2hat_pixeltype *pixelization_scheme, s2hat_scandef *scan_sky_structure_pixel, s2hat_pixparameters *pixpar){
    /* Get s2hat structure which must be distributed by all processors
    All processors must have those same s2hat structure 
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    
    double *mask;
    // int *mask_binary;
    int f_sky=0;
    long npix = nside*nside*12;

    int pixchoice = PIXCHOICE_HEALPIX;   /* Choice  of convention to use for the pixelization scheme, here HEALPIX is chosen, 
    see https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/CparsAndStruct.html for more details*/

    pixpar->par1 = nside; /* NSIDE if HEALPIX convention is used for the transforms, 
    see https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/CsetPixelization.html */


    if (mask_binary == NULL){
        mask_binary = (int *) malloc( npix*sizeof(int));
        int index;
        for (index=0;index<npix;index++){
            mask_binary[index] = 1;
        }
        f_sky = npix;
    }

    set_pixelization(pixchoice, *pixpar, pixelization_scheme);   /* Set C pixelization structure */
    mask2scan(mask_binary, *pixelization_scheme, scan_sky_structure_pixel); /* Set scan pixelization : s2hat structure containing all the info about sky coverage needed by the transforms */
    /* Scan pixelization will be used to avoid doing unecessary calculations */
        
    // printf( "F_sky = %f%%\n", (double)f_sky/(double)npix*100.); fflush(stdout);
    // free(mask_binary); // Problem ?
    return 0;
}


int init_s2hat_global_parameters(Files_path_WIENER_FILTER *Files_WF_struct, int *mask_binary, int lmax, S2HAT_GLOBAL_parameters *Global_param_s2hat){
    /* Create s2hat structure of global parameters of s2hat, which must be distributed to all processors
    All processors must have those same s2hat structure 
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    
    s2hat_pixeltype pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel;
    s2hat_pixparameters pixpar;
    // char *maskfile_path = Files_WF_struct.maskfile_path;
    // bool use_mask_file = Files_WF_struct.use_mask_file;
    int nside = Files_WF_struct->nside;
    
    get_main_s2hat_global_parameters(nside, mask_binary, &pixelization_scheme, &scan_sky_structure_pixel, &pixpar);

    Global_param_s2hat->pixelization_scheme = pixelization_scheme;
    Global_param_s2hat->scan_sky_structure_pixel = scan_sky_structure_pixel;
    Global_param_s2hat->pixpar = pixpar;
    
    Global_param_s2hat->nside = nside;
    Global_param_s2hat->nlmax = lmax-1; // S2HAT will generate alms between 0 and nlmax included, so we have to give lamx decreased by 1 to get the lmax requested
    Global_param_s2hat->nmmax = lmax-1;

    return 0;
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

    if (number_ranks_s2hat < initsize){
        if (initrank==0)
            printf("--- Creating MPI subset because %d total MPI processes vs %d processes needed by S2HAT \n", initsize, number_ranks_s2hat);

        mpi_create_subset(number_ranks_s2hat, initcomm, &s2hat_comm);
        s2hat_rank = -1;
        s2hat_size = 0;
        initroot = 0;
        // if (s2hat_comm == MPI_COMM_NULL){
        if (initrank < number_ranks_s2hat){
            MPI_Comm_rank ( s2hat_comm, &s2hat_rank ); 
            MPI_Comm_size ( s2hat_comm, &s2hat_size );
        }
    }
    else{
        if (initrank==0)
            printf("--- Using MPI default structure because %d total MPI processes sufficient for %d processes needed by S2HAT \n", initsize, number_ranks_s2hat);
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


int init_s2hat_local_parameters_struct(S2HAT_GLOBAL_parameters *Global_param_s2hat, int nstokes, S2HAT_LOCAL_parameters *Local_param_s2hat){
    /* Create s2hat structure of local parameters of s2hat, which will differ for all processors
    
    !!! BEWARE : MPI structure is assumed to already have been assigned using init_MPI_struct_s2hat_local_parameters !!!

    Those local parameters are used to improve the computation efficiency
    Full documentation here : https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/S2HATdocs.html */ 
    int *pixel_numbered_ring;
    s2hat_pixeltype pixelization_scheme = Global_param_s2hat->pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel = Global_param_s2hat->scan_sky_structure_pixel;
    int nlmax = Global_param_s2hat->nlmax;
    int nmmax = Global_param_s2hat->nmmax;
    int number_pixel_total = 12*((Global_param_s2hat->nside)*(Global_param_s2hat->nside));
    int first_pixel_number_ring, last_pixel_number_ring, number_pixels_local;

    int i, j, plms=0, nmvals, first_ring, last_ring, map_size;
    int first_pixel_south_hermisphere;
    int correction_mid_ring = 0; // Correction to take into account the middle ring doesn't contribute twitce

    long int nplm;
    int *mvals;
    
    if (Local_param_s2hat->gangrank != -1){
        get_local_data_sizes( plms, pixelization_scheme, scan_sky_structure_pixel, nlmax, nmmax, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize,
                &nmvals, &first_ring, &last_ring, &map_size, &nplm, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);
        
        // int nmvals_2;
        // nmvals_2 = nummvalues(Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, nmmax);
        // fflush(stdout);
        /* Estimates size of local data object, and set first_ring, last_ring, map_size
        see more info on https://apc.u-paris.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/docs/Cmanual/CgetLocalDataSizes.html */

        mvals = (int *) calloc( nmvals, sizeof( int));    /// TO FREE LATER !!!!!!
        find_mvalues( Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, nmmax, nmvals, mvals);

        Local_param_s2hat->first_pixel_number = Global_param_s2hat->pixelization_scheme.fpix[first_ring-1];
        Local_param_s2hat->last_pixel_number = Global_param_s2hat->pixelization_scheme.fpix[last_ring-1] + Global_param_s2hat->pixelization_scheme.nph[last_ring-1]-1;
        // The last expression correspond to the indice of the first pixel of the last ring (given by fpix[last_ring - 1]),
        // to which we add the number of pixels within the last ring (given by nph[last_ring - 1])
        // BEWARE : we assume the last_pixel_number to be INCLUDED in the last ring

        first_pixel_number_ring = Global_param_s2hat->pixelization_scheme.fpix[first_ring-1];
        // Getting first pixel number of ring probed by local proc, in S2HAT convention
        last_pixel_number_ring = Global_param_s2hat->pixelization_scheme.fpix[last_ring-1] + Global_param_s2hat->pixelization_scheme.nph[last_ring-1]; 
        // Getting last pixel number of ring probed by local proc, in S2HAT convention
        number_pixels_local = last_pixel_number_ring - first_pixel_number_ring;

        if (map_size){
            pixel_numbered_ring = (int *)malloc(map_size*nstokes * sizeof(int));

            for(i=0; i<number_pixels_local; i++)
            {
                for (j=0; j<nstokes; j++)
                    pixel_numbered_ring[i + j*map_size] = i + first_pixel_number_ring + j*number_pixel_total;
            }
            first_pixel_south_hermisphere =  number_pixel_total -  last_pixel_number_ring;
            if (last_ring == Global_param_s2hat->pixelization_scheme.nringsall)
            {
                first_pixel_south_hermisphere += Global_param_s2hat->pixelization_scheme.nph[last_ring-1];
                // If last_ring is equal to the nringsall, which includes the northern rings and the equatorial ring,
                // Then the last ring corresponds to the equatorial ring, which has already been viewed, and we need to avoid redistributing its pixel numberings

                correction_mid_ring = Global_param_s2hat->pixelization_scheme.nph[last_ring-1];
            }

            for(i=0; i<number_pixels_local-correction_mid_ring; i++)
            {
                for (j=0; j<nstokes; j++)
                    pixel_numbered_ring[ i + number_pixels_local + j*map_size] = i + first_pixel_south_hermisphere + j*number_pixel_total;
            }
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
    // S2HAT_GLOBAL_parameters Global_param_s2hat;
    // S2HAT_GLOBAL_parameters *Global_param_s2hat = (S2HAT_GLOBAL_parameters *) malloc( 1 * sizeof(S2HAT_GLOBAL_parameters));
    // init_s2hat_global_parameters(*Files_WF_struct, mask_binary, Files_WF_struct->lmax_Wiener_Filter, Global_param_s2hat);

    init_s2hat_global_parameters(Files_WF_struct, mask_binary, Files_WF_struct->lmax_Wiener_Filter, &(S2HAT_params->Global_param_s2hat));
    // Initialization of Global_param_s2hat structure, for sky pixelization scheme and lmax_WF choice
    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params->Global_param_s2hat);

    int rank;
    MPI_Comm_rank( world_comm, &rank);

    // S2HAT_LOCAL_parameters *Local_param_s2hat = (S2HAT_LOCAL_parameters *) malloc( 1 * sizeof(S2HAT_LOCAL_parameters));
    // init_MPI_struct_s2hat_local_parameters(Local_param_s2hat, Global_param_s2hat->scan_sky_structure_pixel.nringsobs, world_comm);
    init_MPI_struct_s2hat_local_parameters(&(S2HAT_params->Local_param_s2hat), Global_param_s2hat->scan_sky_structure_pixel.nringsobs, world_comm); 
    
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params->Local_param_s2hat);
    if (Local_param_s2hat->gangrank >= 0){
        init_s2hat_local_parameters_struct(Global_param_s2hat, nstokes, Local_param_s2hat);
        }
    // Initialization of Local_param_s2hat structure, including MPI parameters, first/last rings studied, size of pixels cut sky per rank, etc. -- see Wiener filter extension directory for more details
    
    int first_ring = Local_param_s2hat->first_ring;

    // S2HAT_params->Global_param_s2hat = Global_param_s2hat;
    // S2HAT_params->Local_param_s2hat = Local_param_s2hat;
    S2HAT_params->Files_WF_struct = *Files_WF_struct;
    // Initialization of final superstructure S2HAT_params

    S2HAT_params->size_alm = (Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax;
    S2HAT_params->nstokes = nstokes;
    return 0;
}
