#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
#include <lapacke.h>
// #include "fitsio.h"
// #include <mpi.h>
#include <unistd.h>
// #include "s2hat.h"
#include <fitsio.h>
#include <chealpix.h>

#include "mappraiser.h"


void apply_Wiener_filter_pixel(int nside_, int lmax_, int nstokes_, double *CMB_map_ring_, double *CMB_map_output_, double *c_ells_, int number_correlations_, double *mask_binary_, MPI_Comm worldcomm)
{

    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/SO_wH.fits";
    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/No_Mask.fits";
    // char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_multiplesim_white_noise_1.fits";
    // char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_band_limited_1_degree.fits";
    // char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_band_limited_1_degree_bis.fits";
    // char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_band_limited_1024_0.fits";
    char *path_CMB_map = "/Users/mag/Documents/PHD1Y/Space_Work/Inference_Sampling/map_files/Map_band_limited_1024_0.fits";
    // int *mask_binary;
    double *CMB_map;

    int f_sky, npix;
    int i_index, ncomp;
    int index;

    int nside = 512;
    // int lmax = 1535;
    // int lmax = 1024;
    int lmax = 2*nside+2; //+10; //3*nside-1 ;//1500;//025;
    // int nstokes = 3;
    npix = 12*nside*nside;

    // S2HAT_GLOBAL_parameters *Global_param_s2hat;
    // S2HAT_LOCAL_parameters *Local_param_s2hat;
    
    int rank, nprocs;
    MPI_Comm gangcomm;
    


    // MPI_Init( &argc, &argv);
    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &nprocs);

    printf("NEW SIM ################################################################################################################################################ \n");
    printf("Initializing MPI %d %d \n", rank, nprocs);
    fflush(stdout);

    int gangrank = rank;
    int gangsize = nprocs;
    gangcomm = MPI_COMM_WORLD;


    printf("Initializing Global_param_s2hat \n");
    // Files_path_WIENER_FILTER Files_path_WF_struct;
    // Files_path_WF_struct = malloc( 1 * sizeof(Files_path_WIENER_FILTER));
    char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits"; //// TO PUT !!!!
    int number_correlations = 4; // 6; // TO MODIFY LATER !!!!!!!
    int nstokes = 3;
    // printf("Getting into init WF \n"); fflush(stdout);    
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    printf("--- Test init %d %d \n", (S2HAT_params.Files_WF_struct).lmax_Wiener_Filter, lmax); fflush(stdout);
    int i;

    double *mask_binary=NULL;// = calloc(12*nside*nside, sizeof(int));
    // int ga0, number_of_pixels_one_ring = 8;//6*512;
    // for(i=0; i<number_of_pixels_one_ring; i++)
    //     mask_binary[i+gap*number_of_pixels_one_ring]=1;
    // mask_binary = calloc(npix,sizeof(int));
    // for(i=10; i<npix/15; i++)
    //     mask_binary[i]=1;

    printf("Initializing S2HAT_param \n"); fflush(stdout);
    
     // = malloc(1*sizeof(S2HAT_parameters));
    // printf("^^^^^^^ %d ----- RANK for MPI_subgroup before init S2HAT_param \n", gangrank); fflush(stdout);
    
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, gangcomm);
    printf("Finish initializing S2HAT_param !! \n"); fflush(stdout);

    // printf("--- Test2 init %d %d \n", Files_path_WF_struct.lmax_Wiener_Filter, lmax);
    // printf("--- Test3 init %d \n", S2HAT_params.Global_param_s2hat->nlmax); fflush(stdout);
    // S2HAT_parameters *S2HAT_params = &S2HAT_params_;

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);
    int first_ring = Local_param_s2hat->first_ring;
    // printf("Test verif3 : %d %ld \n", first_ring, Global_param_s2hat->pixelization_scheme.fpix[first_ring-1]); fflush(stdout);
    // printf("###### Test4 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    // printf("###### Test12 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    printf("Getting CMB maps !!!\n");
    fflush(stdout);
    // printf("###### Test10 - %ld \n", Local_param_s2shat->pixel_numbered_ring[0]);
    npix = 12*nside*nside;
    printf("###### Test nside - %d %d \n", nside, npix);
    fflush(stdout);
    CMB_map = (double *) malloc( 3*npix*sizeof(double));
    printf("###### ??? \n");
    fflush(stdout);
    printf("###### Test13 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    fflush(stdout);
    read_TQU_maps(nside, CMB_map, path_CMB_map, 3);
    printf("###### Test14 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    // printf("Reading map - rank %d \n", gangrank);
    // fflush(stdout);
    printf("CMB map read !!!\n");
    fflush(stdout);
    printf("###### Test11 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    // double *CMB_map_temp  = CMB_map+0; // &(CMB_map[npix]);
    double *CMB_map_temp  = (double *)malloc(nstokes*npix*sizeof(double));
    memcpy(CMB_map_temp, CMB_map, nstokes*npix*sizeof(double));

    // for(i_index=0; i_index<npix; i_index++){
    //     CMB_map_temp[i_index] = 0;
    // }
    // for(i_index=npix; i_index<3*npix; i_index++){
    //     CMB_map_temp[i_index] = 0;
    // }

    // double *new_CMB_map;
    // new_CMB_map = (double *) malloc( nstokes*npix*sizeof(double));
    // for(i_index=0; i_index<npix; i_index++){
        // for (ncomp=0; ncomp<nstokes; ncomp++){
        //     CMB_map_temp[i_index*nstokes + ncomp] = CMB_map[ncomp*npix + i_index];
        // }
    // }
    printf("###### Test9 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    printf("Minute test1 - rank %d - %f %f %f \n", gangrank, CMB_map[0], CMB_map[npix], CMB_map[2*npix]); fflush(stdout);
    printf("Minute test1_bis - rank %d - %d \n", gangrank, Global_param_s2hat->nlmax); fflush(stdout);
    // printf("Minute test1 - rank %d - %f \n", gangrank, CMB_map_polar[0]);
    // printf("Minute test1_bis - rank %d - %d \n", gangrank, Global_param_s2hat->nlmax);
    // free(CMB_map);
    // printf("Minute test2 - rank %d - %f \n", gangrank, new_CMB_map[0]);
    // printf("Changing map - rank %d \n", gangrank);
    // fflush(stdout);

    double *local_map_pix;//, *local_map_pix_E, *local_map_pix_B;
    local_map_pix = (double *) calloc( 3*Local_param_s2hat->map_size, sizeof(double));
    printf("###### Test8 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]); fflush(stdout);
    printf("###### TEST2 --- MPI struct - %d \n", Local_param_s2hat->gangcomm); fflush(stdout);
    printf("###### CMB Temp - %f \n", CMB_map_temp[0]); fflush(stdout);
    // printf("######2 CMB Temp - %f \n", CMB_map_temp[2*npix]); fflush(stdout);
    Global_param_s2hat->pixelization_scheme;
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_temp, local_map_pix, &S2HAT_params);
    // distribute_full_sky_map_into_local_maps_S2HAT(CMB_map, local_map_pix, S2HAT_params);
    // local_map_pix_E = local_map_pix;
    // local_map_pix_B = local_map_pix+Local_param_s2hat->map_size;
    // distribute_map(Global_param_s2hat->pixelization_scheme, 1, 0, nstokes, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size,
	// 	      local_map_pix, CMB_map, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);
    printf("###### Test5 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    printf("Tiny test - rank %d - %f \n", gangrank, CMB_map_temp[0]); fflush(stdout);

    s2hat_dcomplex *local_alm;
    // s2hat_dcomplex *local_alm_E;
    // s2hat_dcomplex *local_alm_B;
    // local_alm_E = (s2hat_dcomplex *) malloc( (Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax*sizeof(s2hat_dcomplex));
    // local_alm_B = (s2hat_dcomplex *) malloc( (Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax*sizeof(s2hat_dcomplex));
    // local_alm = (s2hat_dcomplex *) malloc( (2*Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax*sizeof(s2hat_dcomplex));
    local_alm = (s2hat_dcomplex *) malloc( ((3*(Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax))*sizeof(s2hat_dcomplex));

    printf("## Test dims local_alm : %d, %d, %d \n", Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, nstokes);

    printf("Little test before pix2alm - rank %d - %f %f \n", gangrank, local_map_pix[0], CMB_map_temp[0]);
    // printf("Little test before pix2alm - rank %d - %f %f \n", gangrank, local_map_pix[0], CMB_map[0]);
    fflush(stdout);


    int max_size_test = 20;

    // printf("Map_pix after distribute_map - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    printf("Map_pix after distribute_map - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map_temp[0]);
    for (index=1;index<max_size_test;index++){
        printf("- %f %f -", local_map_pix[index], CMB_map_temp[index]);
    }
    printf(" \n");
    printf("Map_pix -2- after distribute_map - rank %d - %f %f -", gangrank, local_map_pix[Local_param_s2hat->map_size], CMB_map_temp[npix]);
    for (index=Local_param_s2hat->map_size;index<max_size_test+Local_param_s2hat->map_size;index++){
        printf("- %f %f -", local_map_pix[index], CMB_map_temp[index]);
    }
    printf(" \n");
    fflush(stdout);
    printf("###### Test6 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    int nmaps = 1; // We only provide 1 input set of alm coefficient
    // int lda = nstokes;
    int lda = Global_param_s2hat->nlmax;
    
    double *local_w8ring;
    int i_ring;

    int spin;

    
    // local_w8ring = (double *) malloc( nstokes*(Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1)*sizeof(double));
    // local_w8ring = (double *) malloc( nstokes*(Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1)*sizeof(double));
    // for( i_ring=0; i_ring< nstokes*(Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1); i_ring++)
    //         local_w8ring[i_ring]=1.;
    
    printf("Test 11! \n");
    fflush(stdout);
    // s2hat_map2alm(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
    //             Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
    //             Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm, 
    //             Local_param_s2hat->nplm, NULL,
    //             Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0

    // spin=2;
    // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);

    printf("----- Test nstokes : %d \n", S2HAT_params.nstokes);
    apply_pix2alm(local_map_pix, local_alm, &S2HAT_params);
    // spin=2;
    // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    

    // s2hat_map2alm(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
    //             Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
    //             Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm, 
    //             Local_param_s2hat->nplm, NULL,
    //             Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    
    // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm_E,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    // printf("Test between the 2 pix2alm - rank %d \n", gangrank);
    // fflush(stdout);
    // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix+npix, lda, local_alm_B,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    printf("###### Test7 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    printf("Test1 between the 2 pix2alm - rank %d \n", gangrank);
    fflush(stdout);
    // free(local_w8ring);
    printf("Test2 between the 2 pix2alm - rank %d \n", gangrank);
    fflush(stdout);
    // apply_pix2alm(local_map_pix, local_alm, S2HAT_params); //////////////////////////////////////////////////////////////////////
    printf("###### Test3 - %ld \n", Local_param_s2hat->pixel_numbered_ring[0]);
    printf("Map_pix after first apply_pix2alm - rank %d - %f %f - \n", gangrank, local_map_pix[0], CMB_map_temp[0]); fflush(stdout);

    // printf("Map_pix after 1st pix2alm - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map_polar[0]);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_map_pix[index], CMB_map_polar[index]);
    // }
    // printf(" \n");
    
    // printf("Local_alm after 1st pix2alm - rank %d - %f %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // }
    // printf(" \n");

    

    // printf("Test post-free --- %d !\n", Local_param_s2hat->gangrank);
    // fflush(stdout);

    // double *local_map_pix_2;
    // local_map_pix_2 = (double *) malloc( nstokes*Local_param_s2hat->map_size*sizeof(double));

    // lda = nstokes;
    printf("Tset 0 -- %d \n", Local_param_s2hat->map_size); fflush(stdout);
    // Local_param_s2hat->mvals[Local_param_s2hat->nmvals-1];
    printf("Tset 0F \n"); fflush(stdout);
    // s2hat_alm2map(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
    //             Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
    //             Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, 
    //             local_alm, Local_param_s2hat->nplm, NULL, 
    //             Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    apply_alm2pix(local_map_pix, local_alm, &S2HAT_params);
    // spin=2;
    // s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //                 Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
    //                 Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
    //                 Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);



    // s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //                 Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
    //                 Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm_E,
    //                 Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    // s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //                 Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
    //                 Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix+npix, lda, local_alm_B,
    //                 Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    
    // apply_alm2pix(local_map_pix, local_alm, S2HAT_params); //////////////////////////////////////////////////////////////////////
    printf("Map_pix after first apply_alm2pix - rank %d - %f %f - \n", gangrank, local_map_pix[0], CMB_map_temp[0]); fflush(stdout);
    printf("Map_pix pix npix after first apply_alm2pix - rank %d - %f %f - \n", gangrank, local_map_pix[npix], CMB_map_temp[npix]); 
    fflush(stdout);
    // printf("Map_pix after 1st alm2pix - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map_polar[0]);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_map_pix[index], CMB_map_polar[index]);
    // }
    // printf(" \n");

    // printf("Local_alm after 1st alm2pix - rank %d - re %f im %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // }
    // printf(" \n");

    // fflush(stdout);
    double *full_sky_map_2;
    full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
    gather_map(local_map_pix, full_sky_map_2, nstokes, &S2HAT_params);

    // char filename_save[80];
    // if (Local_param_s2hat->gangrank == 0)
    // {
    //     printf("Saving Stokes components maps \n"); fflush(stdout);
    //     char *Stokes_param = "IQU";
    //     int j;
        
    //     float *full_sky_map_Stokes = (float *) calloc(npix,sizeof(float));

    //     npix = 12*nside*nside;
    //     for(i=0;i<nstokes; i++)
    //     {
    //         for(j=0; j<npix; j++)
    //             full_sky_map_Stokes[j] = full_sky_map_2[j + i*npix];
    //         // printf("--- %d -- Saving Stokes components maps \n", i); fflush(stdout);
    //         sprintf(filename_save, "/global/cscratch1/sd/mag/WF_work/WF_Tests/mapfile_BL1024_%d.fits", i);
    //         // printf("--- %d -- Saving Stokes components maps \n", i); fflush(stdout);
    //         write_healpix_map(full_sky_map_Stokes, nside, filename_save, 0, "C");
    //     }
    //     free(full_sky_map_Stokes); fflush(stdout);
    //     printf("Stokes component maps saved ! \n");

    // }

    // double *local_map_pix_2;
    // local_map_pix_2 = (double *) calloc( nstokes*Local_param_s2hat->map_size, sizeof(double));
    // distribute_map(Global_param_s2hat->pixelization_scheme, 1, 0, nstokes, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size,
	// 	      local_map_pix_2, full_sky_map_2, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);

    // apply_pix2alm(local_map_pix_2, local_alm, S2HAT_params); //////////////////////////////////////////////////////////////////////

    // // printf("Map_pix after 2nd pix2alm - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map_polar[0]);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_map_pix[index], CMB_map_polar[index]);
    // // }
    // // printf(" \n");
    
    // // printf("Local_alm after 2nd pix2alm - rank %d - re %f im %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // // }
    // // printf(" \n");
    // // fflush(stdout);

    // apply_alm2pix(local_map_pix, local_alm, S2HAT_params); //////////////////////////////////////////////////////////////////////
    
    // // printf("Map_pix after 2nd alm2pix - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map_polar[0]);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_map_pix[index], CMB_map_polar[index]);
    // // }
    // // printf(" \n");
    
    // // printf("Local_alm after 2nd alm2pix - rank %d - re %f im %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // // }
    // // printf(" \n");

    // double *full_sky_map;
    // full_sky_map = (double *) malloc(nstokes*npix*sizeof(double));
    // gather_map(local_map_pix, full_sky_map, nstokes, S2HAT_params);
    

    if (Local_param_s2hat->gangrank == 0){
    printf("//////////////////////// Comparison new map vs old \n");
    fflush(stdout);

    printf("From pixel 0 to %d - %f %f -", max_size_test, full_sky_map_2[0], CMB_map_temp[0]);
    for (index=0;index<max_size_test;index++){
        printf("- %f %f -", full_sky_map_2[index], CMB_map_temp[index]);
        }
    printf(" \n");
    fflush(stdout);

    // int new_pixel = npix;
    // printf("From pixel %d to %d - %f %f -", new_pixel, new_pixel+max_size_test, full_sky_map_2[new_pixel], CMB_map_temp[new_pixel]);
    // for (index=new_pixel;index<new_pixel+max_size_test;index++){
    //     printf("- %f %f -", full_sky_map_2[index], CMB_map_temp[index]);
    //         }
    // printf(" \n");
    // fflush(stdout);

    // new_pixel = 2*npix;
    // printf("From pixel %d to %d - %f %f -", new_pixel, new_pixel+max_size_test, full_sky_map_2[new_pixel], CMB_map[new_pixel]);
    // for (index=new_pixel;index<new_pixel+max_size_test;index++){
    //     printf("- %f %f -", full_sky_map_2[index], CMB_map[index]);
    //     }
    // printf(" \n");

    double average_relative_error_T_0_1 = 0;
    double average_relative_error_Q_0_1 = 0;
    double average_relative_error_U_0_1 = 0;
    double max_average_relative_error_T_0_1 = 0;
    double max_average_relative_error_Q_0_1 = 0;
    double max_average_relative_error_U_0_1 = 0;
    // double average_relative_error_T_1_2 = 0;
    // double average_relative_error_Q_1_2 = 0;
    // double average_relative_error_U_1_2 = 0;
    // double average_relative_error_T_0_2 = 0;
    // double average_relative_error_Q_0_2 = 0;
    // double average_relative_error_U_0_2 = 0;
    int size_total_probed = npix;
    // size_total_probed = 10;

    int init_index = 0;
    int number_0 = 0;

    double pixel_value_T;
    double pixel_value_Q;
    double pixel_value_U;

    index = init_index;
    
    max_average_relative_error_T_0_1 = fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
    max_average_relative_error_Q_0_1 = fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
    max_average_relative_error_U_0_1 = fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);
    for(index=init_index;index<size_total_probed+init_index;index++){
        pixel_value_T = fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
        average_relative_error_T_0_1 += pixel_value_T;
        if (max_average_relative_error_T_0_1 < pixel_value_T)
            max_average_relative_error_T_0_1 = pixel_value_T;
        // average_relative_error_T_1_2 += fabs((full_sky_map[index]-full_sky_map_2[index])/full_sky_map[index]);
        // average_relative_error_T_0_2 += fabs((full_sky_map[index]-CMB_map[index])/CMB_map[index]);
        // printf("Test %d ---- %f \t", index, fabs((full_sky_map[index]-CMB_map[index])/CMB_map[index]));
        if(CMB_map[index] == 0.000000){
            number_0++;
        }
        if(nstokes>1){
            pixel_value_Q = fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
            average_relative_error_Q_0_1 += pixel_value_Q;
            if (max_average_relative_error_Q_0_1 < pixel_value_Q)
                max_average_relative_error_Q_0_1 = pixel_value_Q;
            if(CMB_map[index+npix] == 0.000000){
                number_0++;
            }
            // average_relative_error_Q_1_2 += fabs((full_sky_map[index+npix]-full_sky_map_2[index+npix])/full_sky_map[index+npix]);
            // average_relative_error_Q_0_2 += fabs((full_sky_map[index+npix]-CMB_map[index+npix])/CMB_map[index+npix]);
            // average_relative_error_Q += fabs((full_sky_map[index+npix]-CMB_map[index+npix])/CMB_map[index+npix]);
            // printf("Test2 %d ---- %f \t", index, (full_sky_map[index+npix]-CMB_map[index+npix])/CMB_map[index+npix]);
            pixel_value_U = fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);
            average_relative_error_U_0_1 += pixel_value_U;
            if (max_average_relative_error_U_0_1 < pixel_value_U)
                max_average_relative_error_U_0_1 = pixel_value_U;
            if(CMB_map[index+2*npix] == 0.000000){
                number_0++;
            }
            // average_relative_error_U_1_2 += fabs((full_sky_map[index+2*npix]-full_sky_map_2[index+2*npix])/full_sky_map[index+2*npix]);
            // average_relative_error_U_0_2 += fabs((full_sky_map[index+2*npix]-CMB_map[index+2*npix])/CMB_map[index+2*npix]);
            // average_relative_error_U += fabs((full_sky_map[index+2*npix]-CMB_map[index+2*npix])/CMB_map[index+2*npix]);
            // printf("Test3 %d ---- %f \t", index, (full_sky_map[index+2*npix]-CMB_map[index+2*npix])/CMB_map[index+2*npix]);
        }
    }
    printf("\n");
    printf("Average 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            average_relative_error_T_0_1/(size_total_probed), average_relative_error_Q_0_1/(size_total_probed), average_relative_error_U_0_1/(size_total_probed));
    printf("Maximum 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            max_average_relative_error_T_0_1, max_average_relative_error_Q_0_1, max_average_relative_error_U_0_1);
    // if(nstokes>1){
    //     printf("Average 1-2 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
    //             average_relative_error_T_1_2/(size_total_probed), average_relative_error_Q_1_2/(size_total_probed), average_relative_error_U_1_2/(size_total_probed));
    //     printf("Average 0-2 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
    //             average_relative_error_T_0_2/(size_total_probed), average_relative_error_Q_0_2/(size_total_probed), average_relative_error_U_0_2/(size_total_probed));
    // }
    printf("Number of 0 : %d \n", number_0);
    fflush(stdout);
    }

    // printf("Test 6 ! --- %d \n", Local_param_s2hat->gangrank);
    // fflush(stdout);
    // printf("Test pointer : %f \n", *new_CMB_map);
    // fflush(stdout);
    // // free(local_map_pix_2);
    // // free(new_CMB_map);
    // printf("Pointer adresses : %p %f - %p %f  -- %d \n", local_map_pix, *local_map_pix, new_CMB_map, *new_CMB_map, local_map_pix - new_CMB_map);
    // fflush(stdout);
    // free(new_CMB_map);
    // printf("Test 7 !\n");
    // fflush(stdout);
    // printf("Pointer adresses : %p - %p \n", local_map_pix, new_CMB_map);
    // fflush(stdout);
    // // printf("Test pointer : %f \n", *new_CMB_map);
    // fflush(stdout);
    // printf("Test pointer 2 : %f \n", *local_map_pix);
    // fflush(stdout);

    printf("Test 8 !\n");
    fflush(stdout);
    printf("###### Test2 - %ld %d \n", Local_param_s2hat->pixel_numbered_ring[0],  Local_param_s2hat->map_size * S2HAT_params.nstokes);
    // free(local_map_pix);
    printf("Test 3 !\n");
    fflush(stdout);
    free(CMB_map);
    free(CMB_map_temp);

    
    printf("Test 9 !\n"); fflush(stdout);
    free_s2hat_parameters_struct(&S2HAT_params);

    printf("Test 10 !\n"); fflush(stdout);
    // MPI_Finalize();

    printf("Test finish ! - %d \n", gangrank);
    fflush(stdout);

}


void apply_Wiener_filter_pixel_wrong(int nside_, int lmax_, int nstokes_, double *CMB_map_ring_, double *CMB_map_output_, double *c_ells_, int number_correlations_, double *mask_binary_, MPI_Comm worldcomm)
{
// int main(int argc, char** argv){
    
    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/SO_wH.fits";
    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/No_Mask.fits";
    // char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_multiplesim_white_noise_1.fits";
    // char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_band_limited_1_degree.fits";
    // char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_band_limited_1_degree_bis.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_band_limited_1024_0.fits";
    char *path_CMB_map = "/Users/mag/Documents/PHD1Y/Space_Work/Inference_Sampling/map_files/Map_band_limited_1024_0.fits";
    // int *mask_binary;
    double *CMB_map;

    int f_sky, npix;
    int i_index, ncomp;
    int index;

    int nside = 512;
    // int lmax = 1535;
    // int lmax = 1024;
    int lmax = 2*nside; //+10; //3*nside-1 ;//1500;//025;
    // int nstokes = 3;
    npix = 12*nside*nside;

    // S2HAT_GLOBAL_parameters *Global_param_s2hat;
    // S2HAT_LOCAL_parameters *Local_param_s2hat;
    
    int rank, nprocs;
    MPI_Comm gangcomm;
    


    // MPI_Init( &argc, &argv);
    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &nprocs);

    printf("Initializing MPI %d %d \n", rank, nprocs);
    fflush(stdout);

    int gangrank = rank;
    int gangsize = nprocs;
    gangcomm = MPI_COMM_WORLD;


    printf("Initializing Global_param_s2hat \n");
    Files_path_WIENER_FILTER Files_path_WF_struct;

    char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits"; //// TO PUT !!!!
    int number_correlations = 6; // TO MODIFY LATER !!!!!!!
    int nstokes = 3; //3;
    // printf("Getting into init WF \n"); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    
    int *mask_binary=NULL;// = calloc(12*nside*nside, sizeof(int));
    int i, ga0, number_of_pixels_one_ring = 8;//6*512;
    // for(i=0; i<number_of_pixels_one_ring; i++)
    //     mask_binary[i+gap*number_of_pixels_one_ring]=1;
    // mask_binary = calloc(npix,sizeof(int));
    // for(i=10; i<npix/15; i++)
    //     mask_binary[i]=1;
    
    printf("Initializing S2HAT_param \n"); fflush(stdout);
    
    // printf("^^^^^^^ %d ----- RANK for MPI_subgroup before init S2HAT_param \n", gangrank); fflush(stdout);
    // init_s2hat_parameters_superstruct(Files_path_WF_struct, mask_binary, nstokes, S2HAT_params, gangcomm);
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("Getting CMB maps !!!\n");
    fflush(stdout);

    npix = 12*nside*nside;    
    CMB_map = (double *) malloc( nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map, path_CMB_map, nstokes);
    // printf("Reading map - rank %d \n", gangrank);
    // fflush(stdout);
    printf("CMB map read !!!\n");
    fflush(stdout);

    // double *new_CMB_map;
    // new_CMB_map = (double *) malloc( nstokes*npix*sizeof(double));
    // for(i_index=0; i_index<npix; i_index++){
    //     for (ncomp=0; ncomp<nstokes; ncomp++){
    //         new_CMB_map[i_index*nstokes + ncomp] = CMB_map[ncomp*npix + i_index];
    //     }
    // }
    printf("Minute test1 - rank %d - %f %f %f \n", gangrank, CMB_map[0], CMB_map[npix], CMB_map[2*npix]);
    printf("Minute test1_bis - rank %d - %d \n", gangrank, Global_param_s2hat->nlmax);
    // free(CMB_map);
    // printf("Minute test2 - rank %d - %f \n", gangrank, new_CMB_map[0]);
    // printf("Changing map - rank %d \n", gangrank);
    // fflush(stdout);

    double *local_map_pix;
    local_map_pix = (double *) calloc( nstokes*Local_param_s2hat->map_size, sizeof(double));
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map, local_map_pix, &S2HAT_params);
    // distribute_map(Global_param_s2hat->pixelization_scheme, 1, 0, nstokes, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size,
	// 	      local_map_pix, CMB_map, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);

    printf("Tiny test - rank %d - %f \n", gangrank, CMB_map[0]);

    s2hat_dcomplex *local_alm;
    local_alm = (s2hat_dcomplex *) malloc( nstokes*(Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax*sizeof(s2hat_dcomplex));


    printf("Little test before pix2alm - rank %d - %f %f \n", gangrank, local_map_pix[0], CMB_map[0]);
    fflush(stdout);

    int max_size_test = 20;

    printf("Map_pix after distribute_map - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    for (index=1;index<max_size_test;index++){
        printf("- %f %f -", local_map_pix[index], CMB_map[index]);
    }
    printf(" \n");

    int nmaps = 1; // We only provide 1 input set of alm coefficient
    int lda = nstokes;
    // int lda = Global_param_s2hat->nlmax;
    
    double *local_w8ring;
    int i_ring;

    int spin;

    
    local_w8ring = (double *) malloc( nstokes*(Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1)*sizeof(double));
    for( i_ring=0; i_ring< nstokes*(Local_param_s2hat->last_ring-Local_param_s2hat->first_ring+1); i_ring++)
            local_w8ring[i_ring]=1.;
    
    s2hat_map2alm(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm, 
                Local_param_s2hat->nplm, NULL,
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
            // The NULL argument correspond to precomputed Legendre polynomials, only relevant if plms != 0

    // spin=2;            
    // s2hat_map2alm_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //     Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring,
    //     local_w8ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
    //     Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    free(local_w8ring);
    // apply_pix2alm(local_map_pix, local_alm, S2HAT_params); //////////////////////////////////////////////////////////////////////
    printf("Map_pix after first apply_pix2alm - rank %d - %f %f - \n", gangrank, local_map_pix[0], CMB_map[0]); fflush(stdout);

    // printf("Map_pix after 1st pix2alm - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_map_pix[index], CMB_map[index]);
    // }
    // printf(" \n");
    
    // printf("Local_alm after 1st pix2alm - rank %d - %f %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // }
    // printf(" \n");

    

    // printf("Test post-free --- %d !\n", Local_param_s2hat->gangrank);
    // fflush(stdout);

    // double *local_map_pix_2;
    // local_map_pix_2 = (double *) malloc( nstokes*Local_param_s2hat->map_size*sizeof(double));

    // lda = nstokes;
    printf("Tset 0 -- %d \n", Local_param_s2hat->map_size); fflush(stdout);
    // Local_param_s2hat->mvals[Local_param_s2hat->nmvals-1];
    printf("Tset 0F \n"); fflush(stdout);
    s2hat_alm2map(Local_param_s2hat->plms, Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax, 
                Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps, nstokes, 
                Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, 
                local_alm, Local_param_s2hat->nplm, NULL, 
                Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    // spin=2;
    // s2hat_alm2map_spin( Global_param_s2hat->pixelization_scheme, Global_param_s2hat->scan_sky_structure_pixel, spin, Global_param_s2hat->nlmax, Global_param_s2hat->nmmax,
    //                 Local_param_s2hat->nmvals, Local_param_s2hat->mvals, nmaps,
    //                 Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, local_map_pix, lda, local_alm,
    //                 Local_param_s2hat->gangsize, Local_param_s2hat->gangrank, Local_param_s2hat->gangcomm);
    // apply_alm2pix(local_map_pix, local_alm, S2HAT_params); //////////////////////////////////////////////////////////////////////
    printf("Map_pix after first apply_alm2pix - rank %d - %f %f - \n", gangrank, local_map_pix[0], CMB_map[0]); fflush(stdout);
    // printf("Map_pix after 1st alm2pix - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_map_pix[index], CMB_map[index]);
    // }
    // printf(" \n");

    // printf("Local_alm after 1st alm2pix - rank %d - re %f im %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // }
    // printf(" \n");

    // fflush(stdout);
    double *full_sky_map_2;
    full_sky_map_2 = (double *) malloc(nstokes*npix*sizeof(double));
    gather_map(local_map_pix, full_sky_map_2, nstokes, &S2HAT_params);

    // char filename_save[80];
    // if (Local_param_s2hat->gangrank == 0)
    // {
    //     printf("Saving Stokes components maps \n"); fflush(stdout);
    //     char *Stokes_param = "IQU";
    //     int j;

    //     float *full_sky_map_Stokes = (float *) calloc(npix,sizeof(float));

    //     npix = 12*nside*nside;
    //     for(i=0;i<nstokes; i++)
    //     {
    //         for(j=0; j<npix; j++)
    //             full_sky_map_Stokes[j] = full_sky_map_2[j + i*npix];
    //         // printf("--- %d -- Saving Stokes components maps \n", i); fflush(stdout);
    //         sprintf(filename_save, "/global/cscratch1/sd/mag/WF_work/WF_Tests/mapfile_BL1024_%d.fits", i);
    //         // printf("--- %d -- Saving Stokes components maps \n", i); fflush(stdout);
    //         write_healpix_map(full_sky_map_Stokes, nside, filename_save, 0, "C");
    //     }
    //     free(full_sky_map_Stokes); fflush(stdout);
    //     printf("Stokes component maps saved ! \n");

    // }

    // double *local_map_pix_2;
    // local_map_pix_2 = (double *) calloc( nstokes*Local_param_s2hat->map_size, sizeof(double));
    // distribute_map(Global_param_s2hat->pixelization_scheme, 1, 0, nstokes, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size,
	// 	      local_map_pix_2, full_sky_map_2, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);

    // apply_pix2alm(local_map_pix_2, local_alm, S2HAT_params); //////////////////////////////////////////////////////////////////////

    // // printf("Map_pix after 2nd pix2alm - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_map_pix[index], CMB_map[index]);
    // // }
    // // printf(" \n");
    
    // // printf("Local_alm after 2nd pix2alm - rank %d - re %f im %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // // }
    // // printf(" \n");
    // // fflush(stdout);

    // apply_alm2pix(local_map_pix, local_alm, S2HAT_params); //////////////////////////////////////////////////////////////////////
    
    // // printf("Map_pix after 2nd alm2pix - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_map_pix[index], CMB_map[index]);
    // // }
    // // printf(" \n");
    
    // // printf("Local_alm after 2nd alm2pix - rank %d - re %f im %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // // for (index=1;index<max_size_test;index++){
    // //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // // }
    // // printf(" \n");

    // double *full_sky_map;
    // full_sky_map = (double *) malloc(nstokes*npix*sizeof(double));
    // gather_map(local_map_pix, full_sky_map, nstokes, S2HAT_params);
    

    if (Local_param_s2hat->gangrank == 0){
    printf("//////////////////////// Comparison new map vs old \n");
    fflush(stdout);

    printf("From pixel 0 to %d - %f %f -", max_size_test, full_sky_map_2[0], CMB_map[0]);
    for (index=0;index<max_size_test;index++){
        printf("- %f %f -", full_sky_map_2[index], CMB_map[index]);
        }
    printf(" \n");
    fflush(stdout);

    int new_pixel = npix;
    printf("From pixel %d to %d - %f %f -", new_pixel, new_pixel+max_size_test, full_sky_map_2[new_pixel], CMB_map[new_pixel]);
    for (index=new_pixel;index<new_pixel+max_size_test;index++){
        printf("- %f %f -", full_sky_map_2[index], CMB_map[index]);
            }
    printf(" \n");
    fflush(stdout);

    new_pixel = 2*npix;
    printf("From pixel %d to %d - %f %f -", new_pixel, new_pixel+max_size_test, full_sky_map_2[new_pixel], CMB_map[new_pixel]);
    for (index=new_pixel;index<new_pixel+max_size_test;index++){
        printf("- %f %f -", full_sky_map_2[index], CMB_map[index]);
        }
    printf(" \n");

    double average_relative_error_T_0_1 = 0;
    double average_relative_error_Q_0_1 = 0;
    double average_relative_error_U_0_1 = 0;
    double max_average_relative_error_T_0_1 = 0;
    double max_average_relative_error_Q_0_1 = 0;
    double max_average_relative_error_U_0_1 = 0;
    // double average_relative_error_T_1_2 = 0;
    // double average_relative_error_Q_1_2 = 0;
    // double average_relative_error_U_1_2 = 0;
    // double average_relative_error_T_0_2 = 0;
    // double average_relative_error_Q_0_2 = 0;
    // double average_relative_error_U_0_2 = 0;
    int size_total_probed = npix;
    // size_total_probed = 10;

    int init_index = 0;
    int number_0 = 0;

    index = init_index;
    max_average_relative_error_T_0_1 = fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
    max_average_relative_error_Q_0_1 = fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
    max_average_relative_error_U_0_1 = fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);
    for(index=init_index;index<size_total_probed+init_index;index++){
        average_relative_error_T_0_1 += fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
        if (max_average_relative_error_T_0_1 < fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]))
            max_average_relative_error_T_0_1 = fabs((CMB_map[index]-full_sky_map_2[index])/CMB_map[index]);
        // average_relative_error_T_1_2 += fabs((full_sky_map[index]-full_sky_map_2[index])/full_sky_map[index]);
        // average_relative_error_T_0_2 += fabs((full_sky_map[index]-CMB_map[index])/CMB_map[index]);
        // printf("Test %d ---- %f \t", index, fabs((full_sky_map[index]-CMB_map[index])/CMB_map[index]));
        if(CMB_map[index] == 0.000000){
            number_0++;
        }
        if(nstokes>1){
            average_relative_error_Q_0_1 += fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
            if (max_average_relative_error_Q_0_1 < fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]))
                max_average_relative_error_Q_0_1 = fabs((CMB_map[index+npix]-full_sky_map_2[index+npix])/CMB_map[index+npix]);
            // average_relative_error_Q_1_2 += fabs((full_sky_map[index+npix]-full_sky_map_2[index+npix])/full_sky_map[index+npix]);
            // average_relative_error_Q_0_2 += fabs((full_sky_map[index+npix]-CMB_map[index+npix])/CMB_map[index+npix]);
            // average_relative_error_Q += fabs((full_sky_map[index+npix]-CMB_map[index+npix])/CMB_map[index+npix]);
            // printf("Test2 %d ---- %f \t", index, (full_sky_map[index+npix]-CMB_map[index+npix])/CMB_map[index+npix]);
            average_relative_error_U_0_1 += fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);
            if (max_average_relative_error_U_0_1 < fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]))
                max_average_relative_error_U_0_1 = fabs((CMB_map[index+2*npix]-full_sky_map_2[index+2*npix])/CMB_map[index+2*npix]);

            // average_relative_error_U_1_2 += fabs((full_sky_map[index+2*npix]-full_sky_map_2[index+2*npix])/full_sky_map[index+2*npix]);
            // average_relative_error_U_0_2 += fabs((full_sky_map[index+2*npix]-CMB_map[index+2*npix])/CMB_map[index+2*npix]);
            // average_relative_error_U += fabs((full_sky_map[index+2*npix]-CMB_map[index+2*npix])/CMB_map[index+2*npix]);
            // printf("Test3 %d ---- %f \t", index, (full_sky_map[index+2*npix]-CMB_map[index+2*npix])/CMB_map[index+2*npix]);
        }
    }
    printf("\n");
    printf("Average 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            average_relative_error_T_0_1/(size_total_probed), average_relative_error_Q_0_1/(size_total_probed), average_relative_error_U_0_1/(size_total_probed));
    printf("Maximum 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            max_average_relative_error_T_0_1, max_average_relative_error_Q_0_1, max_average_relative_error_U_0_1);
    // if(nstokes>1){
    //     printf("Average 1-2 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
    //             average_relative_error_T_1_2/(size_total_probed), average_relative_error_Q_1_2/(size_total_probed), average_relative_error_U_1_2/(size_total_probed));
    //     printf("Average 0-2 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
    //             average_relative_error_T_0_2/(size_total_probed), average_relative_error_Q_0_2/(size_total_probed), average_relative_error_U_0_2/(size_total_probed));
    // }
    printf("Number of 0 : %d \n", number_0);
    fflush(stdout);
    }

    // printf("Test 6 ! --- %d \n", Local_param_s2hat->gangrank);
    // fflush(stdout);
    // printf("Test pointer : %f \n", *new_CMB_map);
    // fflush(stdout);
    // // free(local_map_pix_2);
    // // free(new_CMB_map);
    // printf("Pointer adresses : %p %f - %p %f  -- %d \n", local_map_pix, *local_map_pix, new_CMB_map, *new_CMB_map, local_map_pix - new_CMB_map);
    // fflush(stdout);
    // free(new_CMB_map);
    // printf("Test 7 !\n");
    // fflush(stdout);
    // printf("Pointer adresses : %p - %p \n", local_map_pix, new_CMB_map);
    // fflush(stdout);
    // // printf("Test pointer : %f \n", *new_CMB_map);
    // fflush(stdout);
    // printf("Test pointer 2 : %f \n", *local_map_pix);
    // fflush(stdout);

    printf("Test 8 !\n");
    fflush(stdout);

    // free(local_map_pix);
    printf("Test 3 !\n");
    fflush(stdout);
    free(CMB_map);

    
    printf("Test 9 !\n"); fflush(stdout);
    free_s2hat_parameters_struct(&S2HAT_params);

    printf("Test 10 !\n"); fflush(stdout);
    // MPI_Finalize();

    printf("Test finish ! - %d \n", gangrank);
    fflush(stdout);

}

void apply_Wiener_filter_pixel_0(int nside, int lmax, int nstokes, double *CMB_map_ring, double *CMB_map_output, double *c_ells, int number_correlations, double *mask_binary, MPI_Comm worldcomm){
    /* Apply Wiener filter in pixel space to CMB map */

    printf("############################################################ START !!! ############################################################ \n");
    // char *path_CMB_map = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/Map_test_band_limited.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTTTE_128.fits";
    char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_SO.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_cut_latS.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_cut_equator.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128_pt_src.fits";
    

    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_woTE_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_onlyBB_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_onlyEEBB_128.fits";
    // char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_band_limited_1024_0.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTE.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTTTE.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_1024_woTTTE_v2.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128_woTE.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128_only_BB.fits";
    // char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128_only_EEBB.fits";
    // char *path_mask = "/pscratch/sd/m/mag/Masks_files/mask_SO_64.fits";
    // char *path_mask = "/pscratch/sd/m/mag/Masks_files/mask_cut_latS_64.fits";
    // char *path_mask = "/pscratch/sd/m/mag/Masks_files/mask_cut_equator_64.fits";
    char *path_mask = "/pscratch/sd/m/mag/Masks_files/mask_cut_pt_src_64.fits";

    int rank, size;

    // MPI_Init();
    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &size);

    int i, j, ell_value;

    // int nside = 512;
    // int nside_2 = 64;
    // int lmax_2 = 2*nside_2;
    int npix = 12*nside*nside;
    
    double *mask_binary_=NULL;
    // double *mask_binary= (double *)malloc(npix*sizeof(double));
    // read_fits_mask(nside, mask_binary, path_mask, 1);
    
    int nstokes_2 = 3; // T, Q, U
    int number_correlations_2 = 4; // TT, EE, BB, TE
    // int nstokes = 2; //Q, U
    // int number_correlations = 2; // EE, BB


    
    int index, index_2;
    int order_matrix = nstokes;

    // double *CMB_map_ring_2 = (double *) malloc( 3*npix*sizeof(double));
    // read_TQU_maps(nside_2, CMB_map_ring_2, path_CMB_map, 3);


    int max_size_pixels = 40;
    // printf("%d --- Reading first elems of map \n", rank);
    // for (i=0; i<max_size_pixels; i++)
    //     printf("- %f %f %f -", CMB_map_ring_2[i], CMB_map_ring_2[i+npix], CMB_map_ring_2[i+2*npix]);
    // printf("\n");

    printf("%d --- Initializing S2HAT_param \n", rank); fflush(stdout);
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations_2);
    
    
    
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary_, nstokes_2, &S2HAT_params, worldcomm);
    printf("%d --- End of initializing S2HAT_param \n", rank); fflush(stdout);

    double *local_map_pix = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));
    printf("%d --- Distributing S2HAT_param \n", rank); fflush(stdout);
    
    double *CMB_map_ring_2 = (double *)malloc(nstokes*npix*sizeof(double));
    memcpy(CMB_map_ring_2, CMB_map_ring, nstokes*npix*sizeof(double));
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring_2, local_map_pix, &S2HAT_params);
    // free(CMB_map_ring_2);
    printf("%d --- End of distributing S2HAT_param \n", rank); fflush(stdout);
    

    s2hat_dcomplex *local_alm_s2hat;

    local_alm_s2hat = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);
    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);


    printf("%d ----!!!! ***local_map_output*** nstokes %d ; map_size %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    double *local_map_output = (double *)calloc(nstokes*S2HAT_params.Local_param_s2hat.map_size,sizeof(double));


    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);


    memcpy(CMB_map_output, local_map_output, nstokes*npix*sizeof(double));
    // memcpy(CMB_map_output, CMB_map_ring, nstokes*npix*sizeof(double));

    printf("%d --- Free steps 0 ! \n", rank); fflush(stdout);
    // free(local_map_pix);
    printf("%d --- Free steps 1 ! \n", rank); fflush(stdout);
    free(local_map_output);

    printf("%d --- Done !!! \n", rank); fflush(stdout);
}


void apply_Wiener_filter_pixel_test(int nside, int lmax, int nstokes, double *CMB_map_ring, double *CMB_map_output, double *c_ells, int number_correlations, double *mask_binary, MPI_Comm worldcomm){
    /* Apply Wiener filter in pixel space to CMB map */

    char *path_CMB_map = "/global/homes/m/mag/perl_midapack/midapack/mappraiser/test_wiener_filter/Map_test_band_limited_128.fits";

    
    printf("############################################################ START !!! ############################################################ \n"); fflush(stdout);
    int rank, size;

    int npix = 12*nside*nside;
    // MPI_Init();
    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &size);

    // memcpy(CMB_map_output, CMB_map_ring, nstokes*npix*sizeof(double));

    printf("%d --- Initializing S2HAT_param -- nside %d lmax %d nstokes %d \n", rank, nside, lmax, nstokes); fflush(stdout);
    char *c_ell_path = "/global/homes/m/mag/midapack/mappraiser/src/test_wiener_filter/cls_limited_128.fits";
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);

    double *mask_binary_ = NULL;
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary_, nstokes, &S2HAT_params, worldcomm);

    double *local_map_pix = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));

    double *CMB_map_ring_2 = (double *)malloc(nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map_ring_2, path_CMB_map, 3);
    // printf("%d --- Test map_size %d npix %d \n", rank, S2HAT_params.Local_param_s2hat.map_size, npix);
    // distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring, local_map_pix, &S2HAT_params);
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_ring_2, local_map_pix, &S2HAT_params);
    free(CMB_map_ring_2);
    // memcpy(CMB_map_output, local_map_pix, nstokes*npix*sizeof(double));

    s2hat_dcomplex *local_alm_s2hat = (s2hat_dcomplex *)calloc(nstokes*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);
    printf("%d --- map2harmonic-end --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);


    printf("%d ----!!!! ***local_map_output*** nstokes %d ; map_size %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));


    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);
    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);

    memcpy(CMB_map_output, local_map_output, nstokes*npix*sizeof(double));

    printf("%d --- Free steps 0 ! \n", rank); fflush(stdout);
    free(local_map_pix);
    printf("%d --- Free steps 1 ! \n", rank); fflush(stdout);
    free(local_map_output);

    printf("%d --- Done !!! \n", rank); fflush(stdout);
}


void apply_Wiener_filter_pixel_true(int nside, int lmax, int nstokes, double *CMB_map_ring, double *CMB_map_output, double *c_ells, int number_correlations, double *mask_binary, MPI_Comm worldcomm){
    /* Apply Wiener filter in pixel space to CMB map */

    printf("############################################################ START !!! ############################################################ \n"); fflush(stdout);
    int rank, size;

    // MPI_Init();
    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &size);

    int ell_value;
    // int i, j;

    int npix = 12*nside*nside;
    // int index, index_2;
    // int order_matrix = nstokes;

    printf("%d --- Initializing S2HAT_param -- nside %d lmax %d nstokes %d \n", rank, nside, lmax, nstokes); fflush(stdout);
    char *c_ell_path = "";
    S2HAT_parameters S2HAT_params;
    init_files_struct_WF(&(S2HAT_params.Files_WF_struct), nside, lmax, c_ell_path, number_correlations);
    
    init_s2hat_parameters_superstruct(&(S2HAT_params.Files_WF_struct), mask_binary, nstokes, &S2HAT_params, worldcomm);

    double *local_map_pix = (double *)calloc(nstokes*S2HAT_params.Local_param_s2hat.map_size,sizeof(double));

    double *CMB_map_output_copy_0 = (double *)calloc(nstokes*npix,sizeof(double));
    memcpy(CMB_map_output_copy_0, CMB_map_ring, nstokes*npix*sizeof(double));
    distribute_full_sky_map_into_local_maps_S2HAT(CMB_map_output_copy_0, local_map_pix, &S2HAT_params);
    free(CMB_map_output_copy_0);

    s2hat_dcomplex *local_alm_s2hat;

    local_alm_s2hat = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    int i, max_size_pixels=40;
    printf("%d --- Reading first elems of local map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", local_map_pix[i], CMB_map_ring[i]);
    printf("\n");

    // S2HAT_GLOBAL_parameters *Global_param_s2hat = &(S2HAT_params.Global_param_s2hat);
    // S2HAT_LOCAL_parameters *Local_param_s2hat = &(S2HAT_params.Local_param_s2hat);

    printf("%d --- map2harmonic --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    apply_pix2alm(local_map_pix, local_alm_s2hat, &S2HAT_params);
    printf("%d --- map2harmonic-end --- size_alm %d !!! \n", rank, S2HAT_params.size_alm); fflush(stdout);
    // free(local_map_pix);

    double **inverse_covariance_matrix;
    inverse_covariance_matrix = malloc((lmax+1)*sizeof(double *));
    for(ell_value=0; ell_value<lmax+1; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(nstokes*nstokes,sizeof(double));
    }
    printf("%d --- getting inverse of covariance matrix !!! \n", rank); fflush(stdout);

    get_covariance_matrix_block_diagonal(c_ells, number_correlations, inverse_covariance_matrix, &S2HAT_params);

    int index, index_2, order_matrix=nstokes, size_ell_max=10;
    for (ell_value=0;ell_value<size_ell_max; ell_value++){
        printf("\n #### ell_init= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
        printf("\n");
    }
    }
    for(ell_value=0; ell_value<lmax+1; ell_value++){
        get_cholesky_decomposition_inverted(nstokes, inverse_covariance_matrix[ell_value], 'L');
    }

    for (ell_value=0;ell_value<size_ell_max; ell_value++){
        printf("\n #### ell_Inverted= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
        printf("\n");
    }
    }
    // double normalization_factor = 0;
    // for (ell_value=0; ell_value<lmax; ell_value++){
    //     normalization_factor += (2*ell_value + 1)/(4*M_PI);
    // }

    // printf("Value of normalization factor %f \n", normalization_factor);
    // normalization_factor = normalization_factor;
    // char cholesky_part = 'L';
    // int info;

    // for(ell_value=0; ell_value<lmax; ell_value++){
    //     for (index=0; index < order_matrix; index++){
    //         for (index_2=0; index_2<order_matrix; index_2++){
    //             inverse_covariance_matrix[ell_value][index*order_matrix + index_2] *= 1/normalization_factor;
    //         }
    //     }
    //     dpotrf_(&cholesky_part, &nstokes, inverse_covariance_matrix[ell_value], &nstokes, &info);

    //     inverse_covariance_matrix[ell_value][nstokes] = 0;
    //     if (nstokes == 3){
    //         inverse_covariance_matrix[ell_value][6] = 0;
    //         inverse_covariance_matrix[ell_value][7] = 0;
    //     }
    // }

    printf("%d ----!!!! ***local_map_output*** nstokes %d ; map_size %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    double *local_map_output = (double *)malloc(nstokes*S2HAT_params.Local_param_s2hat.map_size*sizeof(double));

    s2hat_dcomplex *local_alm_s2hat_inverted;

    local_alm_s2hat_inverted = (s2hat_dcomplex *)calloc(3*S2HAT_params.size_alm,sizeof(s2hat_dcomplex));

    printf("%d --- applying inverse of covariance matrix !!! \n", rank); fflush(stdout);
    // apply_inv_block_diag_covariance_matrix_to_alm(local_alm_s2hat, local_alm_s2hat_inverted, inverse_covariance_matrix, &S2HAT_params);

    printf("%d --- harmonic2map !!! --- %d %d \n", rank, nstokes, S2HAT_params.Local_param_s2hat.map_size); fflush(stdout);
    
    // apply_alm2pix(local_map_output, local_alm_s2hat_inverted, &S2HAT_params);
    apply_alm2pix(local_map_output, local_alm_s2hat, &S2HAT_params);

    printf("%d --- end of harmonic2map ! \n", rank); fflush(stdout);



    double *CMB_map_output_copy = (double *)calloc(nstokes*npix,sizeof(double));
    gather_map(local_map_output, CMB_map_output_copy, nstokes, &S2HAT_params);
    memcpy(CMB_map_output, CMB_map_output_copy, nstokes*npix*sizeof(double));
    free(CMB_map_output_copy);
    
    printf("%d --- Reading first elems of local map \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", CMB_map_output[i], CMB_map_ring[i]);
    printf("\n");
    printf("%d --- Reading first elems of local map + npix \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", CMB_map_output[i+npix], CMB_map_ring[i+npix]);
    printf("\n");
    printf("%d --- Reading first elems of local map + 2*npix \n", rank);
    for (i=0; i<max_size_pixels; i++)
        printf("- %f %f -", CMB_map_output[i+2*npix], CMB_map_ring[i+2*npix]);
    printf("\n");


    double average_relative_error_T_0_1 = 0;
    double average_relative_error_Q_0_1 = 0;
    double average_relative_error_U_0_1 = 0;
    double max_average_relative_error_T_0_1 = 0;
    double max_average_relative_error_Q_0_1 = 0;
    double max_average_relative_error_U_0_1 = 0;

    int size_total_probed = npix;
    // size_total_probed = 10;

    int init_index = 0;
    int number_0 = 0;

    index = init_index;
    max_average_relative_error_T_0_1 = fabs((CMB_map_ring[index]-CMB_map_output[index])/CMB_map_ring[index]);
    max_average_relative_error_Q_0_1 = fabs((CMB_map_ring[index+npix]-CMB_map_output[index+npix])/CMB_map_ring[index+npix]);
    max_average_relative_error_U_0_1 = fabs((CMB_map_ring[index+2*npix]-CMB_map_output[index+2*npix])/CMB_map_ring[index+2*npix]);
    for(index=1;index<size_total_probed+init_index;index++){
        average_relative_error_T_0_1 += fabs((CMB_map_ring[index]-CMB_map_output[index])/CMB_map_ring[index]);
        if (max_average_relative_error_T_0_1 < fabs((CMB_map_ring[index]-CMB_map_output[index])/CMB_map_ring[index]))
            max_average_relative_error_T_0_1 = fabs((CMB_map_ring[index]-CMB_map_output[index])/CMB_map_ring[index]);

        if(CMB_map_ring[index] == 0.000000){
            number_0++;
        }
        if(nstokes>1){
            average_relative_error_Q_0_1 += fabs((CMB_map_ring[index+npix]-CMB_map_output[index+npix])/CMB_map_ring[index+npix]);
            if (max_average_relative_error_Q_0_1 < fabs((CMB_map_ring[index+npix]-CMB_map_output[index+npix])/CMB_map_ring[index+npix]))
                max_average_relative_error_Q_0_1 = fabs((CMB_map_ring[index+npix]-CMB_map_output[index+npix])/CMB_map_ring[index+npix]);

            average_relative_error_U_0_1 += fabs((CMB_map_ring[index+2*npix]-CMB_map_output[index+2*npix])/CMB_map_ring[index+2*npix]);
            if (max_average_relative_error_U_0_1 < fabs((CMB_map_ring[index+2*npix]-CMB_map_output[index+2*npix])/CMB_map_ring[index+2*npix]))
                max_average_relative_error_U_0_1 = fabs((CMB_map_ring[index+2*npix]-CMB_map_output[index+2*npix])/CMB_map_ring[index+2*npix]);
        }
    }
    printf("\n");
    printf("Average 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            average_relative_error_T_0_1/(size_total_probed), average_relative_error_Q_0_1/(size_total_probed), average_relative_error_U_0_1/(size_total_probed));
    printf("Maximum 0-1 relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            max_average_relative_error_T_0_1, max_average_relative_error_Q_0_1, max_average_relative_error_U_0_1);

    printf("Number of 0 : %d \n", number_0);
    fflush(stdout);
    

    printf("%d --- Free step 2 ! \n", rank); fflush(stdout);
    if (S2HAT_params.Local_param_s2hat.gangrank >= 0){
        printf("%d --- Free step 1b ! \n", rank); fflush(stdout);
        free(local_alm_s2hat_inverted);
    }
    printf("%d --- Free step 2b ! \n", rank); fflush(stdout);
    // free(local_map_pix);
    free(local_map_output);
    printf("%d --- Free step 3 ! \n", rank); fflush(stdout);
    free_covariance_matrix(inverse_covariance_matrix, lmax);
    printf("%d --- Free step 4 ! \n", rank); fflush(stdout);
    free_s2hat_parameters_struct(&S2HAT_params);

    printf("%d --- Done !!! \n", rank); fflush(stdout);
}
