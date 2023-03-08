
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
#include <chealpix.h>
#include "midapack.h"
#include "../../include/spherical_harmonics/s2hat_tools.h"


// S2HAT_GLOBAL_parameters;
// /* Get global s2hat structures which must be distributed to all processors*/
// int get_main_s2hat_global_parameters(int nside, int *maskfile_binary, s2hat_pixeltype *pixelization_scheme, s2hat_scandef *scan_sky_structure_pixel, s2hat_pixparameters *pixpar);

// /* Create wrapper structure s2hat of local parameters of s2hat, which will differ for all processors */
// int init_s2hat_global_parameters(Files_path_WIENER_FILTER Files_WF_struct, int *mask_binary, int lmax, S2HAT_GLOBAL_parameters *Global_param_s2hat);

int main_S2HAT_GLOBAL_parameters(){
    int *mask_binary;

    int f_sky, npix;
    int i_index;

    char *path_mask = "/global/cscratch1/sd/mag/Masks_files/SO_wH.fits";
    int nside = 512;
    s2hat_pixeltype pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel;
    s2hat_pixparameters pixpar;

    get_main_s2hat_global_parameters(nside, maskfile_binary, &pixelization_scheme, &scan_sky_structure_pixel, &pixpar);

    printf("Test get_main_s2hat_global_parameters :\n");
    printf(" pixelization_scheme :  type %d, npixsall %d, nphmx %d, nringsall %d \n", pixelization_scheme.type, pixelization_scheme.npixsall, pixelization_scheme.nphmx, pixelization_scheme.nringsall);
    // nph %d, fpix %d, kphi %f, qwqht %f,   :\n", );
    printf(" scan_sky_structure_pixel :  npixsobs %d, nringsobs %d \n", scan_sky_structure_pixel.npixsobs, scan_sky_structure_pixel.nringsobs);
    printf(" pixpar :  par1 %d   par2 %d \n", pixpar.par1, pixpar.par2);
    fflush(stdout);

    destroy_pixelization(pixelization_scheme);
    destroy_scan(scan_sky_structure_pixel);

    int lmax = 1535;

    Files_path_WIENER_FILTER *Files_path_WF_struct;
    Files_path_WF_struct = malloc( 1 * sizeof(Files_path_WIENER_FILTER));
    char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits"; //// TO PUT !!!!
    int number_correlations = 6; // TO MODIFY LATER !!!!!!!
    init_files_struct_WF(Files_path_WF_struct, nside, lmax, c_ell_path, number_correlations);
    
    S2HAT_GLOBAL_parameters *Global_param_s2hat;
    Global_param_s2hat = malloc( 1 * sizeof(S2HAT_GLOBAL_parameters));
    printf("Initializing Global_param_s2hat \n");
    init_s2hat_global_parameters(Files_WF_struct, mask_binary, nside, lmax, Global_param_s2hat);

    s2hat_pixeltype pixelization_scheme_2 = Global_param_s2hat->pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel_2 = Global_param_s2hat->scan_sky_structure_pixel;
    s2hat_pixparameters pixpar_2 = Global_param_s2hat->pixpar;
    printf("###### Test init_s2hat_global_parameters :\n");
    printf(" pixelization_scheme :  type %d, npixsall %d, nphmx %d, nringsall %d \n", pixelization_scheme_2.type, pixelization_scheme_2.npixsall, pixelization_scheme_2.nphmx, pixelization_scheme_2.nringsall);
    // nph %d, fpix %d, kphi %f, qwqht %f,   :\n", );
    printf(" scan_sky_structure_pixel :  npixsobs %d, nringsobs %d \n", scan_sky_structure_pixel_2.npixsobs, scan_sky_structure_pixel_2.nringsobs);
    printf(" pixpar :  par1 %d   par2 %d \n", pixpar_2.par1, pixpar_2.par2);
    fflush(stdout);
    free_s2hat_GLOBAL_parameters_struct(Global_param_s2hat);
    return 0;
}




// void make_mask_binary(double* mask, int* mask_binary, int *f_sky, int npix);
int main_test_make_binary_mask(){
    double *test_mask;
    int *test_mask_binary;
    int npix;
    int i_index;
    int f_sky = 0;

    npix = 30;
    npix = 15;
    test_mask = calloc( npix, sizeof(double));
    for(i_index=0;i_index<npix; i_index++){
        test_mask[i_index] = i_index;
        // test_mask[i_index] = 1;
    }
    test_mask[3] = 0;
    test_mask[10] = 0;
    // test_mask[20] = 0;
    // test_mask[21] = 0;

    test_mask_binary = malloc( npix*sizeof(int));
    make_mask_binary(test_mask, test_mask_binary, &f_sky, npix);

    printf("Binary mask with f_sky= %d \n", f_sky);
    for(i_index=0;i_index<npix; i_index++){
        printf(" %d %f %d \t", i_index, test_mask[i_index], test_mask_binary[i_index]);
    }
    printf("\n");
    fflush(stdout);
    free(test_mask);
    free(test_mask_binary);
    return 0;
}






// S2HAT_LOCAL_parameters;

// /* Create wrapper structure of s2hat of local parameters of s2hat, which will differ for all processors */
// int init_s2hat_local_parameters_struct(S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters *Local_param_s2hat, int *mvals, int gangrank, int gangsize, int gangroot, MPI_Comm gangcomm);


// /* Use s2hat routines to broadcast s2hat global structures */
// void mpi_broadcast_s2hat_global_struc(S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat, int gangroot);

int main_S2HAT_LOCAL_parameters(int argc, char** argv){
// int main(int argc, char** argv){
    
    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/SO_wH.fits";
    char *path_mask = "/global/cscratch1/sd/mag/Masks_files/No_Mask.fits";
    int *mask_binary;

    int f_sky, npix;
    int i_index;

    int nside = 512;
    // int lmax = 1535;
    int lmax = 1024;
    
    S2HAT_GLOBAL_parameters *Global_param_s2hat;
    S2HAT_LOCAL_parameters *Local_param_s2hat;
    
    int rank, nprocs;
    int root, gangroot;
    MPI_Comm gangcomm;
    MPI_Comm root_comm;

    root=0;
    gangroot=0;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

    printf("Initializing MPI %d %d \n", rank, nprocs);
    fflush(stdout);

    int gangrank = rank;
    int gangsize = nprocs;
    gangcomm = MPI_COMM_WORLD;

    
    Global_param_s2hat = malloc( 1 * sizeof(S2HAT_GLOBAL_parameters));
    printf("Initializing Global_param_s2hat \n");
    fflush(stdout);

    if (rank == root){
        init_s2hat_global_parameters(path_mask, nside, lmax, Global_param_s2hat, true);
    }
    
    Local_param_s2hat = malloc( 1 * sizeof(S2HAT_LOCAL_parameters));
    printf("Initializing Local_param_s2hat - rank %d \n", gangrank);
    fflush(stdout);
    init_MPI_struct_s2hat_local_parameters(Local_param_s2hat, gangrank, gangsize, gangroot, gangcomm);
    printf(" MPI  : gangrank %d, gangsize %d, gangroot %d \n", gangrank, gangsize, gangroot);
    printf(" %d - MPI struct in s2hat_local : gangrank %d, gangsize %d, gangroot %d \n", gangrank, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot);
    fflush(stdout);
    
    mpi_broadcast_s2hat_global_struc(Global_param_s2hat, *Local_param_s2hat);
    printf(" %d - S2HAT structures broadcasted ! \n", gangrank);
    fflush(stdout);

    s2hat_pixeltype pixelization_scheme_2 = Global_param_s2hat->pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel_2 = Global_param_s2hat->scan_sky_structure_pixel;
    s2hat_pixparameters pixpar_2 = Global_param_s2hat->pixpar;
    printf("###### Test init_s2hat_global_parameters - %d \n", gangrank);
    printf(" %d - pixelization_scheme :  type %d, npixsall %d, nphmx %d, nringsall %d \n", gangrank, pixelization_scheme_2.type, pixelization_scheme_2.npixsall, pixelization_scheme_2.nphmx, pixelization_scheme_2.nringsall);
    // nph %d, fpix %d, kphi %f, qwqht %f,   :\n", );
    printf(" %d - scan_sky_structure_pixel :  npixsobs %d, nringsobs %d \n", gangrank, scan_sky_structure_pixel_2.npixsobs, scan_sky_structure_pixel_2.nringsobs);
    printf(" %d - pixpar :  par1 %d   par2 %d \n", gangrank, pixpar_2.par1, pixpar_2.par2);
    fflush(stdout);

    init_s2hat_local_parameters_struct(*Global_param_s2hat, Local_param_s2hat);

    printf("###### Test init_s2hat_local_parameters :\n");
    printf(" MPI  : gangrank %d, gangsize %d, gangroot %d \n", gangrank, gangsize, gangroot);
    printf(" %d - MPI struct in s2hat_local : gangrank %d, gangsize %d, gangroot %d \n", gangrank, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot);
    printf(" %d - pixelization_scheme :  plms %d, nmvals %d, first_ring %d, last_ring %d, map_size %d, nplm %ld, mvals[0] %d \n", gangrank,
            Local_param_s2hat->plms, Local_param_s2hat->nmvals, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, Local_param_s2hat->nplm,
            Local_param_s2hat->mvals[0]);
    fflush(stdout);

    free_s2hat_GLOBAL_parameters_struct(Global_param_s2hat);
    printf("Test free 1 - %d \n", gangrank);
    fflush(stdout);
    free_s2hat_LOCAL_parameters_struct(Local_param_s2hat);
    printf("Test free 2 - %d \n", gangrank);
    fflush(stdout);

    MPI_Finalize();
    printf("Test finish ! - %d \n", gangrank);
    fflush(stdout);

    return 0;
}


int main_S2HAT_LOCAL_parameters_v2(int argc, char** argv){        
    char *path_mask = "/global/cscratch1/sd/mag/Masks_files/SO_wH.fits";
    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/No_Mask.fits";
    int *mask_binary;

    int f_sky, npix;
    int i_index;

    int nside = 512;
    // int lmax = 1535;
    int lmax = 1024;
    
    S2HAT_GLOBAL_parameters *Global_param_s2hat;
    S2HAT_LOCAL_parameters *Local_param_s2hat;
    
    int rank, nprocs;
    int root, gangroot;
    MPI_Comm gangcomm;
    MPI_Comm root_comm;

    root=0;
    gangroot=0;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

    printf("Initializing MPI %d %d \n", rank, nprocs);
    fflush(stdout);

    int gangrank = rank;
    int gangsize = nprocs;
    gangcomm = MPI_COMM_WORLD;

    
    Global_param_s2hat = malloc( 1 * sizeof(S2HAT_GLOBAL_parameters));
    printf("Initializing Global_param_s2hat \n");
    fflush(stdout);

    init_s2hat_global_parameters(path_mask, nside, lmax, Global_param_s2hat, true);
    
    
    Local_param_s2hat = malloc( 1 * sizeof(S2HAT_LOCAL_parameters));
    printf("Initializing Local_param_s2hat - rank %d \n", gangrank);
    fflush(stdout);
    init_MPI_struct_s2hat_local_parameters(Local_param_s2hat, gangrank, gangsize, gangroot, gangcomm);
    printf(" MPI  : gangrank %d, gangsize %d, gangroot %d \n", gangrank, gangsize, gangroot);
    printf(" %d - MPI struct in s2hat_local : gangrank %d, gangsize %d, gangroot %d \n", gangrank, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot);
    fflush(stdout);
    
    // mpi_broadcast_s2hat_global_struc(Global_param_s2hat, *Local_param_s2hat);
    // printf(" %d - S2HAT structures broadcasted ! \n", gangrank);
    // fflush(stdout);

    s2hat_pixeltype pixelization_scheme_2 = Global_param_s2hat->pixelization_scheme;
    s2hat_scandef scan_sky_structure_pixel_2 = Global_param_s2hat->scan_sky_structure_pixel;
    s2hat_pixparameters pixpar_2 = Global_param_s2hat->pixpar;
    printf("###### Test init_s2hat_global_parameters - %d \n", gangrank);
    printf(" %d - pixelization_scheme :  type %d, npixsall %d, nphmx %d, nringsall %d \n", gangrank, pixelization_scheme_2.type, pixelization_scheme_2.npixsall, pixelization_scheme_2.nphmx, pixelization_scheme_2.nringsall);
    // nph %d, fpix %d, kphi %f, qwqht %f,   :\n", );
    printf(" %d - scan_sky_structure_pixel :  npixsobs %d, nringsobs %d \n", gangrank, scan_sky_structure_pixel_2.npixsobs, scan_sky_structure_pixel_2.nringsobs);
    printf(" %d - pixpar :  par1 %d   par2 %d \n", gangrank, pixpar_2.par1, pixpar_2.par2);
    fflush(stdout);

    init_s2hat_local_parameters_struct(*Global_param_s2hat, Local_param_s2hat);

    printf("###### Test init_s2hat_local_parameters :\n");
    printf(" MPI  : gangrank %d, gangsize %d, gangroot %d \n", gangrank, gangsize, gangroot);
    printf(" %d - MPI struct in s2hat_local : gangrank %d, gangsize %d, gangroot %d \n", gangrank, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot);
    printf(" %d - pixelization_scheme :  plms %d, nmvals %d, first_ring %d, last_ring %d, map_size %d, nplm %ld, mvals[0] %d \n", gangrank,
            Local_param_s2hat->plms, Local_param_s2hat->nmvals, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size, Local_param_s2hat->nplm,
            Local_param_s2hat->mvals[10]);
    fflush(stdout);

    free_s2hat_GLOBAL_parameters_struct(Global_param_s2hat);
    printf("Test free 1 - %d \n", gangrank);
    fflush(stdout);
    free_s2hat_LOCAL_parameters_struct(Local_param_s2hat);
    printf("Test free 2 - %d \n", gangrank);
    fflush(stdout);

    MPI_Finalize();
    printf("Test finish ! - %d \n", gangrank);
    fflush(stdout);

    return 0;
}


// /* Transform alm coefficients local_alm into a pixel map local_map_pix */
// int apply_alm2pix(s2hat_dcomplex *local_alm, double *local_map_pix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);

// /* Transform local pixel map into local alm coefficients */
// int apply_pix2alm(double *local_map_pix, s2hat_dcomplex *local_alm, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);

int main_alm_pix_tools(int argc, char** argv){        
// int main(int argc, char** argv){
    
    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/SO_wH.fits";
    char *path_mask = "/global/cscratch1/sd/mag/Masks_files/No_Mask.fits";
    char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_multiplesim_white_noise_1.fits";
    int *mask_binary;
    double *CMB_map;

    int f_sky, npix;
    int i_index, ncomp;
    int index;

    int nside = 512;
    // int lmax = 1535;
    // int lmax = 1024;
    int lmax = 3*nside-1 ;//1500;//025;
    int nstokes = 3;
    
    S2HAT_GLOBAL_parameters *Global_param_s2hat;
    S2HAT_LOCAL_parameters *Local_param_s2hat;
    
    int rank, nprocs;
    int root, gangroot;
    MPI_Comm gangcomm;
    MPI_Comm root_comm;

    root=0;
    gangroot=0;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

    printf("Initializing MPI %d %d \n", rank, nprocs);
    fflush(stdout);

    int gangrank = rank;
    int gangsize = nprocs;
    gangcomm = MPI_COMM_WORLD;



    Global_param_s2hat = (S2HAT_GLOBAL_parameters *) malloc( 1 * sizeof(S2HAT_GLOBAL_parameters));
    printf("Initializing Global_param_s2hat \n");
    init_s2hat_global_parameters(path_mask, nside, lmax, Global_param_s2hat, true);

    Local_param_s2hat = (S2HAT_LOCAL_parameters *) malloc( 1 * sizeof(S2HAT_LOCAL_parameters));
    printf("Initializing Local_param_s2hat - rank %d \n", gangrank);
    init_MPI_struct_s2hat_local_parameters(Local_param_s2hat, gangrank, gangsize, gangroot, gangcomm);
    init_s2hat_local_parameters_struct(*Global_param_s2hat, Local_param_s2hat);

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
    printf("Minute test1 - rank %d - %f \n", gangrank, CMB_map[0]);
    // free(CMB_map);
    // printf("Minute test2 - rank %d - %f \n", gangrank, new_CMB_map[0]);
    // printf("Changing map - rank %d \n", gangrank);
    // fflush(stdout);

    double *local_map_pix;
    local_map_pix = (double *) calloc( nstokes*Local_param_s2hat->map_size, sizeof(double));
    distribute_map(Global_param_s2hat->pixelization_scheme, 1, 0, nstokes, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size,
		      local_map_pix, CMB_map, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);
    
    printf("Tiny test - rank %d - %f \n", gangrank, CMB_map[0]);
    
    s2hat_dcomplex *local_alm;
    local_alm = (s2hat_dcomplex *) malloc( nstokes*(Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax*sizeof(s2hat_dcomplex));
    
    
    printf("Little test before pix2alm - rank %d - %f %f \n", gangrank, local_map_pix[0], CMB_map[0]);
    fflush(stdout);

    int max_size_test = 10;

    printf("Map_pix after distribute_map - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    for (index=1;index<max_size_test;index++){
        printf("- %f %f -", local_map_pix[index], CMB_map[index]);
    }
    printf(" \n");

    apply_pix2alm(local_map_pix, local_alm, *Global_param_s2hat, *Local_param_s2hat); //////////////////////////////////////////////////////////////////////
    
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
    apply_alm2pix(local_alm, local_map_pix, *Global_param_s2hat, *Local_param_s2hat); //////////////////////////////////////////////////////////////////////

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
    gather_map(local_map_pix, full_sky_map_2, *Global_param_s2hat, *Local_param_s2hat);
    double *local_map_pix_2;
    local_map_pix_2 = (double *) calloc( nstokes*Local_param_s2hat->map_size, sizeof(double));
    distribute_map(Global_param_s2hat->pixelization_scheme, 1, 0, nstokes, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size,
		      local_map_pix_2, CMB_map, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);

    apply_pix2alm(local_map_pix_2, local_alm, *Global_param_s2hat, *Local_param_s2hat); //////////////////////////////////////////////////////////////////////

    // printf("Map_pix after 2nd pix2alm - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_map_pix[index], CMB_map[index]);
    // }
    // printf(" \n");
    
    // printf("Local_alm after 2nd pix2alm - rank %d - re %f im %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // }
    // printf(" \n");
    // fflush(stdout);
    
    apply_alm2pix(local_alm, local_map_pix, *Global_param_s2hat, *Local_param_s2hat); //////////////////////////////////////////////////////////////////////
    
    // printf("Map_pix after 2nd alm2pix - rank %d - %f %f -", gangrank, local_map_pix[0], CMB_map[0]);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_map_pix[index], CMB_map[index]);
    // }
    // printf(" \n");
    
    // printf("Local_alm after 2nd alm2pix - rank %d - re %f im %f -", gangrank, local_alm[0].re, local_alm[0].im);
    // for (index=1;index<max_size_test;index++){
    //     printf("- %f %f -", local_alm[index].re, local_alm[index].im);
    // }
    // printf(" \n");

    double *full_sky_map;
    full_sky_map = (double *) malloc(nstokes*npix*sizeof(double));
    gather_map(local_map_pix, full_sky_map, *Global_param_s2hat, *Local_param_s2hat);
    

    if (Local_param_s2hat->gangrank == 0){
    printf("//////////////////////// Comparison new map vs old \n");
    fflush(stdout);

    printf("From pixel 0 to %d - %f %f -", max_size_test, full_sky_map[0], CMB_map[0]);
    for (index=0;index<max_size_test;index++){
        printf("- %f %f -", full_sky_map[index], CMB_map[index]);
        }
    printf(" \n");
    fflush(stdout);

    int new_pixel = npix;
    printf("From pixel %d to %d - %f %f -", new_pixel, new_pixel+max_size_test, full_sky_map[new_pixel], CMB_map[new_pixel]);
    for (index=new_pixel;index<new_pixel+max_size_test;index++){
        printf("- %f %f -", full_sky_map[index], CMB_map[index]);
            }
    printf(" \n");
    fflush(stdout);

    new_pixel = 2*npix;
    printf("From pixel %d to %d - %f %f -", new_pixel, new_pixel+max_size_test, full_sky_map[new_pixel], CMB_map[new_pixel]);
    for (index=new_pixel;index<new_pixel+max_size_test;index++){
        printf("- %f %f -", full_sky_map[index], CMB_map[index]);
        }
    printf(" \n");

    double average_relative_error_T = 0;
    double average_relative_error_Q = 0;
    double average_relative_error_U = 0;
    int size_total_probed = npix;
    // size_total_probed = 10;

    int init_index = 0;
    int number_0 = 0;

    for(index=init_index;index<size_total_probed+init_index;index++){
        average_relative_error_T += fabs((full_sky_map[index]-full_sky_map_2[index])/full_sky_map_2[index]);
        // printf("Test %d ---- %f \t", index, fabs((full_sky_map[index]-CMB_map[index])/CMB_map[index]));
        if(full_sky_map[index] == 0.000000){
            number_0++;
        }
        average_relative_error_Q += fabs((full_sky_map[index+npix]-full_sky_map_2[index+npix])/full_sky_map_2[index+npix]);
        // printf("Test2 %d ---- %f \t", index, (full_sky_map[index+npix]-CMB_map[index+npix])/CMB_map[index+npix]);
        average_relative_error_U += fabs((full_sky_map[index+2*npix]-full_sky_map_2[index+2*npix])/full_sky_map_2[index+2*npix]);
        // printf("Test3 %d ---- %f \t", index, (full_sky_map[index+2*npix]-CMB_map[index+2*npix])/CMB_map[index+2*npix]);
    }
    printf("\n");
    printf("Average relative error on all pixels with %d %d : T %f ; Q %f ; U %f \n", nstokes, npix, 
            average_relative_error_T/(size_total_probed), average_relative_error_Q/(size_total_probed), average_relative_error_U/(size_total_probed));
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

    
    printf("Test 9 !\n");
    fflush(stdout);
    free_s2hat_GLOBAL_parameters_struct(Global_param_s2hat);
    printf("Test 10 !\n");
    fflush(stdout);
    free_s2hat_LOCAL_parameters_struct(Local_param_s2hat);

    MPI_Finalize();
    printf("Test finish ! - %d \n", gangrank);
    fflush(stdout);

    return 0;
}




// /* General function to inverse matrix using LAPACK */
// int get_inverse_matrix(int order_matrix, double* matrix_to_be_inverted);

// /* Read c_ell to generate covariance matrix which will be in the form : covariance_matrix_3x3[lmax][9] with 9 being [TT, TE, TB, ET, EE, EB, BT, BE, BB] (with TE=ET, TB=BT and BE=EB) */
// int get_covariance_matrix_3x3(char *c_ell_path, int number_correl, double **covariance_matrix_3x3, S2HAT_GLOBAL_parameters Global_param_s2hat);

// /* Function to obtain inverse of covariance matrix in harmonic domain, from given c_ells */
// int get_inverse_covariance_matrix_3x3(char *c_ell_path, int number_correl, double **inverse_covariance_matrix, S2HAT_GLOBAL_parameters Global_param_s2hat);

// /* Free covariance matrix */
// void free_covariance_matrix(double ** covariance_matrix_3x3, int lmax);


int main_covariance_matrix_tools(int argc, char** argv){
    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/SO_wH.fits";
    char *path_mask = "/global/cscratch1/sd/mag/Masks_files/No_Mask.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits"; //// TO PUT !!!!
    int *mask_binary;

    int f_sky, npix;
    int i_index, ncomp;
    int index, index_2, ell_value;

    int nside = 512;
    // int lmax = 3*nside-1 ;//1500;//025;
    int lmax = 4;
    int nstokes = 3;
    
    S2HAT_GLOBAL_parameters *Global_param_s2hat;
    S2HAT_LOCAL_parameters *Local_param_s2hat;
    
    int rank, nprocs;
    int root, gangroot;
    // MPI_Comm gangcomm;
    // MPI_Comm root_comm;

    root=0;
    gangroot=0;

    // MPI_Init( &argc, &argv);
    // MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    // MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

    // printf("Initializing MPI %d %d \n", rank, nprocs);
    // fflush(stdout);

    // int gangrank = rank;
    // int gangsize = nprocs;
    // gangcomm = MPI_COMM_WORLD;

    Global_param_s2hat = (S2HAT_GLOBAL_parameters *) malloc( 1 * sizeof(S2HAT_GLOBAL_parameters));
    printf("Initializing Global_param_s2hat \n");
    init_s2hat_global_parameters(path_mask, nside, lmax, Global_param_s2hat, true);

    // Local_param_s2hat = (S2HAT_LOCAL_parameters *) malloc( 1 * sizeof(S2HAT_LOCAL_parameters));
    // printf("Initializing Local_param_s2hat - rank %d \n", gangrank);
    // init_MPI_struct_s2hat_local_parameters(Local_param_s2hat, gangrank, gangsize, gangroot, gangcomm);
    // init_s2hat_local_parameters_struct(*Global_param_s2hat, Local_param_s2hat);


    int order_matrix = 3;
    double matrix_to_be_inverted[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    // matrix_to_be_inverted = (double *)calloc(order_matrix*order_matrix, sizeof(double));
    matrix_to_be_inverted ;
    get_inverse_matrix(order_matrix, matrix_to_be_inverted);
    printf("First inversion test, with identity : \n");
    for (index=0; index < order_matrix; index++){
        for (index_2=0; index_2<order_matrix; index_2++){
            printf("%f \t", matrix_to_be_inverted[index*order_matrix + index_2]);
        }
        printf("\n");
    }
    // matrix_to_be_inverted = {0, 0, 2, 2, 0, 0, 0, 2, 0};
    double matrix_to_be_inverted_2[9] = {0, 0, 2, 2, 0, 0, 0, 2, 0};
    get_inverse_matrix(order_matrix, matrix_to_be_inverted_2);
    printf("Second inversion test, with permutation : \n");
    for (index=0; index < order_matrix; index++){
        for (index_2=0; index_2<order_matrix; index_2++){
            printf("%f \t", matrix_to_be_inverted_2[index*order_matrix + index_2]);
        }
        printf("\n");
    }
    double matrix_to_be_inverted_3[9] = {.1, 0, .1, .1, .1, 0, 0, .1, .1};
    get_inverse_matrix(order_matrix, matrix_to_be_inverted_3);
    printf("Third inversion test, with .1 : \n");
    for (index=0; index < order_matrix; index++){
        for (index_2=0; index_2<order_matrix; index_2++){
            printf("%f \t", matrix_to_be_inverted_3[index*order_matrix + index_2]);
        }
        printf("\n");
    }
    double matrix_to_be_inverted_4[9] = {1, 1, .1, 1, 1, .1, 1, 1, 1};
    get_inverse_matrix(order_matrix, matrix_to_be_inverted_4);
    printf("Fourth inversion test, with non-inversible : \n");
    for (index=0; index < order_matrix; index++){
        for (index_2=0; index_2<order_matrix; index_2++){
            printf("%f \t", matrix_to_be_inverted_4[index*order_matrix + index_2]);
        }
        printf("\n");
    }
    // free(matrix_to_be_inverted);

    /////////////////////////////////////////////////////
    int n_correl_to_get = 4; // or 6
    double **covariance_matrix_3x3;

    
    
    covariance_matrix_3x3 = calloc(lmax+1, sizeof(double *));
    for(ell_value=0; ell_value<lmax+1; ell_value++){
        covariance_matrix_3x3[ell_value] = calloc(9,sizeof(double));
    }
    printf("Covariance matrix computed from file : \n");
    get_covariance_matrix_3x3(c_ell_path, n_correl_to_get, covariance_matrix_3x3, *Global_param_s2hat);
    
    order_matrix = 3;
    for (ell_value=0;ell_value<lmax;ell_value++){
        printf("\n #### ell= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", covariance_matrix_3x3[ell_value][index*order_matrix + index_2]);
            }
            printf("\n");
        }
        
    }
    
    fflush(stdout);
    


    printf("Getting inverse of covariance matrix \n");
    fflush(stdout);
    double **inverse_covariance_matrix;
    inverse_covariance_matrix = calloc(lmax+1, sizeof(double *));
    for(ell_value=0; ell_value<lmax+1; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(9,sizeof(double));
    }

    get_inverse_covariance_matrix_3x3(c_ell_path, n_correl_to_get, inverse_covariance_matrix, *Global_param_s2hat);
    printf("Inverse of covariance matrix computed from file : \n");
    order_matrix = 3;
    for (ell_value=0;ell_value<lmax;ell_value++){
        printf("\n #### ell= %d \n", ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.7f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
            printf("\n");
        }
        
    }
    fflush(stdout);
    
    free_covariance_matrix(covariance_matrix_3x3, lmax);
    free_covariance_matrix(inverse_covariance_matrix, lmax);

    free_s2hat_GLOBAL_parameters_struct(Global_param_s2hat);
    // free_s2hat_LOCAL_parameters_struct(Local_param_s2hat);
    // MPI_Finalize();

    return 0;
}



// /* Apply inverse of covariance matrix to local_alm */
// int apply_inv_covariance_matrix_to_alm(s2hat_dcomplex *local_alm, double **inv_covariance_matrix, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);


// /* Transform alm to c_ell coefficients */
// int alm2cls(s2hat_dcomplex *local_alm, double *c_ell_array, S2HAT_GLOBAL_parameters Global_param_s2hat, S2HAT_LOCAL_parameters Local_param_s2hat);


int main_covariance_alms(int argc, char** argv){
// int main(int argc, char** argv){
    // char *path_mask = "/global/cscratch1/sd/mag/Masks_files/SO_wH.fits";
    char *path_mask = "/global/cscratch1/sd/mag/Masks_files/No_Mask.fits";
    char *c_ell_path = "/global/homes/m/mag/midapack/test/spherical_harmonics/test_functions/c_ell_file_lmax_4.fits";
    char *path_CMB_map = "/global/cscratch1/sd/mag/xPure_data/Files_Launch/Map_multiplesim_white_noise_1.fits";
    int *mask_binary;
    
    int f_sky, npix;
    int i_index, ncomp;
    int index, index_2, ell_value;

    int nside = 512;
    // int lmax = 3*nside-1 ;//1500;//025;
    int lmax = 4;
    int nstokes = 3;
    
    S2HAT_GLOBAL_parameters *Global_param_s2hat;
    S2HAT_LOCAL_parameters *Local_param_s2hat;
    
    int rank, nprocs;
    int root, gangroot;
    MPI_Comm gangcomm;
    MPI_Comm root_comm;

    root=0;
    gangroot=0;

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

    printf("Initializing MPI %d %d \n", rank, nprocs);
    fflush(stdout);

    int gangrank = rank;
    int gangsize = nprocs;
    gangcomm = MPI_COMM_WORLD;

    Global_param_s2hat = (S2HAT_GLOBAL_parameters *) malloc( 1 * sizeof(S2HAT_GLOBAL_parameters));
    printf("Initializing Global_param_s2hat \n");
    init_s2hat_global_parameters(path_mask, nside, lmax, Global_param_s2hat, true);

    Local_param_s2hat = (S2HAT_LOCAL_parameters *) malloc( 1 * sizeof(S2HAT_LOCAL_parameters));
    printf("Initializing Local_param_s2hat - rank %d \n", gangrank);
    init_MPI_struct_s2hat_local_parameters(Local_param_s2hat, gangrank, gangsize, gangroot, gangcomm);
    init_s2hat_local_parameters_struct(*Global_param_s2hat, Local_param_s2hat);
    
    /////////////////////////////////////////
    npix = 12*nside*nside;   
    double *CMB_map; 
    CMB_map = (double *) malloc( nstokes*npix*sizeof(double));
    read_TQU_maps(nside, CMB_map, path_CMB_map, nstokes);
    // printf("Reading map - rank %d \n", gangrank);
    // fflush(stdout);

    printf("Minute test1 - rank %d - %f \n", gangrank, CMB_map[0]);
    // free(CMB_map);
    // printf("Minute test2 - rank %d - %f \n", gangrank, new_CMB_map[0]);
    // printf("Changing map - rank %d \n", gangrank);
    // fflush(stdout);

    double *local_map_pix;
    local_map_pix = (double *) calloc( nstokes*Local_param_s2hat->map_size, sizeof(double));
    distribute_map(Global_param_s2hat->pixelization_scheme, 1, 0, nstokes, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size,
		      local_map_pix, CMB_map, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);
    
    s2hat_dcomplex *local_alm;
    local_alm = (s2hat_dcomplex *) malloc( nstokes*(Global_param_s2hat->nlmax+1)*Global_param_s2hat->nmmax*sizeof(s2hat_dcomplex));
    
    apply_pix2alm(local_map_pix, local_alm, *Global_param_s2hat, *Local_param_s2hat); //////////////////////////////////////////////////////////////////////
    
    printf("%d - Local_alm : \n", gangrank);
    int nmvals = Local_param_s2hat->nmvals;
    int index_stokes, m_value;
    
    for (ell_value=0;ell_value<lmax;ell_value++){
        printf("\n %d -alms - #### ell= %d \n", gangrank, ell_value);
        for (m_value=0; m_value < nmvals; m_value++){
            for(index_stokes=0;index_stokes<nstokes;index_stokes++){
                // printf("%.9f %.9f - \t", local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + index_stokes].re, local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + index_stokes].im);
                printf("%.9f %.9f - \t", local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + index_stokes].re, local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + index_stokes].im);
            }
            printf(" ## ");
        }
    }
    printf("\n");
    fflush(stdout);

    ////////////////////////////////////////

    printf("Getting inverse of covariance matrix - %d\n", gangrank);
    fflush(stdout);
    int n_correl_to_get = 4; // or 6

    double **inverse_covariance_matrix;
    // inverse_covariance_matrix = calloc(lmax+1, sizeof(double *));
    inverse_covariance_matrix = calloc(lmax, sizeof(double *));
    // for(ell_value=0; ell_value<lmax+1; ell_value++){
    //     inverse_covariance_matrix[ell_value] = calloc(9,sizeof(double));
    // }
    for(ell_value=0; ell_value<lmax; ell_value++){
        inverse_covariance_matrix[ell_value] = calloc(9,sizeof(double));
    }

    get_inverse_covariance_matrix_3x3(c_ell_path, n_correl_to_get, inverse_covariance_matrix, *Global_param_s2hat);
    printf("Inverse of covariance matrix computed from file : %d \n", gangrank);
    fflush(stdout);

    int order_matrix = 3;

    // double identity[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    // for(ell_value=0;ell_value<lmax;ell_value++){
    //     for (index=0; index<order_matrix ; index++){
    //         for (index_2=0; index_2<order_matrix ; index_2++){
    //             inverse_covariance_matrix[ell_value][index*order_matrix+index_2] = identity[index*order_matrix+index_2];
    //         }
    //     }
    // }

    
    for (ell_value=0;ell_value<lmax;ell_value++){
        printf("\n %d #### -inv_cov- ell= %d \n", gangrank, ell_value);
        for (index=0; index < order_matrix; index++){
            for (index_2=0; index_2<order_matrix; index_2++){
                printf("%.9f \t", inverse_covariance_matrix[ell_value][index*order_matrix + index_2]);
            }
            printf("\n");
        }
    }

    apply_inv_covariance_matrix_to_alm(local_alm, inverse_covariance_matrix, *Global_param_s2hat, *Local_param_s2hat);

    printf("%d - Local_alm after applying inv covariance matrix : \n", gangrank);
    // int nmvals = Local_param_s2hat.nmvals;
    
    
    for (ell_value=0;ell_value<lmax;ell_value++){
        printf("\n %d - new local_alm - #### ell= %d \n", gangrank, ell_value);
        for (m_value=0; m_value < nmvals; m_value++){
            for(index_stokes=0;index_stokes<nstokes;index_stokes++){
                // printf("%.9f %.9f - \t", local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + index_stokes].re, local_alm[m_value*(lmax+1)*nstokes + ell_value*nstokes + index_stokes].im);
                printf("%.9f %.9f - \t", local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + index_stokes].re, local_alm[m_value*(lmax)*nstokes + ell_value*nstokes + index_stokes].im);
            }
            printf(" ## ");
        }
    }
    printf("\n");
    fflush(stdout);
    
    double *c_ell_array;

    printf("Collecting c_ells - %d \n", gangrank);
    fflush(stdout);
    int index_correl, nspec = 6;
    c_ell_array = (double *)calloc( nspec*(lmax), sizeof(double));
    
    alm2cls(local_alm, c_ell_array, nspec, *Global_param_s2hat, *Local_param_s2hat);

    if (gangrank==0){
        for(index_correl=0;index_correl<nspec;index_correl++){
            printf("\n %d - CELLS #### correl = %d \n", gangrank, index_correl);
            for (ell_value=0;ell_value<lmax;ell_value++){
                // printf("\n %d - CELLS #### ell= %d correl = %d \n", gangrank, ell_value, index_correl);
                printf("%.9f \t", c_ell_array[index_correl*(lmax) + ell_value]);
            }
        }
        fflush(stdout);
    }
    printf("\n");
    free(c_ell_array);

    // free_covariance_matrix(covariance_matrix_3x3, lmax);
    free_covariance_matrix(inverse_covariance_matrix, lmax);

    free_s2hat_GLOBAL_parameters_struct(Global_param_s2hat);
    // free_s2hat_LOCAL_parameters_struct(Local_param_s2hat);
    // MPI_Finalize();

    return 0;
}
