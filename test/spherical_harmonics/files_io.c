#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
#include <fitsio.h>
#include <unistd.h>
#include "s2hat.h"
#include "midapack.h"
#include "s2hat_tools.h"

#define EXIT_INFO(Y,Z,args...) { FILE *X=stdout; fprintf( X, "[%s:%d] "Z,__func__, __LINE__, ##args); fflush(X); MPI_Abort( MPI_COMM_WORLD, Y); exit(Y); }
#define INFO(Y,args...)        { FILE *X=stdout; fprintf( X, Y, ##args); fflush(X); }



void read_fits_mask(int nside, double *mask, char *path_mask_file, int col)
/* Function from Xpure, to obtain array from mask path
    MAYBE TO CHANGE (maybe simplify) */
{
  int status = 0, hdutyp, anynul;
  long pixel_index, ipix, npix;
  char ordering[80];
  char comment[81];
  fitsfile *fptr;
  double *tmp;
  char errbuf[31];

  npix = nside*nside*(long)12;

  ffopen( &fptr, path_mask_file, 0, &status);
  if( status) {
    fits_get_errstatus( status, errbuf); 
    printf("\nFITS ERROR : %s\n", errbuf);
    exit(-1);
  }

  fits_movabs_hdu( fptr, 2, &hdutyp, &status);

  fits_read_key( fptr, TSTRING, "ORDERING", ordering, comment, &status);
  if( status) {
    fits_get_errstatus( status, errbuf); 
    printf("\nFITS ERROR : %s (ORDERING) : assumed RING\n", errbuf);
    sprintf( ordering, "RING");
    status=0;
  }

  tmp = (double *) malloc( npix*sizeof(double));

  fits_read_col_dbl( fptr, col, 1, 1, npix, DBL_MAX, tmp, &anynul, &status);
  if( status) {
    fits_get_errstatus( status, errbuf); 
    printf( "\nFITS ERROR: %s\n", errbuf);
  }
  ffclos(fptr, &status);  
  /* revert if NESTED */
  // if( !strcmp( ordering, "NESTED")) {
  //   printf( "NEST -> RING\n");
  //   for( pixel_index=0; pixel_index<npix; pixel_index++) {
  //     nest2ring( nside, pixel_index, &ipix);
  //     mask[ipix] = (double)tmp[pixel_index];
  //   }
  // } else {
  //   for( pixel_index=0; pixel_index<npix; pixel_index++) mask[pixel_index] = (double)tmp[pixel_index];
  // }
  
  for( pixel_index=0; pixel_index<npix; pixel_index++){
    mask[pixel_index] = tmp[pixel_index];
  }
  // printf("Tmp test 3 : %ld %f \t", pixel_index, mask[pixel_index]);
  // fflush(stdout);

  free( tmp);

}

void read_TQU_maps( int nside, double *map, char *infile, int nstokes)
{
  long nele = 12*(long)nside*(long)nside;
  double *mapI = map;
  double *mapQ = map + nele;
  double *mapU = map + nele + nele;

  //READ MAPS
  printf( "           data   : %s ", infile);

  read_fits_mask( nside, mapI, infile, 1);
  printf( "( T ");
  if( nstokes == 3) {
    read_fits_mask( nside, mapQ, infile, 2);
    printf( "Q ");
    read_fits_mask( nside, mapU, infile, 3);
    printf( "U ");
    }
  printf( ")\n");
  fflush( stdout);
}



void make_mask_binary(double* mask, int* mask_binary, int *f_sky, long npix){
  /* Function to transform the mask into binary (composed of 0 and 1 on pixel sky)*/
  long pixel;
  // mask_binary = (int *) calloc( npix, sizeof( int));
  // printf("Tmp test %ld \n", npix);
  for( pixel=0; pixel<npix; pixel++)
      if( mask[pixel] != 0) {
          mask_binary[pixel] = 1;
          // printf("Tmp test %ld \t", pixel,  );
          *f_sky = *f_sky + 1;
      }
  }




void read_fits_cells(int lmax, int number_correl, double *c_ell_array, char *path_list_file, int col)
/* Obtain c_ell array from c_ell path */
{
  int status = 0, hdutyp, anynul;
  long size_c_ell;
  char ordering[80];
  char comment[81];
  fitsfile *fptr;
  double *tmp;
  char errbuf[31];

  size_c_ell = lmax * number_correl; // The size of c_ell array correspond to ell_max * the number of auto-correl and crosscorrelations included [TT, TE, TB, etc.]


  ffopen( &fptr, path_list_file, 0, &status);
  if( status) {
    fits_get_errstatus( status, errbuf); 
    printf("\nFITS ERROR : %s\n", errbuf);
    exit(-1);
  }

  /* move HDU */
  fits_movabs_hdu( fptr, 2, &hdutyp, &status);
  if( status != 0) {
    fits_get_errstatus( status, errbuf);
    EXIT_INFO( status, "%s\n", errbuf);
  }

  /* read size of column */
  // fits_read_key( fptr, TLONG, "NAXIS2", &col, comment, &status);
  // if( status) {
  //   fits_get_errstatus( status, errbuf); 
  //   EXIT_INFO( status, "%s\n", errbuf);
  // }

  fits_read_col_dbl( fptr, col, 1, 1, size_c_ell, DBL_MAX, c_ell_array, &anynul, &status);
  if( status) {
    fits_get_errstatus( status, errbuf); 
    printf( "\nFITS ERROR: %s\n", errbuf);
  }
  ffclos(fptr, &status);
  if( status != 0) {
    fits_get_errstatus( status, errbuf);
    EXIT_INFO( status, "%s\n", errbuf);
  }
}

