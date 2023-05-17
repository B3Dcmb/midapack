#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
#include <unistd.h>

#include <chealpix.h>
#include "../mapmat/csort.h"
// #include "s2hat.h"
// #include "midapack.h"
#include "s2hat_tools.h"


void make_mask_binary(double* mask, int* mask_binary, int *f_sky, long npix){
  /* Function to transform the mask into binary (composed of 0 and 1 on pixel sky)*/
  long pixel;
  // mask_binary = (int *) calloc( npix, sizeof( int));
  for( pixel=0; pixel<npix; pixel++)
      if( mask[pixel] != 0) {
          mask_binary[pixel] = 1;
          *f_sky = *f_sky + 1;
      }
      else{
        mask_binary[pixel] = 0;
      }
  }


int convert_indices_nest2ring(int *indices_nest, int *indices_ring, long int number_of_indices, int nstokes, int nside){
  // Compute S2HAT ring version of MAPPRAISER nest map map_nest
  // Expect map of nstokes*npix, with (nstokes, npix) for MAPPPRAISER convention (in column-wise order, so [TQUTQUTQUTQU])
  // Expects as well number_of_indices to be the number of pixels * nstokes
  // Return ring map with (npix, nstokes) for S2HAT convention (so [TTTTTTQQQQQUUUUU])
  long int ipix, indice_transformed, npix = 12*nside*nside;
  for( ipix=0; ipix<number_of_indices; ipix++) {

    nest2ring(nside, indices_nest[ipix]/nstokes, &indice_transformed);
    indices_ring[ipix] = indice_transformed + (indices_nest[ipix]%nstokes)*npix;
    // Change indices in nest ordering into indices in ring ordering
  }
  return 0;
}

int convert_indices_ring2nest(int *indices_ring, int *indices_nest, long int number_of_indices, int nstokes, int nside)
{
  // Compute MAPPRAISER nest version of S2HAT ring map map_ring
  // Expect map of nstokes*npix, with (npix, nstokes) for S2HAT convention (in column-wise order, so [TTTTTTQQQQQUUUUU])
  // Expects as well number_of_indices to be the number of pixels * nstokes
  // Return nest map with (nstokes, npix) for MAPPPRAISER convention (so [TQUTQUTQUTQU])

  long int ipix, indice_transformed, npix = 12*nside*nside;
  for( ipix=0; ipix<number_of_indices; ipix++) {
      // ring2nest(nside, indices_ring[ipix]%npix, &indice_transformed);
      // // indices_nest[ipix] = indice_transformed + (indices_ring[ipix]%nstokes)*npix;
      // indices_nest[ipix] = indice_transformed + (indices_ring[ipix]/npix)*nstokes;
      ring2nest(nside, indices_ring[ipix]%npix, &indice_transformed);
      indices_nest[ipix] = indice_transformed*nstokes + (indices_ring[ipix]/npix);
    // Change indices in ring ordering into indices in nest ordering
  }
  return 0;
}


int get_projectors_indices(int *indices_nest, int *ordered_indices_ring, int size_indices, int nstokes, int nside, int *projector_ring2nest, int *projector_nest2ring)
{
  /* Build projectors to convert indices_nest from the nest pixel distribution to ring pixel distribution,
     then reorder them so that the indices in ring order are monotonous 
     The indices_nest are expected to be not have any redundancy */

  int i, j;
  int *indices_ring = (int *)malloc(size_indices*sizeof(int));

  convert_indices_nest2ring(indices_nest, indices_ring, size_indices, nstokes, nside);

  memcpy(ordered_indices_ring, indices_ring, size_indices*sizeof(int));

  ssort(ordered_indices_ring, size_indices, 0); 
  // Argument flag=0 to use quicksort to sort the indices

  for (i=0; i<size_indices; i++) 
  {
    for (j=0; j<size_indices; j++)
    {
      if (ordered_indices_ring[i] == indices_ring[j])
      {
        projector_ring2nest[i] = j;
        projector_nest2ring[j] = i;
        break ;
        // The indices_nest are expected to be not have any redundancy
      }
    }
  }

  free(indices_ring);
  // free(ordered_indices_ring);
  return 0;
}

int project_values_into_different_scheme(double *values_in, int number_values, int *projector_in2out, double *values_out)
{
  int i;
  for (i=0; i<number_values; i++)
    values_out[i] = values_in[projector_in2out[i]];
  
  return 0;
}

void convert_full_map_nest2ring(double *map_nest, double *map_ring, int nside, int nstokes){
  // Compute S2HAT ring version of MAPPRAISER nest map map_nest
  // Expect map of nstokes*npix
  int pixel_index, nstokes_index;
  long ipix;
  long npix = 12*nside*nside;

  for( pixel_index=0; pixel_index<npix; pixel_index++) {
    nest2ring( nside, pixel_index, &ipix);
    for (nstokes_index=0; nstokes_index<nstokes; nstokes_index++){
      map_ring[ipix + nstokes_index*npix] = map_nest[pixel_index*nstokes + nstokes_index];
      // Change map in ring ordering with S2HAT convention [npix, stokes]
      // into map in nest ordering, with MAPPRAISER convention [nstokes, npix] in column-wise ordering
    }
  }
}

void convert_full_map_ring2nest(double *map_ring, double *map_nest, int nside, int nstokes){
  // Compute MAPPRAISER nest version of S2HAT ring map map_ring
  // Expect map of nstokes*npix
  int pixel_index, nstokes_index;
  long ipix;
  long npix = 12*nside*nside;

  for( pixel_index=0; pixel_index<npix; pixel_index++) {
    ring2nest( nside, pixel_index, &ipix);
    for (nstokes_index=0; nstokes_index<nstokes; nstokes_index++){
      map_nest[ipix*nstokes + nstokes_index] = map_ring[pixel_index + nstokes_index*npix];
      // Change map in nest ordering, with MAPPRAISER convention [nstokes, npix]
      // into map in ring ordering with S2HAT convention [npix, stokes] in column-wise ordering
    }
  }
}



int gather_map(double *local_map_pix, double *full_sky_map, int nstokes, S2HAT_parameters *S2HAT_params)
{
    // Gather all S2HAT processes local maps into a full_sky_map

    S2HAT_GLOBAL_parameters *Global_param_s2hat = S2HAT_params->Global_param_s2hat;
    S2HAT_LOCAL_parameters *Local_param_s2hat = S2HAT_params->Local_param_s2hat;
    // int nstokes = 3;
    collect_map(Global_param_s2hat->pixelization_scheme, 1, 0, nstokes, full_sky_map, Local_param_s2hat->first_ring, Local_param_s2hat->last_ring, Local_param_s2hat->map_size,
            local_map_pix, Local_param_s2hat->gangrank, Local_param_s2hat->gangsize, Local_param_s2hat->gangroot, Local_param_s2hat->gangcomm);
    return 0;
}
