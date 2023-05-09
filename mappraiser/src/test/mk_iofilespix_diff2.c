// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// Utilitary code to build the difference map output in pixel gif between two results maps stored in binaries
// files.

/** @file   mk_iofilespix_diff2.c
    @author Frederic Dauvergne
    @date   November 2012 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include "midapack.h"

#include <stdbool.h>
#include <errno.h>
#include <unistd.h>
#include <chealpix.h>


extern void mapoutput_mp_map2gif_(int* nside, double* map, int* filenamelen,char* output_file);


//ioReadbinWritePixfile
int main(int argc, char *argv[])
{

  int i,k;

  int nside=256;
  int npix = 12*pow(nside,2);

  double *map01, *map02, *mapDiff;
  map01    = (double *) calloc(npix, sizeof(double));
  map02    = (double *) calloc(npix, sizeof(double));
  mapDiff  = (double *) calloc(npix, sizeof(double));


  ioReadPixbinfile( npix, map01);

  ioReadPixbinfile2( npix, map02);


  //Compute de difference between both maps:
   for(i=0; i<npix; i++)
    mapDiff[i]= map01[i]-map02[i];



  ioWritePixfileDiff(mapDiff, nside);

  free(map01);
  free(map02);
  free(mapDiff);

  return 0;
}


int ioWritePixfileDiff( double *map, int nside)
{

  char fn[256];
  char *bn = "mapoutAllDiff";

  sprintf(fn, "%s.gif", bn);
  printf(" output map file name: %s\n", fn);

  int filenamelen = strlen(fn);

  mapoutput_mp_map2gif_((int*)(&nside), map, &filenamelen, fn);


  return 0;
}


int ioReadPixbinfile( int mapsize, double *map)
{


        char map_vectorFile[256];
        char *map_vectorFileNameBasis = "mapoutAll";

        FILE *fp;

        sprintf(map_vectorFile, "%s.dat", map_vectorFileNameBasis);


        fp=fopen(map_vectorFile, "rb");
        fread(map, sizeof(double), mapsize, fp);
        fclose(fp);


        return 0;
}


int ioReadPixbinfile2( int mapsize, double *map)
{


        char map_vectorFile[256];
        char *map_vectorFileNameBasis = "mapoutAll2";

        FILE *fp;

        sprintf(map_vectorFile, "%s.dat", map_vectorFileNameBasis);


        fp=fopen(map_vectorFile, "rb");
        fread(map, sizeof(double), mapsize, fp);
        fclose(fp);


        return 0;
}
