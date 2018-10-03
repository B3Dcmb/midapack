// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// Utilitary code to build the map output in gif and binary directly from the input data (point_data and
// pure_signal files).

/** @file   mk_unpointing.c
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


//cluster Adamis:
extern const char *WORKDIR="/data/dauvergn/Test_mapmaking/new_data_set/";
double FKNEE=1.00;
//extern const char *WORKDIR="/data/dauvergn/Test_mapmaking/fred_pack_data/";
//double FKNEE=0.25;
//Hoppper:
//extern const char *WORKDIR="/global/homes/d/dauvergn/data/fred_pack_data/";


int main(int argc, char *argv[])
{

  int i,k;

//  int rank0, idp, m, Alcount, Nb_t_Intervals, t_Interval_length, LambdaBlock, Nb_t_Intervals_loc;
  int nside=256;
  int npix = 12*pow(nside,2);

  double *map;
  map    = (double *) calloc(npix, sizeof(double));


  int Nb_t_Intervals, t_Interval_length;
  int Nnz=1;

//Read data source files:
  int part_id;      // stationnaly period id number
  int *point_data;  // scann strategy input data for the pointing matrix
  double *signal;   // signal input data


  Nb_t_Intervals=8;
  t_Interval_length = pow(2,20);

  point_data = (int *) malloc(Nnz*t_Interval_length * sizeof(int));
  signal     = (double *) malloc(t_Interval_length * sizeof(double));


  for (k=0; k<Nb_t_Intervals; k++) {

  part_id=k;

  ioReadfilePure(t_Interval_length, part_id, point_data, signal);

  printf("#part_id=%d\n", part_id );  //interval id number
  printf("#   indice \t point_data \t pure_signal\n");

//  for(i=0; i<t_Interval_length; i++)
//  printf("%d \t %d \t %lf\n", i, point_data[i], signal[i] );

/*
  int point_data_prev=-1;
  for(i=0; i<t_Interval_length; i++) {
    if (point_data[i] != point_data_prev)
      printf("%10d \t %7d \t %lf\n", i, point_data[i], signal[i] );
  point_data_prev=point_data[i];
  }
*/

  for(i=0; i<t_Interval_length; i++)
    map[point_data[i]]= signal[i];

  }//end of the loop

  free(point_data);
  free(signal);


//write gif output file
  ioWritePixfile(map, nside);

//write binarie output file
  ioWritePixbinfile( npix, map);



 free(map);

  return 0;
}



int ioReadtxtfile( int *size0, int *idp, int *m, int *Alcount, int *Nb_t_Intervals, int *t_Interval_length, int *LambdaBlock, int *Nb_t_Intervals_loc, int map_id)
{

//  int orank0, oidp, om, oAlcount;

  FILE* file;
  char filename [1024];
  sprintf(filename,"mapout%d.txt", map_id);
  printf(" file name: %s\n", filename);

  file = fopen(filename, "r");
  fscanf(file, "%d", size0);
  fscanf(file, "%d", idp );
  fscanf(file, "%d", m );
  fscanf(file, "%d", Alcount );
  fscanf(file, "%d", Nb_t_Intervals );
  fscanf(file, "%d", t_Interval_length );
  fscanf(file, "%d", LambdaBlock );
  fscanf(file, "%d", Nb_t_Intervals_loc );

  fclose(file);


  return 0;
}


int x2map( double *map, int npix, double *x, int *lstid, int xsize)
{

  int i;
  int mask=0;


 if (mask==1) {
  for(i=0; i<npix; i++)
    map[i]= 1.e30;
 }


  for(i=0; i<xsize; i++)
    map[lstid[i]]= x[i];


  return 0;
}


int ioWritePixfile( double *map, int nside)
{

  char fn[256];
  char *bn = "mapoutAll";

  sprintf(fn, "%s.gif", bn);
  printf(" output map file name: %s\n", fn);

  int filenamelen = strlen(fn);

  mapoutput_mp_map2gif_((int*)(&nside), map, &filenamelen, fn);


  return 0;
}


int ioWritePixbinfile( int mapsize, double *map)
{


        char map_vectorFile[256];
        char *map_vectorFileNameBasis = "mapoutAll";

        FILE *fp;

        sprintf(map_vectorFile, "%s.dat", map_vectorFileNameBasis);


        fp=fopen(map_vectorFile, "wb");
        fwrite(map, sizeof(double), mapsize, fp);
        fclose(fp);


        return 0;
}
