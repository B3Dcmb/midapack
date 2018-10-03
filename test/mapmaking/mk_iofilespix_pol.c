// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// Utilitary code to build the map output in binary format using distributed binary results
// This contains healpix dependency

/** @file   mk_iofilespix_pol.c
    @author Hamza El Bouhargani
    @date   October 2018 */


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


//ioReadbinWritePixfile
int main(int argc, char *argv[])
{

  int i,k;

  int part_id;
  int *point_data;
  double *signal;

  int rank0, idp, m, Alcount, Nb_t_Intervals, t_Interval_length, LambdaBlock, Nb_t_Intervals_loc;

  int nside=512;
  int npix = 12*pow(nside,2);

  double *mapI;
  mapI    = (double *) calloc(npix, sizeof(double));
  double *mapQ;
  mapQ    = (double *) calloc(npix, sizeof(double));
  double *mapU;
  mapU    = (double *) calloc(npix, sizeof(double));
  int *hits;
  hits = (int *) calloc(npix, sizeof(int));
  double *Cond;
  Cond = (double *) calloc(npix, sizeof(double));


//read output from binaries files:
  int xsize;//= Alcount;
  int map_id=0;
  double *x;
  int *lstid;
  int *lhits;
  double *cond;

  int size0;
  ioReadtxtfile( &size0, &idp, &m, &Alcount, &Nb_t_Intervals, &t_Interval_length, &LambdaBlock, &Nb_t_Intervals_loc, 0);

  int nbproc=size0;
  printf(" nbproc=%d, m=%d, Alcount=%d, Nb_t_Intervals=%d, t_Interval_length=%d, LambdaBlock=%d\n", nbproc, m, Alcount, Nb_t_Intervals, t_Interval_length, LambdaBlock);

  for (k=0; k<nbproc; k++) {
  map_id=k;

  if (k != 0) // already read otherwize
    ioReadtxtfile( &size0, &idp, &m, &Alcount, &Nb_t_Intervals, &t_Interval_length, &LambdaBlock, &Nb_t_Intervals_loc, map_id);


  xsize=Alcount;
  x   = (double *) malloc(xsize*sizeof(double));
  lstid = (int *) malloc(xsize*sizeof(int));
  lhits = (int *) malloc((int)(xsize/3)*sizeof(int));
  cond = (double *) malloc((int)(xsize/3)*sizeof(double));

  ioReadbinfile( xsize, map_id, lstid, lhits, cond, x);

  x2map_pol( mapI, mapQ, mapU, hits, Cond, npix, x, lstid, lhits, cond, xsize);

  free(x);
  free(lstid);
  }//end loop over processors

//binarie one
  ioWritePixbinfile_pol( npix, mapI, mapQ, mapU, hits, Cond);

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


int x2map_pol( double *mapI, double *mapQ, double *mapU, int *hits, double *Cond, int npix, double *x, int *lstid, int *lhits, double *cond, int xsize)
{

  int i;
  int mask=0;


 if (mask==1) {
  for(i=0; i<npix; i++)
    mapI[i]= 1.e30;
    mapQ[i]= 1.e30;
    mapU[i]= 1.e30;
 }


  for(i=0; i<xsize; i++){
    if(i%3 == 0){
      mapI[(int)(lstid[i]/3)]= x[i];
      hits[(int)(lstid[i]/3)]= lhits[(int)(i/3)];
      Cond[(int)(lstid[i]/3)]= cond[(int)(i/3)];
    }
    else if (i%3 == 1)
      mapQ[(int)(lstid[i]/3)]= x[i];
    else
      mapU[(int)(lstid[i]/3)]= x[i];
  }


  return 0;
}


int ioWritePixbinfile_pol( int mapsize, double *mapI, double *mapQ, double *mapU, int *hits, double *Cond)
{


        char mapI_vectorFile[256];
        char *mapI_vectorFileNameBasis = "mapoutAll_I";
        char mapQ_vectorFile[256];
        char *mapQ_vectorFileNameBasis = "mapoutAll_Q";
        char mapU_vectorFile[256];
        char *mapU_vectorFileNameBasis = "mapoutAll_U";
        char hits_vectorFile[256];
        char *hits_vectorFileNameBasis = "mapoutAll_hits";
        char Cond_vectorFile[256];
        char *Cond_vectorFileNameBasis = "mapoutAll_cond";

        FILE *fp;

        sprintf(mapI_vectorFile, "%s.dat", mapI_vectorFileNameBasis);
        sprintf(mapQ_vectorFile, "%s.dat", mapQ_vectorFileNameBasis);
        sprintf(mapU_vectorFile, "%s.dat", mapU_vectorFileNameBasis);
        sprintf(hits_vectorFile, "%s.dat", hits_vectorFileNameBasis);
        sprintf(Cond_vectorFile, "%s.dat", Cond_vectorFileNameBasis);


        fp=fopen(mapI_vectorFile, "wb");
        fwrite(mapI, sizeof(double), mapsize, fp);
        fclose(fp);

        fp=fopen(mapQ_vectorFile, "wb");
        fwrite(mapQ, sizeof(double), mapsize, fp);
        fclose(fp);

        fp=fopen(mapU_vectorFile, "wb");
        fwrite(mapU, sizeof(double), mapsize, fp);
        fclose(fp);

        fp=fopen(hits_vectorFile, "wb");
        fwrite(hits, sizeof(int), mapsize, fp);
        fclose(fp);

        fp=fopen(Cond_vectorFile, "wb");
        fwrite(Cond, sizeof(double), mapsize, fp);
        fclose(fp);


        return 0;
}
