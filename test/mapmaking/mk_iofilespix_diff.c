// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012 
// Utilitary code to do the difference between a distributed binaries result and a one processor result
// stored in the Prev directory.

/** @file   mk_iofilespix_diff.c
    @author Frederic Dauvergne
    @date   November 2012 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <midapack.h>

#include <stdbool.h>
#include <errno.h>
#include <unistd.h>
#include <chealpix.h>


extern void __mapoutput_MOD_map2gif(int* nside, double* map, int* filenamelen,char* output_file);


//ioReadbinWritePixfile
int main(int argc, char *argv[])
{

  int i,k;

  int part_id;
  int *point_data;
  double *signal;

  int rank0, idp, m, Alcount, Nb_t_Intervals, t_Interval_length, LambdaBlock, Nb_t_Intervals_loc;

  int nside=256;
  int npix = 12*pow(nside,2);

  double *map, *mapPrev;
  mapPrev    = (double *) calloc(npix, sizeof(double));
  map    = (double *) calloc(npix, sizeof(double));

  int xsize;//= Alcount;
  int map_id=0;  double *x;
  int *lstid;

  int size0;
  int nbproc;

//Read Prev output for 1 proc
  ioReadtxtfilePrev( &size0, &idp, &m, &Alcount, &Nb_t_Intervals, &t_Interval_length, &LambdaBlock, &Nb_t_Intervals_loc, 0);

  nbproc=1;
  printf(" nbproc=%d, m=%d, Alcount=%d, Nb_t_Intervals=%d, t_Interval_length=%d, LambdaBlock=%d\n", nbproc, m, Alcount, 
Nb_t_Intervals, t_Interval_length, LambdaBlock);  

  xsize=Alcount;
  x   = (double *) malloc(xsize*sizeof(double));
  lstid = (int *) malloc(xsize*sizeof(int));

  ioReadbinfilePrev( xsize, map_id, lstid, x);

  x2map( mapPrev, npix, x, lstid, xsize);



//read output from binaries files:
  ioReadtxtfile( &size0, &idp, &m, &Alcount, &Nb_t_Intervals, &t_Interval_length, &LambdaBlock, &Nb_t_Intervals_loc, 0);

  nbproc=size0;
  printf(" nbproc=%d, m=%d, Alcount=%d, Nb_t_Intervals=%d, t_Interval_length=%d, LambdaBlock=%d\n", nbproc, m, Alcount, Nb_t_Intervals, t_Interval_length, LambdaBlock);

  for (k=0; k<nbproc; k++) {
  map_id=k;

  if (k != 0) // already read otherwize
    ioReadtxtfile( &size0, &idp, &m, &Alcount, &Nb_t_Intervals, &t_Interval_length, &LambdaBlock, &Nb_t_Intervals_loc, map_id);


  xsize=Alcount;
  x   = (double *) malloc(xsize*sizeof(double));
  lstid = (int *) malloc(xsize*sizeof(int));

  ioReadbinfile( xsize, map_id, lstid, x);

  x2map( map, npix, x, lstid, xsize);

  free(x);
  free(lstid);
  }//end loop over processors


  //Compute de difference between both maps:
//  mapdiff(mapPrev, map, npix);

   for(i=0; i<npix; i++)
    map[i]= map[i]-mapPrev[i];

 

  ioWritePixfile(map, nside);

  free(map);
  free(mapPrev);

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

  __mapoutput_MOD_map2gif((int*)(&nside), map, &filenamelen, fn);


  return 0;
}


//=====================
//For diff

int ioReadtxtfilePrev( int *size0, int *idp, int *m, int *Alcount, int *Nb_t_Intervals, int *t_Interval_length, int *LambdaBlock, int *Nb_t_Intervals_loc, int map_id)
{

//  int orank0, oidp, om, oAlcount;

  FILE* file;
  char filename [1024];
  sprintf(filename,"Prev/mapout%d.txt", map_id);
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


int ioReadbinfilePrev( int mapsize, int mappart_id, int *lstid, double *map)
{

        int i;

        char lstid_vectorFile[256];
        char x_vectorFile[256];
        char *lstid_vectorFileNameBasis = "Prev/mapout_lstid_";
        char *x_vectorFileNameBasis = "Prev/mapout_";

        FILE *fp;

        sprintf(lstid_vectorFile, "%s%01d.dat", lstid_vectorFileNameBasis, mappart_id);
        sprintf(x_vectorFile, "%s%01d.dat", x_vectorFileNameBasis, mappart_id);

        printf(" Map id file name: %s\n", lstid_vectorFile);
        printf(" Map file name: %s\n", x_vectorFile);


        fp=fopen(lstid_vectorFile, "rb");
        fread(lstid, sizeof(int), mapsize, fp);
        fclose(fp);

        fp=fopen(x_vectorFile, "rb");
        fread(map, sizeof(double), mapsize, fp);
        fclose(fp);


        return 0;
}


