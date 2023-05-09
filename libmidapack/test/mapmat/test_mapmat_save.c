/** @file   test_mapmat_save.c
    @brief  program to save matrix
  
    @author Pierre Cargemel
    @date   December 2011 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "mapmat.h"


extern char *optarg;

/** @brief display common available options

    Usage: test_mapmat_save -M [number of rows] -N [number of columns] -Z [number of non-zeros] -f[filename]*/
void usage(){
    fprintf(stderr, "\n*#  Usage: test_mat [options][values] ");
    fprintf(stderr, "\n*#  -M [value] , to specify global number of rows]");
    fprintf(stderr, "\n*#  -N [value] , to specify global number of columns]");
    fprintf(stderr, "\n*#  -Z [value] , to specify global number non-zero values per column]");
    fprintf(stderr, "\n*#  -f [filename] , output filename");
    fprintf(stderr, "\n*#  -h to print usage\n");
}


/** @brief check and return error code

    */
int check(int err){
    if(err==1){
      printf("..............test failed\n");
    return 1;
    }
    printf("......successfully tested\n");
    return 0;
}


/** @brief display common available options

    */
int main(int argc, char *argv[]){
  int rows, cols;
  int Nnz, Nbls;
  int err;
  MAPMAT At;
  double t0, t1;
  char ch;
  char *filename;
  int rank;

  /*mpi init*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*set options*/
  while((ch = getopt( argc, argv, "M:N:Z:f:h" )) != EOF) {
    switch(ch) {
      case 'N':
 	cols = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && cols <= 0)) {
	  usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
	}
      break;
      case 'M':
	rows = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && rows <= 0)) {
          usage();
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}
      break;
      case 'Z':
        Nnz = strtol(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && Nnz <= 0)) {
	usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
	}
      break;
      case 'f':
        filename = strdup(optarg);
        if (rank==0)
           printf("\nfilename=%s ", filename);
      break;
      case 'h':
        MPI_Finalize();
	if(rank==0)
          usage();
        return 0;
      break;
    }
  }

  /*test starts*/
  MapMatGen(&At, MPI_COMM_WORLD, rows, cols, Nnz);
  MapMatInfo(&At);

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    t0 = MPI_Wtime();
  
  err = MapMatSave(&At, filename);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    t1 = MPI_Wtime();
    printf("\nMapMatSave execution time =\t%lf sec \n", t1-t0);
  }
  
  MapMatFree(&At);
  
  MPI_Finalize();

  return check(err);
}


