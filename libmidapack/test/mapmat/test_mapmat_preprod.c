/** @file   mapmat_preprod.c
    @brief  program to test mapper-matrix/vector preproduct
  
    PreProduct function is used as a pre-step before performing matrix-vector product. It allocates sending and receiving arrays if it is not already done. 
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

    Usage: matmap_save -M [number of rows] -N [number of columns] -Z [number of non-zeros] -B [number of blocks] -f[filename]*/
void usage(){
    fprintf(stderr, "\n*#  Usage: test_mat [options][values] ");
    fprintf(stderr, "\n*#  -M [number of rows]");
    fprintf(stderr, "\n*#  -N [number of columns]");
    fprintf(stderr, "\n*#  -Z [number non-zero values per column]");
    fprintf(stderr, "\n*#  -f [filename]");
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
  int Nnz;
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
  while((ch = getopt( argc, argv, "M:N:Z:f:" )) != EOF) {
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
    }
  }

  /*test starts*/

  err = MapMatGen(&At, MPI_COMM_WORLD, rows, cols, Nnz);
 
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    t0 = MPI_Wtime();
  
  err = MapMatPreProd(&At);
  
  MPI_Barrier(MPI_COMM_WORLD);

  MapMatInfo(&At);

  if(rank==0){
    t1 = MPI_Wtime();
    printf("\nMapMatPreProd execution time =\t%lf sec \n", t1-t0);
  }

  MPI_Finalize();

  return check(err);
}

