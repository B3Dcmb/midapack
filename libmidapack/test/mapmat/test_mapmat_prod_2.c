/*  @file
    @brief  program to test mapmatrix/vector product
  
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

    Usage: matmap_prod -M [number of rows] -N [number of columns] -Z [number of non-zeros] -f[filename]*/
void usage(){
    fprintf(stderr, "\n*#  Usage: test_mat [options][values] ");
    fprintf(stderr, "\n*#  -M [number of rows]");
    fprintf(stderr, "\n*#  -n [number of local columns]");
    fprintf(stderr, "\n*#  -Z [number non-zero values per column]");
    fprintf(stderr, "\n*#  -f [filename], to specify input file");
    fprintf(stderr, "\n*#  -h help, print usage ");
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
  int rows, cols, begin;
  int indices;
  double values=1.0;
  int Nnz, width;
  int err, i;
  MAPMAT At;
  double t0, t1;
  char ch;
  char *filename;
  char fn[100];
  FILE *in;
  int rank, size;

  /*mpi init*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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
        printf("%s\n", filename);    
      break;
      case 'h':
	usage();
      break;
    }
  }

  /*test starts*/
  int *res_indices;
  res_indices= (int *) calloc(rows, sizeof(int));              /*< allocate full vector indices on each proc and init to 0 */
  double *res_val;
  res_val  = (double *) calloc(rows, sizeof(double));            /*< allocate full vector values  on each proc and init to 0.0 */

  double *vec_val;
  partition(&begin, &width, cols, rank, size); 
  vec_val  = (double *) calloc(width, sizeof(double));           /*< allocate vector values  according partitioning and init to 0.0 */
 

  MapMatCreate(&At, MPI_COMM_WORLD);             //init distribuated matrix
  
  MapMatSetSize(&At, width, Nnz);           //set diffezent sizes  

  MapMatMalloc(&At);                             //table allocation;

  MapMatInfo(&At);

  sprintf(fn,"%s_%d.dat", filename, rank);    
  printf("%s\n", fn);    
  in = fopen(fn, "r");
  i=0;
  while(feof(in)== 0 && i<width){
    fscanf(in, "%ld ", &indices);
    err = MapMatSetIndices(&At, 1, i, 1, &indices);
  //  printf(" %d ", At.indices[i]);
    err = MapMatSetValues(&At, 1, i, 1, &values);
    //printf(" %lf ", At.values[i]);
    i++;  
  }
  fclose(in);
  printf("i=%d\n", i);

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    t0 = MPI_Wtime();

    err = MapMatVecProd(&At, vec_val, res_indices, res_val);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    t1 = MPI_Wtime();
    printf("\n$ MapMatVecProd execution time :\t%lf sec \n", t1-t0);
  }

  /*testing consistancy*/
  for(i=0; i<rows; i++){
    if(res_val[i] != 0.0){
      MPI_Finalize();
      return check(1);
    }
  }
  
  MapMatFree(&At);                                                /*free memory*/  
  free(vec_val);
  free(res_val);
  free(res_indices);
    
  MPI_Finalize();
  
  return check(err);
}


