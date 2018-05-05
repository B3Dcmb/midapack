/** @file   mapmat_prod.c
    @brief  program to test mapper-matrix/vector product
  
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

    Usage: matmap_save -M [number of rows] -N [number of columns] -Z [number of non-zeros] -f[filename]*/
void usage(){
    
    fprintf("\n*#  Usage: mpirun -n [number of processes] test_mat [option[value]] ");
    fprintf("\n*#  -M [int] set the global number of rows. Rows are uniformly distributed over processes");
    fprintf("\n*#  -N [int] set the global number of columns. ");
    fprintf("\n*#  -Z [int] set the number non-zero values per column (default = 1)");
    fprintf("\n*#  -o [filename] specify an output filename");
    fprintf("\n*#  -t use timers and print execution to benchmark"); 
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
  int Nnz, width;
  int err, i;
  MAPMAT At;                                                //matrix
  int *y1_i, *y2_i, *y3_i;                                 //resulting indices vectors
  double *y1_v, *y2_v, *y3_v;                               //resulting values vectors
  double *x1, *x2, *x3;                                     //parameter vectors
  double lambda;                                            //scalaire
  double t0, t1;                                            //timers
  char ch;
  char *filename;
  int rank, size;

  /*mpi init*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /*set options*/
  while((ch = getopt( argc, argv, "M:N:Z:" )) != EOF) {
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
    }
  }

  /*test procedure :
  generate a random matrix (At) and random vectors (x1, x2 and x3)
  do y1 = At x1,  y2= At x2
  then x3 = lambda x1 + x2
  compute y3 At x3 
  check y3 = lambda y1 + y2*/
  err = MapMatGen(&At, MPI_COMM_WORLD, rows, cols, Nnz);         /*building matrix*/
  err = MapMatInfo(&At);                                       /*print info*/

  y1_i= (int *) calloc(rows, sizeof(int));                /*< allocate full vector indices on each proc and init to 0 */
  y2_i= (int *) calloc(rows, sizeof(int));                /*< allocate full vector indices on each proc and init to 0 */
  y3_i= (int *) calloc(rows, sizeof(int));                /*< allocate full vector indices on each proc and init to 0 */
  y1_v  = (double *) calloc(rows, sizeof(double));          /*< allocate full vector values  on each proc and init to 0.0 */
  y2_v  = (double *) calloc(rows, sizeof(double));          /*< allocate full vector values  on each proc and init to 0.0 */
  y3_v  = (double *) calloc(rows, sizeof(double));          /*< allocate full vector values  on each proc and init to 0.0 */

  lambda=2.0;
  partition(&begin, &width, cols, rank, size); 
  x1  = (double *) malloc(width * sizeof(double));           /*< allocate vector values  according partitioning */
  x2  = (double *) malloc(width * sizeof(double));           /*< allocate vector values  according partitioning */
  x3  = (double *) malloc(width * sizeof(double));           /*< allocate vector values  according partitioning */
  for(i=0; i<width; i++){
    x1[i] = rand();
    x2[i] = rand();
    x3[i] = lambda * x1[i] + x2[i];
  } 

  MPI_Barrier(MPI_COMM_WORLD);                              // PreProd  
  if(rank==0)
    t0 = MPI_Wtime();
  err = MapMatPreProd(&At);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    t1 = MPI_Wtime();
    printf("\n$ MapMatPreProd execution time :\t%lf sec \n", t1-t0);
  }
  MPI_Barrier(MPI_COMM_WORLD);                              // first matrix-vector prouct y1= At x1
  if(rank==0)
    t0 = MPI_Wtime();
  err = MapMatVecProd(&At, y1_v, y1_i, x1);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    t1 = MPI_Wtime();
    printf("\n$ MapMatVecProd execution time :\t%lf sec \n", t1-t0);
  }
  MPI_Barrier(MPI_COMM_WORLD);                              // first matrix-vector prouct y2= At x2
  if(rank==0)
    t0 = MPI_Wtime();
  err = MapMatVecProd(&At, y2_v, y2_i, x2);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    t1 = MPI_Wtime();
    printf("\n$ MapMatVecProd execution time :\t%lf sec \n", t1-t0);
  }
  MPI_Barrier(MPI_COMM_WORLD);                              // first matrix-vector prouct y3= At x3
  if(rank==0)
    t0 = MPI_Wtime();
  err = MapMatVecProd(&At, y3_v, y3_i, x3);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){
    t1 = MPI_Wtime();
    printf("\n$ MapMatVecProd execution time :\t%lf sec \n", t1-t0);
  }


  /*testing consistancy*/
  for(i=0; i<rows; i++){
    if(y1_i[i] != y3_i[i]){
      printf("\n$ ERROR : wrong index at rows %d, (%d %d %d)\n", i, y1_i[i], y2_i[i], y3_i[i]);
      err=1;
    }
    if(y3_v[i] != lambda * y1_v[i] * y2_v[i]){
      printf("\n$ ERROR : wrong value at rows %d, %lf instead %lf \n", i, y3_v[i], lambda * y1_v[i] * y2_v[i]);
      err=1;
    }
  }
  
  MapMatFree(&At);                                                /*free memory*/  
  free(y1_i);
  free(y1_v);
  free(x1);
  free(y2_i);
  free(y2_v);
  free(x2);
  free(y3_i);
  free(y3_v);
  free(x3);
  MPI_Finalize();
  
  return check(err);
}


