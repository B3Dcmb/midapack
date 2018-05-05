/** @file   multiI.c
    @brief  test matrix-vector product with specific matrix.

    tested matrix looks like :
    @todo add comments 
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
#include "midapack.h"


extern char *optarg;

/** @brief display common available options

    Usage: multiI -M [number of rows] -N [number of columns] */
void usage(){
    fprintf(stderr, "\n*#  Usage: fem1D [options][values] ");
    fprintf(stderr, "\n*#  -M [number of rows]");
    fprintf(stderr, "\n*#  -N [number of columns]");
    fprintf(stderr, "\n*#  number of rows should divide number of columns");
    fprintf(stderr, "\n*#  -f [mat flag]");
    fprintf(stderr, "\n*#  -o [filename]");
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
  int width;
  int err, i;
  double occ;
  Mat At;                                  //matrix
  int *y1_i;                                 //resulting indices vectors
  double *y1_v;                               //resulting values vectors
  double *x1;                                     //parameter vectors
  double *values;
  int *indices;
  double t0, t1;                                            //timers
  char ch;
  char *filename;
  int rank, size;
  int flag=0;
  int out=0;
  /*mpi init*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /*set options*/
  while((ch = getopt( argc, argv, "N:M:f:o:" )) != EOF) {
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
      case 'f':
 	flag = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && flag < 0)) {
	  usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
	}
      break;
      case 'o':
        filename = strdup(optarg);
        out=1;
        if (rank==0)
           printf("\nfilename=%s ", filename);
      break;
    }
  }

  /*test procedure :*/
  
  if(cols%rows!=0){
    usage();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  occ=(double) cols/rows;
 
  double one =1.0;         

  partition(&begin, &width, cols, rank, size); 

  
  x1    = (double *) malloc(width* sizeof(double));          /*< allocate vector values  according partitioning ant init to 0.0*/

  
  for(i=0; i<(width); i++){   
      x1[i]=(1.0);
  }

  values   = (double *) malloc(width * sizeof(double));      // allocate vector 
  indices    = (int *) malloc(width * sizeof(int));          // allocate vector 

  for(i=0; i<width; i++){   
    indices[i]=(begin+i)%rows;                   
    values[i]=1.0;                   
  }
 
  MatInit(&At, width, 1, indices, values, flag, MPI_COMM_WORLD);             //init distribuated matrix


  y1_i= (int *) calloc(At.lcount, sizeof(int));                /*< allocate full vector indices on each proc and init to 0 */
  y1_v  = (double *) calloc(At.lcount, sizeof(double));          /*< allocate full vector values  on each proc and init to 0.0 */
 
  MPI_Barrier(MPI_COMM_WORLD);                              // first matrix-vector prouct y1= At x1
  if(rank==0)
    t0 = MPI_Wtime();
  err = TrMatVecProd_Naive(&At, x1, y1_v, 0);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    t1 = MPI_Wtime();

  for(i=0; i<At.lcount; i++){
    if( y1_v[i] != occ){
      err=1;
      if(rank==0){
        printf("\n$ ERROR : bad result rank %d (%d, %lf) \n", rank, y1_i[i], y1_v[i]);
      }
    }
    y1_v[i]=0.0;
  }
  if (rank==0)
    printf("\n$ TrMatVecProd_Naive execution time :\t%lf sec \n", t1-t0);

  if (err==0){
    printf("completed");
  }
  else{
    printf("failed");
  }

  
  MatInfo(&At, 0, "MultiI");

  MPI_Barrier(MPI_COMM_WORLD);                              // first matrix-vector prouct y1= At x1
  if(rank==0){t0 = MPI_Wtime();}

  err = TrMatVecProd(&At, x1, y1_v, 0);
  if(rank==0){t1 = MPI_Wtime();}
    
  /*testing consistancy*/
  for(i=0; i<At.lcount; i++){
    if( y1_v[i] != occ){
      err=1;
      if(rank==0){printf("\n$ ERROR : bad result rank %d (%d, %lf) \n",rank, y1_i[i], y1_v[i]);}
    }
  }
    if(rank==0){printf("\n$ TrMatVecProd execution time :\t%lf sec \n", t1-t0);}
  

  if (err==0){
    printf("completed");
  }
  else{
    printf("failed");
  }
 
  MatFree(&At);                                                /*free memory*/  
  free(y1_i);
  free(y1_v);
  free(x1);  
  MPI_Finalize();
  
  return check(err);
}


