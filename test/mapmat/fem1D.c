/** @file   fem1D.c
    @brief  test matrix-vector product with specific matrix.
 
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


void usage(){
  fprintf(stderr, "\n*#  Usage: fem1D [options][values] ");
  fprintf(stderr, "\n*#  -N [number of nodes]");
  fprintf(stderr, "\n*#  -f [mat flag]");
  fprintf(stderr, "\n*#  -o [filename]");
}


int main(int argc, char *argv[]){
  int rows, cols, begin, N;
  int Nnz, width;
  int err, i;
  double h;                                   //diamètre élément
  double *K;                                  //raideur
  int *indices;                                 
  Mat At;                                  //matrix
  int *y1_i;                                 //resulting indices vectors
  double *y1_v;                               //resulting values vectors
  double *x1;                                     //parameter vectors
  double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;                                            //timers
  char ch;
  char *filename;
  int rank, size;
  int flag=0;
  int output=0;

  /*mpi init*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /*set options*/
  while((ch = getopt( argc, argv, "N:f:o:" )) != EOF) {
    switch(ch) {
      case 'N':
 	N = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && N <= 0)) {
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
        output=1;
        if (rank==0)
          printf("\nfilename=%s ", filename);
      break;
    }
  }

  /*test procedure :*/
  err=0;
  rows = N; 
  cols = rows;
  Nnz  = 3;
  h=1.0/N;



  partition(&begin, &width, cols, rank, size); 

  
  x1    = (double *) malloc(width* sizeof(double));          /*< allocate vector values  according partitioning ant init to 0.0*/
  for(i=0; i<(width); i++)   
    x1[i]=(1.0);
  

  K  = (double *) malloc(Nnz * width * sizeof(double));           /*< matrice de raideur locale*/ 
  indices  = (int *) malloc(Nnz * width * sizeof(int));       
  for(i=0; i<width; i++){                      //assemblage
    indices[i*Nnz] =(begin+i-1+N)%N;
    indices[i*Nnz+1] = begin +i;
    indices[i*Nnz+2] = (begin+i+1)%N;
    K[i*Nnz] = -1.0 * h;
    K[i*Nnz+1] = 2.0 * h;
    K[i*Nnz+2] = -1.0 * h;
  }
  //for(i=0; i<width*Nnz; i++)                      //assemblage
  //  printf("%d ",indices[i]);
  //printf("\n ");

  MatInit(&At, width, Nnz, indices, K, flag, MPI_COMM_WORLD);             //init distribuated matrix
  //for(i=0; i<width*Nnz; i++)                      //assemblage
  //  printf("%d ",indices[i]);
  //printf("\n ");

  y1_i= (int *) calloc(At.lcount, sizeof(int));                // allocate full vector indices on each proc and init to 0 /
  y1_v  = (double *) calloc(At.lcount, sizeof(double));          // allocate full vector values  on each proc and init to 0.0 /

  MPI_Barrier(MPI_COMM_WORLD);                              // Pre  
  t2 = MPI_Wtime();
  TrMatVecProd_Naive(&At, x1, y1_v, 0);
  t3 = MPI_Wtime();

  for(i=0; i<At.lcount; i++){
    if(y1_v[i] < 0.0){
      err=1;
      if(rank==0)
        printf("\n$ ERROR : rank %d (%d, %10.9f) \n",rank, At.lindices[i], y1_v[i]);
    }
    y1_v[i]=0.0;
  }
  if(err==0){
     printf("$ completed\n");
  }
  else{
   printf("$ failed\n");
  }  

  MatInfo(&At, 0, "FEM");

  MPI_Barrier(MPI_COMM_WORLD);                 
  t6 = MPI_Wtime();
  TrMatVecProd(&At, x1, y1_v, 0);
  t7 = MPI_Wtime();

  for(i=0; i<At.lcount; i++){
    if(  y1_v[i] < -0.0000000000000001){
      err=1;
      printf("\n$ ERROR : rank %d (%d, %15.14f) \n",rank,  y1_i[i], y1_v[i]);
    }
  }
  if(err==0){
     printf("$ completed\n");
  }
  else{
   printf("$ failed\n");
  }  
  
  for(i=0; i<(width); i++){   
      x1[i]=(0.0);
  }

  MPI_Barrier(MPI_COMM_WORLD);                 
  t8 = MPI_Wtime();
  err = MatVecProd(&At, y1_v, x1, 0);
  t9 = MPI_Wtime();
  
  if(rank==0)
    printf("\nInit\t%lf\nTrMatVec_Naive\t%lf\nCommunication\t%lf\nTrMatVec\t%lf\nMatVec\t%lf\n", t1-t0, t3-t2, t5, t7-t6, t9-t8);
  
  
    free(y1_i);
    free(y1_v);
    free(x1); 
    free(K);  
    free(indices);
    MatFree(&At);                                                //free memory/  
  MPI_Finalize();
  
  return err;
}


