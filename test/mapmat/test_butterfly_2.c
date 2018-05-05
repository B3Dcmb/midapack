/** @file   test_butterfly_2.c
    Test butterfly-like communication scheme in the worst case : 
    thats means distributed vectors are fully overlapped.
    Consequently operation consists in summing all vector and get the result on each processor.
    That's equivalent to the so-called MPI allreduce.
    Input : you can specify the size of the vector, N.
    Output : print performances of buuterfly init, and butterfly reduce.
    @author Pierre Cargemel
    @date Septemeber 2012*/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <mpi.h>
#include "../../src/mapmat/butterfly.h"
#include "../../src/mapmat/bitop.h"


extern char *optarg;


int main(int argc, char *argv[]){
  int rank, size, steps;
  char ch;
  int i, k; 
  int *nR, **R, *nS, **S;
  int nRmax, nSmax;
  int *A, nA, *com_A, ncom_A;
  double *vA, *vcom_A;
  FILE *out;
  char fn [100];
  double t0, t1, tmax, tcom, tcommax ;
 
  nA=1;

  MPI_Init(&argc, &argv);                                   /*mpi init*/
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  while((ch = getopt(argc, argv, "N:")) != EOF){          /*set options*/
    switch(ch) {
      case 'N':
 	nA = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nA <= 0)) {
          return 1;
	}
      break;
    }
  }
  //==========================================================
  if(is_pow_2(size)!=0){
    MPI_Finalize();
    return 1;
  }
  
  steps= log_2(size);
  A = (int *) malloc(nA * sizeof(int));
  vA = (double *) malloc(nA * sizeof(double));

  for(i=0; i<nA; i++){
    A[i]=i;
    vA[i]=1.0;
  }

  //==========================================================

  R = (int **) malloc(steps*sizeof(int *));
  nR = (int *) malloc(steps*sizeof(int));
  S = (int **) malloc(steps*sizeof(int *));
  nS = (int *) malloc(steps*sizeof(int));
  
  MPI_Barrier(MPI_COMM_WORLD);
  t0=MPI_Wtime();
  butterfly_init(A, nA, R, nR, S, nS, &com_A, &ncom_A, steps, MPI_COMM_WORLD);
  t1=MPI_Wtime()-t0;

  MPI_Reduce(&t1, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(rank==0)
    printf("pre : %d\t%d\t%lf\n", size, nA, tmax);
  

  //=============================================================================
  vcom_A = (double *) malloc(ncom_A * sizeof(double));
  for(i=0; i<ncom_A; i++){
    vcom_A[i]=0.0;
  } 

  int shared;
 
  shared = m2m(vA, A, nA, vcom_A, com_A, ncom_A);


  nRmax=0;
  nSmax=0;
  for(k=0; k<steps; k++)
    if(nR[k]>nRmax)
      nRmax=nR[k];
  for(k=0; k<steps; k++)
    if(nS[k]>nSmax)
      nSmax=nS[k];



  MPI_Barrier(MPI_COMM_WORLD);
  t0=MPI_Wtime();
  tcom=butterfly_reduce(R, nR, nRmax, S, nS, nSmax, vcom_A, steps, MPI_COMM_WORLD);
  t1=MPI_Wtime()-t0;

  MPI_Reduce(&t1, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(rank==0)
    printf("butter : %d\t%d\t%lf\n", size, nA, tmax);

  MPI_Reduce(&tcom, &tcommax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(rank==0)
    printf("com : %d\t%d\t%lf\n", size, nA, tcommax);
  

  free(nR);
  free(nS);
  free(S);
  free(R);
  free(vcom_A);
  free(com_A);
  free(vA);
  free(A);
  MPI_Finalize();
  return 0;
}

