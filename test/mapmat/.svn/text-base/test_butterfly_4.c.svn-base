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
  double *vA, *vcom_A, *rcom_A;
  FILE *out;
  char fn [100];
  double t0, t1, tmax, tcom, tcommax ;
 
  MPI_Group world_group;
  MPI_Group node_group;
  MPI_Comm node_comm;
  int core_per_node = 24;  // hopper
  int *node_lead_grank, nnodes, worldsize;

  nA=1;

  MPI_Init(&argc, &argv);                                   /*mpi init*/
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &worldsize);

  while((ch = getopt(argc, argv, "N:S:")) != EOF){          /*set options*/
    switch(ch) {
      case 'N':
 	nA = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nA <= 0)) {
          return 1;
	}
      break;
      case 'S':
 	core_per_node = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nA <= 0)) {
          return 1;
	}
      break;
    }
  }

  MPI_Barrier( MPI_COMM_WORLD);

  nnodes = worldsize/core_per_node;
  node_lead_grank = (int *)calloc( nnodes, sizeof( int));
  for( i = 0; i<nnodes; i++) node_lead_grank[i] = i*core_per_node;

  MPI_Comm_group( MPI_COMM_WORLD, &world_group);
  MPI_Group_incl( world_group, nnodes, node_lead_grank, &node_group);
  MPI_Comm_create( MPI_COMM_WORLD, node_group, &node_comm);

  free( node_lead_grank);

  size =nnodes;

  // fprintf(stdout, " nA = %d [%d < %d]\n", nA, size, worldsize); fflush( stdout);

  //==========================================================
  if(is_pow_2(size)!=0){
    MPI_Finalize();
    return 1;
  }
  
  if( rank == (rank/core_per_node)*core_per_node) {

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
  
  MPI_Barrier(node_comm);
  t0=MPI_Wtime();
  butterfly_init(A, nA, R, nR, S, nS, &com_A, &ncom_A, steps, node_comm);
  t1=MPI_Wtime()-t0;

  MPI_Reduce(&t1, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, node_comm);
  if(rank==0)
    printf("pre : %d\t%d\t%lf\n", size, nA, tmax);
  

  //=============================================================================
  vcom_A = (double *) malloc(ncom_A * sizeof(double));
  for(i=0; i<ncom_A; i++){
    vcom_A[i]=0.0;
  } 

  int shared;
 
  shared = m2m(vA, A, nA, vcom_A, com_A, ncom_A);

  //  printf(" shared = %d\n", shared);


  nRmax=0;
  nSmax=0;
  for(k=0; k<steps; k++)
    if(nR[k]>nRmax)
      nRmax=nR[k];
  for(k=0; k<steps; k++)
    if(nS[k]>nSmax)
      nSmax=nS[k];


  //  fprintf(stdout, "in vA[%d] = %lf\n", rank, vA[100]); fflush( stdout);


  MPI_Barrier(node_comm);
  t0=MPI_Wtime();
  tcom=butterfly_reduce(R, nR, nRmax, S, nS, nSmax, vcom_A, steps, node_comm);
  t1=MPI_Wtime()-t0;

  MPI_Reduce(&t1, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, node_comm);
  if(rank==0)
    printf("butter : %d\t%d\t%lf\n", size, nA, tmax);

  MPI_Reduce(&tcom, &tcommax, 1, MPI_DOUBLE, MPI_MAX, 0, node_comm);
  if(rank==0)
    printf("com : %d\t%d\t%lf\n", size, nA, tcommax);
  
  MPI_Barrier( node_comm);

  if (rank == 0) fprintf(stdout, "STD: result: in <-> out = %lf <-> %lf\n", vA[100], vcom_A[100]); fflush( stdout);

  free(nR);
  free(nS);
  free(S);
  free(R);
  free(com_A);
  free(vcom_A);
  free(vA);
  free(A);

  // true butterfly

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
  
  MPI_Barrier(node_comm);
  t0=MPI_Wtime();
  truebutterfly_init(A, nA, R, nR, S, nS, &com_A, &ncom_A, steps, node_comm);
  t1=MPI_Wtime()-t0;

  MPI_Reduce(&t1, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, node_comm);
  if(rank==0)
    printf("true pre : %d\t%d\t%lf\n", size, nA, tmax);
  

  //=============================================================================
  vcom_A = (double *) malloc(ncom_A * sizeof(double));
  for(i=0; i<ncom_A; i++){
    vcom_A[i]=0.0;
  } 

  shared = m2m(vA, A, nA, vcom_A, com_A, ncom_A);

  //  printf(" shared = %d\n", shared);


  nRmax=0;
  nSmax=0;
  for(k=0; k<steps; k++)
    if(nR[k]>nRmax)
      nRmax=nR[k];
  for(k=0; k<steps; k++)
    if(nS[k]>nSmax)
      nSmax=nS[k];


  //  fprintf(stdout, "in vA[%d] = %lf\n", rank, vA[100]); fflush( stdout);


  MPI_Barrier(node_comm);
  t0=MPI_Wtime();
  tcom=truebutterfly_reduce(R, nR, nRmax, S, nS, nSmax, vcom_A, steps, node_comm);
  t1=MPI_Wtime()-t0;

  MPI_Reduce(&t1, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, node_comm);
  if(rank==0)
    printf("true butter : %d\t%d\t%lf\n", size, nA, tmax);

  MPI_Reduce(&tcom, &tcommax, 1, MPI_DOUBLE, MPI_MAX, 0, node_comm);
  if(rank==0)
    printf("true com : %d\t%d\t%lf\n", size, nA, tcommax);
  
  MPI_Barrier( node_comm);

  if (rank == 0) fprintf(stdout, "TRUE: result: in <-> out = %lf <-> %lf\n", vA[100], vcom_A[100]); fflush( stdout);

  free(nR);
  free(nS);
  free(S);
  free(R);
  free(com_A);
  free(vA);
  free(A);

  // now MPI_ALLREDUCE

  rcom_A = (double *)calloc( nA, sizeof( double));

  for(i=0; i<ncom_A; i++){
    vcom_A[i]=1.0;
  } 

  MPI_Barrier(node_comm);
  t0=MPI_Wtime();
  MPI_Allreduce( vcom_A, rcom_A, nA, MPI_DOUBLE, MPI_SUM, node_comm);
  t1=MPI_Wtime()-t0;

  MPI_Reduce(&t1, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, node_comm);
  if(rank==0)
    printf("allreduce : %d\t%d\t%lf\n", size, nA, tmax);

  MPI_Barrier( node_comm);

  if( rank == 0) fprintf(stdout, "ALL: result: in vA[%d] <-> out vcommA[%d] = %lf <-> %lf\n", rank, rank, vcom_A[100], rcom_A[100]); fflush( stdout);

  free(vcom_A);
  free(rcom_A);

  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
  return 0;
}

