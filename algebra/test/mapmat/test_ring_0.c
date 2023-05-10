/** @file   test_butterfly_0.c
    @brief <b> Test butterfly-like communication scheme : butterfly_init .</b>
    @n Usage : mpirun -n [nb procs] test_butterfly_0 -n [elements in the local set] -d [derive]. 
    The number of elements is the size of the local array.
    The option -d enable to specify a derive parameter to generate random ordered sets. Bigger is d, faster element increases.
    @author Pierre Cargemel
    @date April 2012*/
 
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
  int i, k, d;
  int nA; 
  int *A, *nR, **R, *nS, **S;
  int *com_A, ncom_A;
  FILE *out;
  char fn [100];
  double t0, t1;
  int rmax, ravg, rmin;
 
  nA=1;
  d=2;

  MPI_Init(&argc, &argv);                                   /*mpi init*/
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  while((ch = getopt(argc, argv, "n:d:")) != EOF){          /*set options*/
    switch(ch) {
      case 'n':
 	nA = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nA <= 0)) {
          return 1;
	}
      break;
      case 'd':
 	d = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && d <= 0)) {
          return 1;
	}
      break;
    }
  }
  //==========================================================
  if(is_pow_2(size)!=0){
    if(rank==0)
      printf("\nnumber of processes is not a power of 2\n"); 
    MPI_Finalize();
    return 0;
  }
  
  steps= size;
  srand(nA*rank);
  A = (int *) malloc(nA * sizeof(int));
  sprintf(fn,"output_%d.dat", rank);    
  out=fopen(fn,"w");
  if(out==NULL)
     printf("cannot open file %s", fn);
  
  fprintf(out, "\nA : ");
  for(i=0; i<nA; i++){
    if(i==0){
      A[i]=rand()%d;
    }
    else{
      A[i]=A[i-1]+1+rand()%d;
    } 
    fprintf(out, " %d ",A[i]);
  }
  fclose(out);
  //==========================================================

  R = (int **) malloc(steps*sizeof(int *));
  nR = (int *) malloc(steps*sizeof(int));
  S = (int **) malloc(steps*sizeof(int *));
  nS = (int *) malloc(steps*sizeof(int));
  
  if(rank==0)
    printf("\nblillc");
  MPI_Barrier(MPI_COMM_WORLD);
  t0=MPI_Wtime();
  ring_init(A, nA, R, nR, S, nS, steps, MPI_COMM_WORLD);
  t1=MPI_Wtime();

  out=fopen(fn,"a");
  for(k=1; k<steps; k++){
    fprintf(out, "\nS%d : ", k);
    for(i=0; i<nS[k]; i++)
      fprintf(out, " %d ",S[k][i]);
   fprintf(out, "\nR%d : ", k);
    for(i=0; i<nR[k]; i++)
      fprintf(out, " %d ",R[k][i]);
  }
  fclose(out);

  if(rank==0)
    printf("\nstep\trmax\travg\trmin");
  for(k=1; k<steps; k++){
    MPI_Allreduce(&nR[k], &rmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&nR[k], &ravg, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nR[k], &rmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if(rank==0)
      printf("\n%d\t%d\t%d\t%d", k, rmax, ravg/size, rmin);
  } 

  free(nR);
  free(A);
  free(R);
  MPI_Finalize();
  return 0;
}

