/** @file   test_cindex_0.c
    @author Pierre Cargemel
    @date   April 2012 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "../../src/mapmat/csort.h"


extern char *optarg;


int main(int argc, char *argv[]){
  char ch;
  int i, err;
  int nA, nT;
  int *A, *T;
  double t0, t1;
   
  MPI_Init(&argc,&argv); 
  
  nA=1;
  nT=1;

  while((ch = getopt(argc, argv, "n:m:")) != EOF){          /*set options*/
    switch(ch) {
      case 'n':
 	nA = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nA <= 0)) {
          return 1;
	}
      break;
      case 'm':
 	nT = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nT <= 0)) {
          return 1;
	}
      break;
    }
  }
  
  A = (int *) malloc(nA * sizeof(int));
  T = (int *) malloc(nT * sizeof(int));
  for(i=0; i<nA; i++)
    A[i]=rand()%nT;
  
  printf("\nA :");
  for(i=0; i<nA; i++)
    printf("\t%d",A[i]);
  
  for(i=0; i<nT; i++)
    T[i]=i;

  printf("\nT :");
  for(i=0; i<nT; i++)
    printf("\t%d",T[i]);
  
  
  t0 = MPI_Wtime();
//=======================Modification to compile without openmp
#ifdef W_OPENMP  
  omp_pindex(T, nT, A, nA);
#else
  sindex(T, nT, A, nA);
#endif
//=======================End of modification

  t1 =MPI_Wtime();

  printf("\nA :");
  for(i=0; i<nA; i++)
    printf("\t%d",A[i]);
  
  
  MPI_Finalize();

  return 0; 
}


