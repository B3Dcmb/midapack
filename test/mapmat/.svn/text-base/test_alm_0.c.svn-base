/** @file   test_alm_0.c
    @brief  <b>Test subroutines handling maps and submaps. </b>
    @n Usage : ./test_alm_0 -n [elements in the first map] -m [elements in the second map] -d [derive]. 
    The option -d enable to specify a derive parameter. Bigger is d, faster set elements grow.
    @author Pierre Cargemel
    @date   April 2012 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "../../src/mapmat/als.h"
#include "../../src/mapmat/alm.h"


extern char *optarg;


int main(int argc, char *argv[]){
  char ch;
  int i,d;
  int nA, nB, nAandB;
  int *A, *B, *AandB;
  double *vA, *vB, *vAandB;

  nA=1;
  nB=1;
  d=2;
  while((ch = getopt(argc, argv, "n:m:d:")) != EOF){          /*set options*/
    switch(ch) {
      case 'n':
 	nA = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nA <= 0)) {
          return 1;
	}
      break;
      case 'm':
 	nB = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nB <= 0)) {
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

  
  A = (int *) malloc(nA * sizeof(int));                  //array A ascending ordered without redundancy
  vA = (double *) malloc(nA * sizeof(double));                  //associated values
  printf("\n(A, vA) :");
  for(i=0; i<nA; i++){
    if(i==0){
      A[i]=rand()%d;
    }
    else{
      A[i]=A[i-1]+1+rand()%d;
    } 
    vA[i]=1.0;
    printf(" (%d, %.2f) ",A[i], vA[i]);
  }

  B = (int *) malloc(nB * sizeof(int));                  //array B ascending ordered without redundancy
  vB = (double *) malloc(nB * sizeof(double));                  //associated values
//  printf("\n(B, vB) :");
  for(i=0; i<nB; i++){
    if(i==0){
      B[i]=rand()%d;
    }
    else{
      B[i]=B[i-1]+1+rand()%d;
    }
    vB[i]=1.0;
//    printf(" (%d, %.2f) ",B[i], vB[i]);
  }

  nAandB = card_and(A, nA, B, nB);
  printf("\ncard(AandB) = %d", nAandB);
  
  AandB = (int *) malloc(nAandB * sizeof(int));
  vAandB = (double *) malloc(nAandB * sizeof(double));                  //associated values

  set_and(A, nA, B, nB, AandB);

  printf("\n(AandB, vAandB) : ");
  for(i=0; i<nAandB; i++){
    vAandB[i]=1.0; 
    printf(" (%d,%.2f) ",AandB[i]);
  }

  subset2map(A, nA, AandB, nAandB);
  
  printf("\nmap AandB / A : ");
  for(i=0; i<nAandB; i++){
    printf(" %d ",AandB[i]);
  }

  s2m_sum(vA, vAandB, AandB, nAandB);

  printf("\n(A, vA) : ");
  for(i=0; i<nA; i++){
    printf(" (%d, %.2f) ", A[i], vA[i]);
  }

  m2s(vA, vAandB, AandB, nAandB);

  printf("\nvsubmap : ");
  for(i=0; i<nAandB; i++){
    printf(" %.2f ", vAandB[i]);
  }

  printf("\n");
  
  free(A);
  free(vA);
  free(B);
  free(vB);
  free(AandB);
  free(vAandB);
  
  return 0;
}


