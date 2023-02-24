/** @file   test_als_0.c
    @brief  <b>Test algorithms for sets of integer : logical or </b>
    @n Usage : ./test_als_0 -n [elements in the first set] -m [elements in the second sets] -d [derive]. 
    The number of elements is the size of the array.
    The option -d enable to specify a derive parameter. Bigger is d, faster set elements grow.
    @author Pierre Cargemel
    @date   April 2012 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "../../src/mapmat/als.h"


extern char *optarg;


int main(int argc, char *argv[]){
  char ch;
  int i,d;
  int nA, nB, nAorB;
  int *A, *B, *AorB;

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
  printf("\nA :");
  for(i=0; i<nA; i++){
    if(i==0){
      A[i]=rand()%d;
    }
    else{
      A[i]=A[i-1]+1+rand()%d;
    } 
    printf("\t%d",A[i]);
  }

  B = (int *) malloc(nB * sizeof(int));                  //array B ascending ordered without redundancy
  printf("\nB :");
  for(i=0; i<nB; i++){
    if(i==0){
      B[i]=rand()%d;
    }
    else{
      B[i]=B[i-1]+1+rand()%d;
    }
    printf("\t%d",B[i]);
  }

  nAorB = card_or(A, nA, B, nB);
  printf("\ncard(AorB) = %d", nAorB);
  
  AorB = (int *) malloc(nAorB * sizeof(int));

  set_or(A, nA, B, nB, AorB);

  printf("\nAorB : ");
  for(i=0; i<nAorB; i++){
    printf("\t%d",AorB[i]);
  }
  printf("\n");
  
  free(A);
  free(B);
  free(AorB);
  
  return 0;
}


