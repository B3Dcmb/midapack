/** @file   test_als_2.c
    @brief  <b>Test algorithms for sets of integers : logical and </b>
    @note       (c) Copyright 2012 APC CNRS Université Paris Diderot
    @n Usage : ./test_als_2 -n [elements in the first set] -m [elements in the second sets] -d [derive]. 
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
  int nA, nB, nAandB;
  int *A, *B, *AandB;

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

  nAandB = card_and(A, nA, B, nB);
  printf("\ncard(AandB) = %d", nAandB);
  
  AandB = (int *) malloc(nAandB * sizeof(int));

  set_and(A, nA, B, nB, AandB);

  printf("\nAandB : ");
  for(i=0; i<nAandB; i++){
    printf("\t%d",AandB[i]);
  }
  printf("\n");
  
  free(A);
  free(B);
  free(AandB);
  
  return 0;
}


