/** @file   test_als_1.c
    @brief  <b>Test algorithms for sets of integer : card, merge </b>
    @n Usage : ./test_als_1 -n [elements in the set] -d [derive]. 
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
  nB=0;
  d=2;
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

  
  A = (int *) malloc(nA * sizeof(int));                  //array A ascending ordered with possible redundancy
  printf("\nA :");
  for(i=0; i<nA; i++){
    if(i==0){
      A[i]=rand()%d;
    }
    else{
      A[i]=A[i-1]+rand()%d;
    } 
    printf("\t%d",A[i]);
  }

  nB = card(A, nA);
  printf("\ncard(B) = %d", nB);
  
  B = (int *) malloc(nB * sizeof(int));

  merge(A, nA, B);

  printf("\nB :");
  for(i=0; i<nB; i++){
    printf("\t%d",B[i]);
  }
  printf("\n");
  
  free(A);
  free(B);
  
  return 0;
}


