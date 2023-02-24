/** @file   test_csort_0.c
    @brief  <b>Test sequential sort </b>
    @note       (c) Copyright 2012 APC CNRS Université Paris Diderot
    Usage : ./test_csort_0 -n [number of elements] -m [max] -f [method]. 
    The number of elements is the size of the array to sort.
    The max specify elements belong to the set [0 max] (elements are randomly generated).
    the option -f enable to specify a sort algorithm (quick, counting, shell, bubble, insetion).
    After the sorting routine has been called, the program checked result and print some detalis.
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
  char ch, *method_name;
  int i, d, m, err, method_flag;
  int nA, n;
  int *A;
  double t0, t1;
   
  MPI_Init(&argc,&argv); 
  
  method_flag=0;
  nA=1;
  d=2;
  while((ch = getopt(argc, argv, "n:m:f:")) != EOF){          /*set options*/
    switch(ch) {
      case 'n':
 	nA = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && nA <= 0)) {
          return 1;
	}
      break;
      case 'm':
 	m = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && m <= 0)) {
          return 1;
	}
      break;
      case 'f':
        method_name = strdup(optarg);
  //      printf("\nmethod name %s", method_name);
        if(strcmp(method_name,"quick")==0){
           method_flag=0;
        }
        else if(strcmp(method_name,"bubble")==0){
           method_flag=1;
        }
        else if(strcmp(method_name,"insertion")==0){
           method_flag=2;
        }
        else if(strcmp(method_name,"counting")==0){
           method_flag=3;
        }
        else if(strcmp(method_name,"shell")==0){
           method_flag=4;
        }
        else{
           method_flag=0;
        } 
      break;
    }
  }
  
  A = (int *) malloc(nA * sizeof(int));
//  printf("\nA :");
  for(i=0; i<nA; i++){
    A[i]=rand()%m;
//    printf("\t%d",A[i]);
  }
  
  t0 = MPI_Wtime();
  
  n = ssort(A, nA, method_flag);

  t1 =MPI_Wtime();

//  printf("\nA :");
//  for(i=0; i<n; i++){
//    printf("\t%d",A[i]);
//  }
  
  err = monotony(A, n);
 
  free(A);
  
  if(err==0){
    printf("%s_sort/nA: %d /n: %d /m: %d /time: %lf /completed\n", method_name, nA, n, m, t1-t0);
  }
  else{
    printf("failed\n");
  }

  MPI_Finalize();

  return err;
}


