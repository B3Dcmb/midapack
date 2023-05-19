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
#include <mpi.h>

// #include "../../src/mapmat/als.h"
// #include "midapack.h"

int old_card_or(int *A1, int n1, int *A2, int n2);
int old_set_or(int *A1, int n1, int *A2, int n2, int *A1orA2);

int new_card_or(int *A1, int n1, int *A2, int n2);
int new_set_or(int *A1, int n1, int *A2, int n2, int *A1orA2);

int main(int argc, char *argv[]){
  char ch;
  int i,d;
  int nA, nB, old_nAorB, new_nAorB;
  int *A, *B, *old_AorB, *new_AorB;

  nA=10;
  nB=8;
  // d=2;
  // while((ch = getopt(argc, argv, "n:m:d:")) != EOF){          /*set options*/
  //   switch(ch) {
  //     case 'n':
 	// nA = strtoul(optarg, 0, 10);
	// if ((errno == ERANGE) || (errno != 0 && nA <= 0)) {
  //         return 1;
	// }
  //     break;
  //     case 'm':
 	// nB = strtoul(optarg, 0, 10);
	// if ((errno == ERANGE) || (errno != 0 && nB <= 0)) {
  //         return 1;
	// }
  //     break;
  //     case 'd':
 	// d = strtoul(optarg, 0, 10);
	// if ((errno == ERANGE) || (errno != 0 && d <= 0)) {
  //         return 1;
	// }
  //     break;
  //   }
  // }

  
  A = (int *) malloc(nA * sizeof(int));                  //array A ascending ordered without redundancy
  printf("A :");
  for(i=0; i<nA; i++){
    // if(i==0){
    //   A[i]=rand()%d;
    // }
    // else{
    //   A[i]=A[i-1]+1+rand()%d;
    // }
    A[i]=i+5;
    printf("\t%d",A[i]);
  }
  printf(" --- nA %d \n", nA);

  B = (int *) malloc(nB * sizeof(int));                  //array B ascending ordered without redundancy
  printf("B :");
  for(i=0; i<nB; i++){
    // if(i==0){
    //   B[i]=rand()%d;
    // }
    // else{
    //   B[i]=B[i-1]+1+rand()%d;
    // }
    
    B[i]=i;// + 100;
    // B[i]=i + 1;
    printf("\t%d",B[i]);
  }
  printf(" --- nB %d \n", nB);

  old_nAorB = old_card_or(A, nA, B, nB);
  printf("card(old_AorB) = %d \n", old_nAorB);
  old_AorB = (int *) malloc(old_nAorB * sizeof(int));
  old_set_or(A, nA, B, nB, old_AorB);
  printf("old_AorB : ");
  for(i=0; i<old_nAorB; i++){
    printf("\t%d",old_AorB[i]);
  }
  printf("\n");
  new_nAorB = new_card_or(A, nA, B, nB);
  printf("card(new_AorB) = %d \n", new_nAorB);
  new_AorB = (int *) malloc(new_nAorB * sizeof(int));
  new_set_or(A, nA, B, nB, new_AorB);
  printf("new_AorB : ");
  for(i=0; i<new_nAorB; i++){
    printf("\t%d",new_AorB[i]);
  }
  printf("\n");
  
  free(A);
  free(B);
  free(old_AorB);
  free(new_AorB);
  
  return 0;
}


int old_card_or(int *A1, int n1, int *A2, int n2){
  int i=0, j=0, k= 0;
  printf("### Entering old_card_or \n");
  while( i<n1 || j<n2){
    // printf("### while : i %d n1 %d --- j %d n2 %d \n", i, n1, j, n2);
    if(A1[i] < A2[j]){
      if(i<n1){ i++; }
      else{ j++; }
    }
    else if(A1[i] > A2[j]){
      if(j<n2){ j++; }
      else{ i++; }
    }
    else{
      if(i<n1){ i++; }
      if(j<n2){ j++; }
    }
    k++; 
  }
  return k;
}

int old_set_or(int *A1, int n1, int *A2, int n2, int *A1orA2)
{
  int i=0, j=0, k=0;
  printf("~~~~ Entering old_set_or \n");
  while( i<n1 || j<n2){
    // printf("~~~~ while : i %d n1 %d --- j %d n2 %d \n", i, n1, j, n2);
    if(A1[i] < A2[j]){
      if(i<n1){
        A1orA2[k]=A1[i];
        i++;
      }
      else{
        A1orA2[k]=A2[j];
        j++;
      }
    }
    else if(A1[i] > A2[j]){
      if(j<n2){
        A1orA2[k]=A2[j];
        j++;  
      }
      else{
        A1orA2[k]=A1[i];
        i++;
      }
    }
    else{
      A1orA2[k]=A1[i];
      i++;  
      j++;    
    }
    k++;
  }
  return k;
}


int new_card_or(int *A1, int n1, int *A2, int n2){
  int i=0, j=0, k= 0;
  printf("<<<< Entering new_card_or \n");
  // while( (i!=n1-1) || (j!=n2-1)){
  while( (i<n1) || (j<n2)){
    // printf("<<<< while : i %d n1 %d --- j %d n2 %d \n", i, n1, j, n2);
    if ((i==n1) || (j==n2)){
      if(i==n1){
        j++;
      }
      else{
        i++;
      }
    }
    else if(A1[i] < A2[j]){
      if(i<n1){ i++; }
      else{ j++; }
    }
    else if(A1[i] > A2[j]){
      if(j<n2){ j++; }
      else{ i++; }
    }
    else{
      if(i<n1){ i++; }
      if(j<n2){ j++; }
    }
    k++;
  }
  // if (A1[i] != A2[j]){
  //   k += 2; // Count the last two elems
  // }
  // else {
  //   k++; // Count the last elem
  // }
  return k;
}

int new_set_or(int *A1, int n1, int *A2, int n2, int *A1orA2)
{
  int i=0, j=0, k=0;
  printf(">>>> Entering new_set_or \n");
  while( (i<n1) || (j<n2)){
    // printf(">>>> while : i %d n1 %d --- j %d n2 %d \n", i, n1, j, n2);
    if ((i==n1) || (j==n2)){
      if(i==n1){
        A1orA2[k]=A2[j];
        j++;
      }
      else{
        A1orA2[k]=A1[i];
        i++;
      }
    }
    else if(A1[i] < A2[j]){
      if(i<n1){
        A1orA2[k]=A1[i];
        i++;
      }
      else{
        A1orA2[k]=A2[j];
        j++;
      }
    }
    else if(A1[i] > A2[j]){
      if(j<n2){
        A1orA2[k]=A2[j];
        j++;  
      }
      else{
        A1orA2[k]=A1[i];
        i++;
      }
    }
    else{
      A1orA2[k]=A1[i];
      i++;  
      j++;    
    }
    k++;
  }
  return k;
}