/** @file als.c
    @brief Implementation of subroutines handling sets of inidices (union, intersection, merge, cardinal ...)
    @note  Copyright (c) 2010-2012 APC CNRS Université Paris Diderot. This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/lgpl.html
    @note For more information about ANR MIDAS'09 project see http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html 
    @note ACKNOWLEDGMENT: This work has been supported in part by the French National  Research Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
    @author Pierre Cargemel
    @date April 2012*/

#include <stdlib.h>

/** Compute cardinal(A) 
    Perform a loop onto elements of a set, counting different elements.
    The array should be in ascending order with possible redundancy. 
    @param n number of elemnets in A
    @param A set of indices 
    @return number of elements in A*/ 
int card(int *A, int nA){
  int i;
  int tmp=A[0];
  int c=1;
  for(i=1; i<nA; i++){
    if(A[i] != tmp){
      c++;  
      tmp=A[i];
    }
  } 
  return c;
}


/** Merge redundant elements.
    Fill a new array containing elements of A, merging redundant elements.
    The array A should be in ascending order with possible redundancy.
    The array B has to be allocated.
    @param nA number of elemnets in A
    @param A set of indices 
    @param B output array
    @return void*/ 
void merge(int *A, int nA, int *B){
  int i=0, j=0;
  B[0]=A[0];
  for(i=1; i<nA; i++){
    if(A[i] != B[j]){
      j++;  
      B[j]=A[i];
    }
  } 
}


/** Compute \f$ card(A_1 \cup A_2) \f$ 
    A1 and A2 should be two ascending ordered monotmony sets. 
    of size n1 and n2.
    @param n1 number of elemnets in A1
    @param A1 set of indices 
    @param n2 number of elemnets in A2
    @param A2 set of indices
    @return size of the union 
    @ingroup matmap_group22*/
int card_or(int *A1, int n1, int *A2, int n2){
  int i=0, j=0, k=0;
  while( (i<n1) || (j<n2)){
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
  return k;
}

/** Compute \f$ card(A_1 \cap A_2) \f$ 
    A1 and A2 should be two ascending ordered monotony sets, 
    of size n1 and n2.
    @param n1 number of elemnets in A1
    @param A1 set of indices 
    @param n2 number of elemnets in A2
    @param A2 set of indices
    @return size of the intersection 
    @ingroup matmap_group22*/
int card_and(int *A1, int n1, int *A2, int n2){
  int i=0, j=0, k= 0;
  while( i<n1 && j<n2){
    if(A1[i] < A2[j]){
      i++;  
    }
    else if(A1[i] > A2[j]){
      j++;
    }
    else{
      k++; 
      i++;
      j++;
    }
  }
  return k;
}

/** Compute \f$  A1 \cup A_2 \f$ 
    A1 and A2 should be two ascending ordered sets. 
    It requires the sizes of these two sets, n1 and n2.
    A1andA2 has to be previouly allocated. 
    @param n1 number of elemnets in A1
    @param A1 set of indices 
    @param n2 number of elemnets in A2
    @param A2 set of indices
    @param address to the set A1orA2
    @return number of elements in A1orA2 
    @ingroup matmap_group22*/
int set_or(int *A1, int n1, int *A2, int n2, int *A1orA2)
{
  int i=0, j=0, k=0;
  while( (i<n1) || (j<n2)){
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

/** Compute \f$ A_1 \cap A_2 \f$
    A1 and A2 should be two monotony sets in ascending order. 
    It requires the sizes of these two sets, n1 and n2.
    A1andA2 has to be previously allocated. 
    @param n1 number of elemnets in A1
    @param A1 set of indices 
    @param n2 number of elemnets in A2
    @param A2 set of indices
    @param address to the set A1andA2
    @return number of elements in A1andA2 
    @ingroup matmap_group22*/
int set_and(int *A1, int n1, int *A2, int n2, int *A1andA2){
  int i=0, j=0, k= 0;
  while( i<n1 && j<n2){
    if(A1[i] < A2[j]){
      i++;  
    }
    else if(A1[i] > A2[j]){
      j++;
    }
    else{
      A1andA2[k]=A1[i];
      k++; 
      i++;
      j++;
    }
  }
  return k;
}

/** @brief Compute map A1 and A2 / A1 
    @param n1 number of elemnets in A1
    @param A1 set of indices 
    @param n2 number of elemnets in A2
    @param A2 set of indices
    @param address to the set A1andA2
    @return number of elements in A1andA2*/ 
int map_and(int *A1, int n1, int *A2, int n2, int *mapA1andA2){
  int i=0, j=0, k= 0;
  while( i<n1 && j<n2){
    if(A1[i] < A2[j]){
      i++;  
    }
    else if(A1[i] > A2[j]){
      j++;
    }
    else{
      mapA1andA2[k]=i;
      k++; 
      i++;
      j++;
    }
  }
  return k;
}


/** Transform a subset into a mapper array 
    Parse a subset and replace element by his index in the larger set.
    A and subA should be two monotony sets in ascending ordered. 
    subA has to belong to A. 
    @param A a set of indices(monotony)
    @param nA number of elemnets in A
    @param subA a subset of A(monotony)
    @param nsubA number of elemnets in A
    @return void*/ 
void subset2map(int *A, int nA, int *subA, int nsubA){
  int i=0, j=0;
  while( i<nA && j<nsubA){
    if(A[i] < subA[j]){
      i++;  
    }
    else{
      subA[j]=i;
      i++;
      j++;
    }
  }
}

