/** @file   cindex.c
    @brief  Indexing subroutines implemetation
    @note  Copyright (c) 2010-2012 APC CNRS Université Paris Diderot. This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/lgpl.html
    @note For more information about ANR MIDAS'09 project see http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html 
    @note ACKNOWLEDGMENT: This work has been supported in part by the French National  Research Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
    @author Pierre Cargemel
    @date   May 2012*/


#include <stdlib.h>
/** Sequential reindexing  
    @param T monotony array
    @param nT number of index
    @param A tab to reindex
    @param nA number of element to reindex
    @return array of indices 
    @ingroup matmap_group22*/
int sindex(int *T, int nT, int *A, int nA){
  int i, tmp;
  i=0;
  for(i=0; i<nA; i++){
    tmp = A[i];
    A[i] =dichotomy(nT, T, tmp);
  }
}


#ifdef W_OPENMP
/** Multithread (OpenMP) reindexing 
    @param T monotony array
    @param nT number of index
    @param A tab to reindex
    @param nA inumber of element to reindex
    @return array of indices 
    @ingroup matmap_group22*/
int omp_pindex(int *T, int nT, int *A, int nA){
//  printf("\nomp_pindex");
  int i;
  int *count, *disp;
  int q, r;
  int tid, nths;

  #pragma omp parallel private(tid) shared(nths)
  {//---fork---just to get the number of threads
    nths = omp_get_num_threads();
    tid = omp_get_thread_num();
//    printf("\ntid %d nths %d", tid, nths);
  }//---join--- 

  q = nA/nths;
  r = nA%nths; 

  count = (int *) malloc(nths *sizeof(int));
  disp = (int *) malloc(nths *sizeof(int));

  for(i=0; i<nths; i++){
    if(i<r){
      count[i] = q+1;
    }
    else{
      count[i] = q;
    }
  }
  
  disp[0] = 0;
  for(i=0; i<nths-1; i++){
    disp[i+1] = disp[i] + count[i];
  }
  
  #pragma omp parallel private(tid) shared(T, nT, A, disp, count)
  {//---fork---1st step, sort on local chunk
    tid = omp_get_thread_num();
    sindex(T, nT, A+disp[tid], count[tid]);
  }//---join---
  free(count);
  free(disp); 
  return 0; 
}
#endif



/** dichotmic search of an integer in a monotony array
    @param number elemnent array of values
    @param monotony array 
    @param element to search
    @return index of searched element*/ 
int dichotomy(int nT, int *T, int e){
  int min, max, pivot;
  min=0;
  max=nT-1;
  pivot=(max-min)/2;
  while(e != T[pivot] && max > min ){
    if(T[pivot]<e){
      min=pivot+1;
    }
    else{
      max=pivot;
    }
    pivot= min + (max-min)/2;
  }
  return pivot;
}

