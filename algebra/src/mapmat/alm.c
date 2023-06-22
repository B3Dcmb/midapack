/** @file alm.c
    @brief Implementation of subroutines handling maps, distributions or functions.
    That means, almost all structures describes as sets of indices associated to sets of values).
    @note  Copyright (c) 2010-2012 APC CNRS Université Paris Diderot. This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/lgpl.html
    @note For more information about ANR MIDAS'09 project see http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
    @note ACKNOWLEDGMENT: This work has been supported in part by the French National  Research Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
    @author Pierre Cargemel
    @date April 2012*/

#include <stdlib.h>

/** Set some map values into a submap values array
    @param mapval array of values
    @param submapval  array of values
    @return array of indices*/
void m2s(double *mapval, double *submapval, int *subset, int count){
  int i;

//#pragma omp parallel for
  for(i=0; i< count; i++){
    submapval[i]=mapval[subset[i]];
  }
}


/** @brief Local mat vec prod
    @param indm table of integers apval array of values
    @param submapval  array of values
    @return void*/
void lmatvecprod(int *ind, double *val, int m, int nnz, double *in, double *out){
  int i, j, k;
  k=0;
  for(i=0; i<m; i++){                               /*<local transform reduce*/
    for(j=0; j<nnz; j++){
      out[i]+=val[k]*in[ind[k]];
      k++;
    }
  }
}


/** @brief Sum  submap values the submap values array
    @param mapval array of values
    @param submapval  array of values
    @return array of indices*/
void s2m_sum(double *mapval, double *submapval, int *subset, int count){
  int i;
//#pragma omp parallel for
  for(i=0; i< count; i++){
    mapval[subset[i]] += submapval[i];
  }
}


/** @brief assign submap values the submap values array
    @param mapval array of values
    @param submapval  array of values
    @return array of indices*/
void s2m(double *mapval, double *submapval, int *subset, int count){
  int i;
  for(i=0; i< count; i++){
    mapval[subset[i]] = submapval[i];
  }
}

/** @brief Sum  submap values the submap values array
    @return void*/
void cnt_nnz_dot_prod(double *out, double *in, int cnt, int *ind, double *val, int nnz){
  int i, j, k;
  k=0;
  for(i=0; i<cnt; i++)                                   /*<local transform reduce*/
    for(j=0; j<nnz; j++)
      out[ind[k]]+=val[k]*in[i];
}

#if OPENMP
/** @brief Sum  submap values the submap values array
    @return void */
void omp_cnt_nnz_dot_prod(double *out, double *in, int cnt, int *ind, double *val, int nnz){
  int i, j, k;
  k=0;
  for(i=0; i<cnt; i++)                                   /*<local transform reduce*/
    for(j=0; j<nnz; j++)
      out[ind[k]]+=val[k]*in[i];
}
#endif

/** Function m2m for "map to map"
    Extract values from one map (A1, vA1), and for each pixel shared with an other map (A2, vA2),
    assign pixel value in vA1 and to pixel value in vA2.
    @return a number of elements shared between A1 and A2
    @sa m2m_sum
    @ingroup matmap_group22*/
int m2m(double *vA1, int *A1, int n1, double *vA2, int *A2, int n2){
  int i=0, j=0, k= 0;
  while( i<n1 && j<n2){
    if(A1[i] < A2[j]){
      i++;
    }
    else if(A1[i] > A2[j]){
      j++;
    }
    else{
      vA2[j]=vA1[i];
      k++;
      i++;
      j++;
    }
  }
  return k;
}

/** Function m2m_sum for "sum map to map"
    Extract values from one map (A1, vA1), and for each pixel shared with an other map (A2, vA2),
    sum pixel value in vA1 to pixel value in vA2.
    @return a number of elements shared between A1 and A2
    @sa m2m
    @ingroup matmap_group22*/
int m2m_sum(double *vA1, int *A1, int n1, double *vA2, int *A2, int n2){
  int i=0, j=0, k= 0;
  while( i<n1 && j<n2){
    if(A1[i] < A2[j]){
      i++;
    }
    else if(A1[i] > A2[j]){
      j++;
    }
    else{
      vA2[j]+=vA1[i];
      k++;
      i++;
      j++;
    }
  }
  return k;
}

/** Function m2m_sum_i for "sum map to map" (integer version)
    Extract values from one integer map (A1, vA1), and for each pixel shared with an other integer map (A2, vA2),
    sum pixel value in vA1 to pixel value in vA2.
    @return a number of elements shared between A1 and A2
    @sa m2m
    @ingroup matmap_group22*/
int m2m_sum_i(int *vA1, int *A1, int n1, int *vA2, int *A2, int n2){
  int i=0, j=0, k= 0;
  while( i<n1 && j<n2){
    if(A1[i] < A2[j]){
      i++;
    }
    else if(A1[i] > A2[j]){
      j++;
    }
    else{
      vA2[j]+=vA1[i];
      k++;
      i++;
      j++;
    }
  }
  return k;
}
