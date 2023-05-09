/** @file:  mapmatc.c
    @brief  Sparse Matrix with Coarse Level 
    @note  Copyright (c) 2010-2012 APC CNRS Université Paris Diderot. This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/lgpl.html

    @note For more information about ANR MIDAS'09 project see http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html 
    @note ACKNOWLEDGMENT: This work has been supported in part by the French National  Research Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
    @author Pierre Cargemel
    @date   October 2012*/
#ifdef W_MPI 
 #include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mapmatc.h"

int CMatInit(CMat *A, int r, int *m, int *nnz, int **indices, double **values, int flag
#ifdef W_MPI
  ,MPI_Comm comm
#endif 
  ){
  int M, k, *tmp_indices;
  A->r    = r;						// set number of local rows
  A->m    = m;						// 
  A->nnz  = nnz;					// 
  A->disp = (int *) malloc((A->r+1)*sizeof(int));		// allocate disp array    
  A->disp[0]=0;
//  printf(" %d\t%d\t%d\n", A->m[0], A->nnz[0], A->disp[0]);
  for(k=1; k<=A->r; k++){
    A->disp[k]=A->disp[k-1]+A->m[k-1]*A->nnz[k-1];
  //  if(k!=A->r)
  //    printf(" %d\t%d\t", A->m[k], A->nnz[k]);
  //  printf(" %d\n", A->disp[k]);
  }
  A->indices = indices;					//
  A->values  = values;	
  /*int i, j;
  for(k=0; k<A->r; k++){
    for(i=0; i<A->m[k]*A->nnz[k]; i+=A->nnz[k]){
      for(j=0; j<A->nnz[k]; j++){
        printf(" %d ", A->indices[k][i+j]);
      }
    }
    printf("\n");
  }*/
  tmp_indices = (int *) malloc(A->disp[A->r]*sizeof(int));	// allocate a tmp copy of indices tab to sort    
  for(k=0; k<A->r; k++){
    memcpy(tmp_indices+A->disp[k], A->indices[k], A->m[k]*A->nnz[k]*sizeof(int));	// copy        
  }

  A->lcount = ssort(tmp_indices, A->disp[A->r], 0);  	  		    	// sequential sort tmp_indices (flag:3=counting sort)
  A->lindices = (int *) malloc((A->lcount)*sizeof(int));      //          
  memcpy(A->lindices, tmp_indices, (A->lcount) *sizeof(int));	// copy tmp_indices into lindices and free
  free(tmp_indices);     			

  for(k=0; k<A->r; k++){
    sindex(A->lindices, A->lcount, A->indices[k], A->nnz[k]*A->m[k]);    // transform indices tab in local indices tab 
  }
  /*for(k=0; k<A->r; k++){
    for(i=0; i<A->m[k]*A->nnz[k]; i+=A->nnz[k]){
      for(j=0; j<A->nnz[k]; j++){
        printf(" %d ", A->indices[k][i+j]);
      }
    }
    printf("\n");
  }*/
  //printf("cmat init 0\n");
#ifdef W_MPI
  A->comm = comm;					// link communivcator
  return CMatComShape(A, flag);  			// build communication scheme
#endif 
}


int CMatFree(CMat *A){
  free(A->disp);
  free(A->lindices);
#ifdef W_MPI
  if(A->flag != NONE){		//if necessary free communication tab
    if(A->R)			//
      free(A->R);		//
    if(A->nR)			//
      free(A->nR);		//
    if(A->S)			//
      free(A->S);		//
    if(A->nS)
      free(A->nS);
  }
#endif 
  return 0;
}



#ifdef W_MPI
int CMatComShape(CMat *mat, int flag){  
  //printf("commshape 0\n");
  int size;
  mat->flag = flag;
  MPI_Comm_size(mat->comm, &size);
  if(flag==BUTTERFLY){
    if(is_pow_2(size)==0){
    mat->flag=flag;
    mat->steps = log_2(size);
    }
    else{
      mat->flag=RING;
      mat->steps = size;
    }    
  }
  else if(flag==NONE){
    mat->flag=flag;
    return 0;
    }
  else{
    mat->flag=flag;
    mat->steps = size;
  }
  mat->S = (int** ) malloc(mat->steps * sizeof(int* ));                 /*<allocate sending maps tab*/
  mat->R = (int** ) malloc(mat->steps * sizeof(int* ));                 /*<allocate receiving maps tab*/
  mat->nS = (int* ) malloc(mat->steps * sizeof(int));                   /*<allocate sending map sizes tab*/
  mat->nR = (int* ) malloc(mat->steps * sizeof(int));                   /*<allocate receiving map size tab*/

  if(mat->flag == BUTTERFLY){
    butterfly_init(mat->lindices, mat->lcount, mat->R, mat->nR, mat->S, mat->nS, &(mat->com_indices), &(mat->com_count), mat->steps, mat->comm);
  } 
  else{
    ring_init(mat->lindices, mat->lcount, mat->R, mat->nR, mat->S, mat->nS, mat->steps, mat->comm);
    mat->com_count = mat->lcount;
    mat->com_indices = mat->lindices;
  }
  //printf("commshape 1\n");
 return 0;
}
#endif 



int CMatVecProd(CMat *A, double *xvalues, double* yvalues, int pflag){
  int i, j, k;
  int l;
  for(i=0; i<A->disp[A->r]; i++)                                                                     
      yvalues[i] = 0.0;
  l=0;   
  for(k=0; k<A->r; k++){						//coarse levels              
    for(i=0; i<A->m[k]; i+=A->nnz[k]){ 					//rows      
      for(j=0; j<A->nnz[k]; j++){					//non-zero per row
        yvalues[l] += A->values[k][i+j] * xvalues[A->indices[k][i+j]];   
      }
      l++;
    }
  }
  return 0;
}




int CTrMatVecProd(CMat *A, double *in_values, double* out_values, int pflag){
  int i, j, k;
  int l;
  int nSmax, nRmax;
  double *lvalues;

  lvalues = (double *) malloc(A->lcount *sizeof(double));    /*<allocate and set to 0.0 local values*/
  for(i=0; i < A->lcount; i++)
    lvalues[i]=0.0;

  l=0;   
  for(k=0; k<A->r; k++){						//coarse levels              
    for(i=0; i<A->m[k]; i+=A->nnz[k]){ 					//rows      
      for(j=0; j<A->nnz[k]; j++){					//non-zero per row
        lvalues[A->indices[k][i+j]] += A->values[k][i+j] * in_values[l];   
      }
      l++;
    }
  }
  memcpy(out_values, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
#ifdef W_MPI
  nRmax=0;
  nSmax=0;

  if(A->flag  == BUTTERFLY){                                  /*<branch butterfly*/
    for(k=0; k< A->steps; k++)                                /*compute max communication buffer size*/
      if(A->nR[k] > nRmax)
        nRmax = A->nR[k];
    for(k=0; k< A->steps; k++)
      if(A->nS[k] > nSmax)
        nSmax = A->nS[k];

    double *com_val;
    com_val=(double *) malloc( A->com_count *sizeof(double)); 
    for(i=0; i < A->com_count; i++){
      com_val[i]=0.0;
    } 
    m2m(lvalues, A->lindices, A->lcount, com_val, A->com_indices, A->com_count);
    butterfly_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
    m2m(com_val, A->com_indices, A->com_count, out_values, A->lindices, A->lcount);
    free(com_val);
  }
  else if(A->flag == RING){
    for(k=1; k< A->steps; k++)                                  /*compute max communication buffer size*/
      if(A->nR[k] > nRmax)
        nRmax = A->nR[k];
    
    nSmax = nRmax;  
    ring_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, out_values, A->steps, A->comm);
  }
  else if(A->flag == NONBLOCKING){
    ring_nonblocking_reduce(A->R, A->nR, A->S, A->nS, lvalues, out_values, A->steps, A->comm);
  }
  else if(A->flag == NOEMPTY){
    int ne=0;
    for(k=1; k< A->steps; k++)
      if(A->nR[k]!=0)
        ne++;  
    ring_noempty_reduce(A->R, A->nR, ne, A->S, A->nS, ne, lvalues, out_values, A->steps, A->comm);
  }
  else{
    return 1;
  }
#endif 
  free(lvalues);
  return 0;
}



