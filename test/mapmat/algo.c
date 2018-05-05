/** @file   algo.c
    @brief  <b>Subroutines implemetation</b>
    @author Pierre Cargemel
    @date   November 2011*/
#include "algo.h"

void sort(int nb, int *array){
  int i, j;
  int tmp;
  for(i=0; i<nb-1 ; i++){
    tmp = array[i+1];
    j=i;
    while(j != -1 && tmp < array[j] ){
      array[j+1]=array[j];
      array[j]=tmp;
      j--;
    } 
  }
}


int count(int nb, int *array){
  int i;
  int tmp=array[0];
  int cpt=1;
  for(i=1; i<nb; i++){
    if(array[i] != tmp){
      cpt++;  
      tmp=array[i];
    }
  } 
  return cpt;
}
 

int shared(int nb1, int *array1, int nb2, int *array2){
  int i=0, j=0 ;
  int cpt = 0;
  while( i<nb1 && j<nb2 ){        //while end of tab is not reached 
    if(array1[i] < array2[j]){    
      i++;                        
    }
    else if(array1[i] > array2[j]){
      j++;
    }
    else{
      cpt++; 
      i++;
      j++;
    }
  } 
  return cpt;
}


void intersection(int nb1, int *array1, int nb2, int *array2, int nb3, int* array3){
  int i=0, j=0, k= 0;
//  printf("here");
  while( i<nb1 && j<nb2 && k<nb3){
    if(array1[i] < array2[j]){
      i++;  
    }
    else if(array1[i] > array2[j]){
      j++;
    }
    else{
      array3[k]=i;
      k++; 
      i++;
      j++;
    }
  } 
}


void reduce(int lcount, int *lindices, double *lvalues, int pcount, int* pindices, double* pvalues){
  int i=0, j=0 ;
  int cpt = 0;
  while( i<lcount && j<pcount ){
    if(lindices[i] < pindices[j]){
      i++;  
    }
    else if(lindices[i] > pindices[j]){
      j++;
    }
    else{
//      printf("(%d %lf) + (%d %lf) = ",lindices[i] ,lvalues[i], pindices[j], pvalues[j]); 
      lvalues[i]+=pvalues[j];
  //    printf("%lf\n",lvalues[i]); 
      i++;
      j++;
    }
  }
}



void condense(int nbin, int *arrayin, int nbout, int *arrayout){
  int i, j;
  j=0;
  arrayout[j]=arrayin[j];
  for(i=0; i<nbin; i++){
    if(arrayin[i]==arrayout[j]){
    }
    else{
      j++;
      arrayout[j]=arrayin[i];
    }
  }
}




void print(int nb, int *array){
  int i;
  for(i=0; i<nb; i++)
    printf("%d\t%d\n", i, array[i]);
}


int partition(int *fi, int *wi, int N, int me, int us){
  int r, k;
  k = N / us;
  r = N - k*us;
  if( me < r){
    *fi = (k+1) * me;
    *wi = k+1;
  }
  else{
    *fi = r*(k+1) + k*(me-r);
    *wi = k;
    }
  return 0;
}
