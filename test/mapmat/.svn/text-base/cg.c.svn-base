/** @file   pcg.c
  
    @author Pierre Cargemel
    @date   December 2011 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include "midapack.h"




/**CGNE: Conjugate Gradient for Normal Equation
 r = A'b - A'A*x;
 if( r'r < tol )
   return;
 w = r;
 for k=1:K
   z = A'A*w;
     a = (r'*r)/(w'*z);
   x = x + a*w;
   r = r - a*z;
   B = (r'*r)/(w'*z);
   if( r'r < tol )
     break;
   w = r + B*w;
 end
*/
int CGNE(Mat A, double *x, double*b, double tol, int K){

  int 		i, j, k, err;
  int		m, n, rank, size;
  double 	*Aw;			//time domain vectors 
  double	*z ,*r, *w, *pond;	//map domain vectors
  double        a, B, C, res, resold;	//coef and residual
  double 	localreduce;		//reduce buf
  double	st, t;				//timer

  m=A.m;					//number of local time samples
  n=A.lcount;					//number of local pixels
  MPI_Comm_rank(A.comm, &rank);			//
  MPI_Comm_size(A.comm, &size);			//
  
  r = (double *) malloc(n*sizeof(double));	//residual,
  w = (double *) malloc(n*sizeof(double));	//direction,   
  z = (double *) malloc(n*sizeof(double));	//A'A direction
  pond = (double *) malloc(n*sizeof(double));	//weigh for dot product,
  Aw = (double *) malloc(m*sizeof(double));	//A direction

  for(j=0; j<n; j++)   			        // 
    pond[j]=1.0;				//

  greedyreduce(&A, pond);
  for(j=0; j<n; j++){   		        // 
    pond[j]=1.0/pond[j];			//
  }
  
  //init CG 
  MatVecProd(&A, x, Aw, 0);		//r = A'b -A'Ax
  for(i=0; i<m; i++)   		        // 
    Aw[i]=b[i]-Aw[i];			//
  TrMatVecProd(&A, Aw, r, 0);		//

  res=0.0;				//lucky exit
  localreduce=0;			//
  for(j=0;j<n ;j++){ 
    localreduce+=r[j]*r[j]*pond[j];
//    printf("rank %d l[%d]= %lf\n", rank, j,localreduce);
  }
  MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(rank==0)				//
     printf("\ninitial residual %lf", res); 
  if(res<tol){				//
    if(rank ==0)			//
      printf("\n--> converged (%lf < %lf)\n", res, tol);
  return 0;				//
  } 
 
  st=MPI_Wtime();			//start timer

  for(j=0; j<n; j++)   		        //w = -r 
    w[j]=r[j];

  //starts loop CG 
  for(k=0; k<K ; k++){

    MatVecProd(&A, w, Aw, 0);		// 
    //for(i=0;i<m ;i++){   
    //  printf("rank %d [%d] Aw=%lf \n", rank, i, Aw[i]);
    //}
    //printf("\n");
    TrMatVecProd(&A, Aw, z, 0);		//z = A'A*w
    C=0.0;				//
    localreduce=0.;			//
    //for(i=0;i<m ;i++)    
    //  localreduce+=Aw[i]*Aw[i];
    for(j=0;j<n ;j++){   
      //printf("rank %d [%d] w=%lf z=%lf\n", rank, j, w[j] , z[j]);
      localreduce+=w[j]*z[j]*pond[j];
    }
    MPI_Allreduce(&localreduce, &C, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    a=res/C;				//

    for(j=0; j<n; j++)			//x = x + a*w
       x[j]= x[j] + a*w[j];		// 

    //TrMatVecProd(&A, Aw, z, 0);		//z = A'A*w

    for(j=0; j<n; j++)			//r = r - a*z
       r[j]= r[j] - a*z[j];		// 

    resold=res;				//
    res=0.0;				//res = n_r'*A*r
    localreduce=0.;			//
    for(j=0;j<n ;j++)    
      localreduce+=r[j]*r[j]*pond[j];
    MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    
    t=MPI_Wtime();			//stop timer
    if(rank==0)				//print iterate info
      printf("\nk=%d res=%lf alpha=%lf t=%lf", k, res, a, t-st); 
    if(res<tol){			//
      if(rank ==0)			//
        printf("\n--> converged (%lf < %lf)\n", res, tol);
      break;				//
    }
    st=MPI_Wtime();			//start timer

    B=res/resold;			//B=res/resold;
    for(j=0; j<n; j++)			// w = -r + B*w;
       w[j]= r[j] + B*w[j];		// 
  }					//end loop CG 

  free(z);				//free vectors
  free(r);				//
  free(Aw);				//
  free(w);				//
  free(pond);				//

  if(k==K){				//check unconverged
    if(rank ==0)
      printf("\n--> unconverged, max iterate reached (%lf > %lf)\n", res, tol);
    return 1;
  }
  return 0;
}

