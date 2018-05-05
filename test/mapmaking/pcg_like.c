// Midapack library
// Part of the Mapmaking example - version 1.1b, September 2012
// new pcg from scratch including the Toeplitz matrix

/** @file   pcg.c
    @author Frederic Dauvergne, Pierre Cargemel
    @date   September 2012 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <midapack.h>


/* Preconditioned Conjugate Gradient algorithm for the Map-Making equation, i.e.
   solve the system A^t N^-1 A x = A^t N-1 b,
   where A is unpointing matix from the mapmat module and N^-1 is a symmetric piecewise Toeplitz 
   from Toeplitz module.

   The algoritm is slightly different from the classical PCG. Matrix are applied succesively
   at each iterate to take advantage of the sparsity. The dot product is perform in the 
   "well-ballanced" time domain. 

   Occupancy  : 1 vector in map domain, 3 vectors in time domain.
   Complexity : at each iterate operation, we use 1 Toeplitz product, 3 dots product,
                4 axpy, 3 unpointing product (MatVecProd) and 1 pointing product (TrMatVecProd). 
*/

int PCG_GLS_like(Mat A, Tpltz Nm1, double *x, double*b, double tol, int K)//, double *c)
{

  int 		i, j, k ;			// some indexes 
  int		m, n, rank, size;
  double 	localreduce;			//reduce buffer
  double	st, t;				//timers

  m=A.m;					//number of local time samples
  n=A.lcount;					//number of local pixels
  MPI_Comm_rank(A.comm, &rank);			//
  MPI_Comm_size(A.comm, &size);			//
  

  double *_g, *ACg, *Ah, *Nm1Ah ;		// time domain vectors 
  double *g, *Cg, *h ;		      		// map domain vectors
  double *AtNm1Ah ;                    		// map domain
  double res, g2, g2p, ro, gamma, coeff ;	// scalars 


  //map domain
  h = (double *) malloc(n*sizeof(double));      //descent direction  

  //time domain
  _g = (double *) malloc(m*sizeof(double));   
  Ah = (double *) malloc(m*sizeof(double));    
  Nm1Ah = (double *) malloc(m*sizeof(double)); 

  ACg=Ah;  
  g = Ah;  //gradient, work because m>n 
  Cg=Nm1Ah;   //ok because m>n


  printf("n=%d, m=%d, A.nnz=%d \n", n, m, A.nnz );

//Init CG descent
  st=MPI_Wtime();			//start timer

  double *c;
  c = (double *) malloc(n*sizeof(double)); 


// for no preconditionner:
//  for(j=0; j<n; j++)                    // 
//    c[j]=1.;


  precondjacobilike( A, Nm1, c);
//  precondjacobilike_avg( A, Nm1, c);


  if (rank==0) {
    for(i=0; i<min(A.lcount, 10); i++)
    printf("[rank %d] c[%d]=%f \n", rank, i, c[i]);
  }


  MatVecProd(&A, x, _g, 0);		//
  for(i=0; i<m; i++)			//
    _g[i] = _g[i] - b[i];		//

  stbmmProd(Nm1, _g);			// _g = Nm1 (Ax-b)
  TrMatVecProd(&A, _g, g, 0);		//  g = At _g 


  for(j=0; j<n; j++)                    // 
    Cg[j]=c[j]*g[j]; 			//  Cg = C g  with C = Id
 
  for(j=0; j<n; j++)   		        //  h = -Cg 
    h[j]=-Cg[j];  


   MatVecProd(&A, Cg, ACg, 0);		//  ACg = A Cg

  g2=0.0;  				//g2 = "res"                          
  localreduce=0.0;                    
  for(i=0; i<m; i++)                    //  g2 = (ACg , _g)    
    localreduce+= ACg[i] * _g[i]; 	  
  MPI_Allreduce(&localreduce, &g2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//Test if already converged
   if(rank==0)                          //
     printf("k=%d res=%e \n", (-1), g2);

   if(g2<tol){                         //
      if(rank==0)                      //
        printf("--> converged (%e < %e)\n", g2, tol);
        k=K;//to not enter inside the loop
   }

  double tol2rel = tol;//*tol*g2;


// PCG Descent Loop
  for(k=0; k<K ; k++){               


    MatVecProd(&A, h, Ah, 0);		// Ah = A h
			
    for(i=0; i<m; i++)			
      Nm1Ah[i] = Ah[i];

    stbmmProd(Nm1, Nm1Ah);		// Nm1Ah = Nm1 Ah


    coeff=0.0;                              //
    localreduce=0.0;                    //
    for(i=0; i<m; i++)                  //         
      localreduce+= Ah[i]*Nm1Ah[i];         //
    MPI_Allreduce(&localreduce, &coeff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    ro=0.0;				//
    localreduce=0.0;			//
    for(i=0; i<m; i++)			//         
      localreduce+= _g[i]*Ah[i];  	//
    MPI_Allreduce(&localreduce, &ro, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

    ro = -ro/coeff ; 


    for(j=0; j<n; j++)			//  x = x + ro*h 
      x[j] = x[j] + ro*h[j];		//

    for(i=0; i<m; i++)			//  _g = _g + ro * Nm1Ah
      _g[i] = _g[i] + ro * Nm1Ah[i] ;	//

    TrMatVecProd(&A, _g, g, 0);


    for(j=0; j<n; j++)                  // 
      Cg[j]=c[j]*g[j];                       //  Cg = C g  with C = Id

    MatVecProd(&A, Cg, ACg, 0);         //  ACg = A Cg

    g2p=g2;                               // g2p = "res"                          
    localreduce=0.0;
    for(i=0; i<m; i++)                    // g2 = (ACg , _g)    
      localreduce+= ACg[i] * _g[i];


    MPI_Allreduce(&localreduce, &g2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


   t=MPI_Wtime();                       //stop timer
   if(rank==0)                          //print iterate info
     printf("k=%d res_g2=%e t=%lf \n", k, g2, t-st);
   if(g2<tol2rel){                         //
      if(rank ==0)                      //
        printf("--> converged (%e < %e) \n", g2, tol2rel);
      break;    
   }
   if(g2>g2p){                         //
      if(rank ==0)                      //
        printf("--> g2>g2p pb (%e < %e) \n", g2, tol2rel);
//      break;                        //
   }
  //  printf(" x[0]=%lf\n", x[0]);
   st=MPI_Wtime();


  gamma = g2/g2p ;

  for(j=0; j<n; j++)			// h = h * gamma - Cg
    h[j] = h[j] * gamma - Cg[j] ; 


  }  //End loop


  if(k==K){				//check unconverged
    if(rank==0)
      printf("--> unconverged, max iterate reached (%lf > %lf)\n", g2, tol2rel);
  }


      if(rank ==0)                      //
        printf("--> res=%e  \n", g2);

  free(_g);
  free(h);				
  free(Ah);
  free(Nm1Ah);


  return 1;
} 


