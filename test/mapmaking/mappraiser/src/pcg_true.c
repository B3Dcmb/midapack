// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// PCG routine applied to the mapmaking equation
// This can use the diagonal or the block-diagonal jacobi preconditionners

/** @file   pcg_true.c
    @author Frederic Dauvergne
    @date   November 2012
    @Last_update February 2019 by Hamza El Bouhargani*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include "midapack.h"
#include "mappraiser.h"



int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz Nm1, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int K)
{
  int 		i, j, k ;			// some indexes
  int		m, n, rank, size;
  double 	localreduce;			//reduce buffer
  double	st, t;				//timers
  double	solve_time;
  double	res, res0, *res_rel;
  double *tmp;
  FILE *fp;

  res_rel = (double *) malloc(1*sizeof(double));

  m=A->m;					//number of local time samples
  n=A->lcount;					//number of local pixels
  MPI_Comm_rank(A->comm, &rank);			//
  MPI_Comm_size(A->comm, &size);			//


  double *_g, *ACg, *Ah, *Nm1Ah ;		// time domain vectors
  double *g, *Cg, *h ;		      		// map domain vectors
  double *AtNm1Ah ;                    		// map domain
  double ro, gamma, coeff ;			// scalars
  double g2pix, g2pixp;
  double norm2b;

/*
  printf("n=%d, m=%d, A->nnz=%d \n", n, m, A->nnz );
*/

//Init CG descent
  // double *c;
  // c = (double *) malloc(n*sizeof(double));
  Mat BJ;


// for no preconditionner:
 // for(j=0; j<n; j++)                    //
 //   c[j]=1.;

  st=MPI_Wtime();
  // precondjacobilike( A, Nm1, lhits, cond, c);
 // precondjacobilike_avg( A, Nm1, c);
 // Compute preconditioner and process degenerate pixels
  precondblockjacobilike(A, Nm1, &BJ, b, cond, lhits);
// Redefine number of pixels in the map
  n=A->lcount-(A->nnz)*(A->trash_pix);
// Reallocate memory for well-conditioned map
  tmp = realloc(x, n * sizeof(double));
  if(tmp !=NULL){
    x = tmp;
  }


  //map domain
  h = (double *) malloc(n * sizeof(double));      //descent direction
  g = (double *) malloc(n * sizeof(double));
  AtNm1Ah = (double *) malloc(n * sizeof(double));

  //time domain
  Ah = (double *) malloc(m * sizeof(double));

  _g = Ah;
  Cg = AtNm1Ah;
  Nm1Ah = Ah;

  double *pixpond;
  pixpond = (double *) malloc(n * sizeof(double));

  //compute pixel share ponderation
  get_pixshare_pond(A, pixpond);

  //if we want to use the true norm to compute the residual
  int TRUE_NORM=1;  //0: No ; 1: Yes


  t=MPI_Wtime();

   if(rank==0)
 printf("Init PCG   t=%lf \n", t-st);
 fflush(stdout);

  st=MPI_Wtime();

  MatVecProd(A, x, _g, 0);		//
  // for(i=0; i<50; i++){//
  //     printf("MatVecProd: _g[%d] = %f\n",i,_g[i]);
  // }

  for(i=0; i<m; i++)	//
    _g[i] = b[i] + noise[i] - _g[i];		//
  // for(i=0; i<50; i++){
  //   printf("b-_g:_g[%d] = %f\n",i,_g[i]);
  // }

  stbmmProd(Nm1, _g);		// _g = Nm1 (b-Ax)
  // for(i=0; i<50; i++){//
  //     printf("Nm1*_g: _g[%d] = %f\n",i,_g[i]);
  // }

  TrMatVecProd(A, _g, g, 0);		//  g = At _g
  // for(i=0; i<50; i++){			//
  //     printf("At*_g: index = %d, g[%d] = %.18f\n", A->lindices[i], i, g[i]);
  // }

  MatVecProd(&BJ, g, Cg, 0);
  // for(j=0; j<n; j++)                    //
  //   Cg[j]=c[j]*g[j]; 			//  Cg = C g  with C = Id
  // for(i=3360; i<3380; i++){
  //     printf("index = %d , Cg[%d]=%.18f\n", A->lindices[i], i, Cg[i]);
  // }


  for(j=0; j<n; j++)   		        //  h = -Cg
    h[j]=Cg[j];


  g2pix=0.0;                               //g2 = "res"
  localreduce=0.0;
  for(i=0; i<n; i++)                    //  g2 = (Cg , g)
    localreduce+= Cg[i] * g[i] * pixpond[i];

  MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  t=MPI_Wtime();

  solve_time += (t-st);
//Just to check with the true norm:
  if (TRUE_NORM==1) {
  res=0.0;                               //g2 = "res"
  localreduce=0.0;
  for(i=0; i<n; i++)                    //  g2 = (Cg , g)
    localreduce+= g[i] * g[i] * pixpond[i];

  MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else {

  res = g2pix;

  }//end if

  double g2pixB = g2pix;
  double tol2rel = tol*tol*res;//tol*tol*g2pixB;//*g2pixB;//tol;//*tol*g2;
  res0=res;
//Test if already converged
   if(rank==0) {
     *res_rel = sqrt(res)/sqrt(res0);
     printf("res=%e \n", res);
     printf("k=%d res_g2pix=%e res_g2pix_rel=%e res_rel=%e t=%lf\n", 0, g2pix , sqrt(g2pix)/sqrt(g2pixB), sqrt(res)/sqrt(res0), t-st);
     char filename[256];
     sprintf(filename,"%s/pcg_residuals_%s.dat",outpath, ref);
     fp=fopen(filename, "wb");
     fwrite(res_rel, sizeof(double), 1, fp);
   }


   if(res<tol){                         //
      if(rank==0)                      //
        printf("--> converged (%e < %e)\n", res, tol);
        k=K;//to not enter inside the loop
   }

  st=MPI_Wtime();

  fflush(stdout);

// PCG Descent Loop
  for(k=1; k<K ; k++){


    MatVecProd(A, h, Ah, 0);		// Ah = A h
    // for(i=0; i<8; i++){//
    //     printf("MatVecProd: Ah[%d] = %f\n",i,Ah[i]);
    // }
    stbmmProd(Nm1, Nm1Ah);		// Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)
    // for(i=0; i<8; i++){//
    //     printf("Nm1Ah: Nm1Ah[%d] = %f\n",i,Nm1Ah[i]);
    // }
    TrMatVecProd(A, Nm1Ah, AtNm1Ah, 0); //  AtNm1Ah = At Nm1Ah
    // for(i=n-1; i>n-9; i--){//
    //     printf("lhs: AtNm1Ah[%d] = %f\n",i,AtNm1Ah[i]);
    // }
    coeff=0.0;
    localreduce=0.0;
    for(i=0; i<n; i++)
      localreduce+= h[i]*AtNm1Ah[i] * pixpond[i] ;
    MPI_Allreduce(&localreduce, &coeff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // printf("Coeff = %f\n",coeff);
/*
    ro=0.0;
    localreduce=0.0;
    for(i=0; i<n; i++)
      localreduce+= g[i]*Cg[i] *pixpond[i] ;
    MPI_Allreduce(&localreduce, &ro, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
*/
    ro = g2pix;

    ro = ro/coeff;
    // printf("ro = %f\n",ro);

    for(j=0; j<n; j++)			//  x = x + ro*h
      x[j] = x[j] + ro*h[j];		//


    for(j=0; j<n; j++)                  //   g = g + ro * (At Nm1 A) h
      g[j] = g[j] - ro * AtNm1Ah[j] ;


    MatVecProd(&BJ, g, Cg, 0);
    // for(j=0; j<n; j++)                  //
    //   Cg[j]=c[j]*g[j];                       //  Cg = C g  with C = Id
    // for(i=n-1; i>n-9; i--){//
    //     printf("g[%d] = %.34f\n",i,g[i]);
    //     printf("Cg[%d] = %.34f\n",i,Cg[i]);
    // }


    g2pixp=g2pix;                               // g2p = "res"
    localreduce=0.0;
    for(i=0; i<n; i++)                    // g2 = (Cg , g)
      localreduce+= Cg[i] * g[i] *pixpond[i];

    MPI_Allreduce(&localreduce, &g2pix, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   t=MPI_Wtime();                       //stop timer
   solve_time += (t-st);

//Just to check with the true norm:
  if (TRUE_NORM==1) {
  localreduce=0.0;
  for(i=0; i<n; i++)                    // g2 = (Cg , g)
    localreduce+= g[i] * g[i] *pixpond[i];

  MPI_Allreduce(&localreduce, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  else {

  res = g2pix;

  }//end if


   if(rank==0){                          //print iterate info
     *res_rel = sqrt(res)/sqrt(res0);
     printf("k=%d res_g2pix=%e res_g2pix_rel=%e res_rel=%e t=%lf \n", k, g2pix, sqrt(g2pix)/sqrt(g2pixB), sqrt(res)/sqrt(res0), t-st);
     fwrite(res_rel, sizeof(double), 1, fp);
   }
//   if(g2pix<tol2rel){                         //
fflush(stdout);

   if(res<tol2rel){
      if(rank ==0) {                     //
        printf("--> converged (%e < %e) \n", res, tol2rel);
        printf("--> i.e. \t (%e < %e) \n", sqrt(res/res0), tol);
        printf("--> solve time = %lf \n", solve_time);
        fclose(fp);
      }
      break;
   }
   if(g2pix>g2pixp){                         //
      if(rank ==0)                      //
        printf("--> g2pix>g2pixp pb (%e > %e) \n", g2pix, g2pixp);
  //    break;                        //
   }

   st=MPI_Wtime();


  gamma = g2pix/g2pixp ;

  for(j=0; j<n; j++)			// h = h * gamma - Cg
    h[j] = h[j] * gamma + Cg[j] ;


  }  //End loop


  if(k==K){				//check unconverged
    if(rank==0){
      printf("--> unconverged, max iterate reached (%lf > %lf)\n", g2pix, tol2rel);
      fclose(fp);
    }
  }

      if(rank ==0)
        printf("--> res_g2pix=%e  \n", g2pix);


  free(h);
  free(Ah);
  free(g);
  free(AtNm1Ah);
  free(res_rel);


  return 0;
}
