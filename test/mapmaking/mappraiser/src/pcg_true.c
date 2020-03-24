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
  // for(i=0; i<10; i++){//
  //     printf("MatVecProd: _g[%d] = %f\n",i,_g[i]);
  // }

  for(i=0; i<m; i++)	//
    _g[i] = b[i] + noise[i] - _g[i];		//
  // for(i=0; i<10; i++){
  //   printf("b-_g:_g[%d] = %f\n",i,_g[i]);
  // }

  stbmmProd(Nm1, _g);		// _g = Nm1 (b-Ax)
  // for(i=0; i<10; i++){//
  //     printf("Nm1*_g: _g[%d] = %f\n",i,_g[i]);
  // }

  TrMatVecProd(A, _g, g, 0);		//  g = At _g
  // for(i=0; i<10; i++){			//
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

int PCG_GLS_templates(char *outpath, char *ref, Mat *A, Tpltz Nm1, TemplateClass *X, double *B, int **sweeptstamps, int npoly, int *nsweeps, int **az_binned, int n_sss_bins, int nces, int nb_blocks_loc, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int K, double sampling_freq)
{
  int 		i, j, k, l, ces_id ;			// some indexes
  int		m, n, rank, size;
  double 	localreduce;			//reduce buffer
  double	st, t,t2,st2;				//timers
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

  int ndet = nb_blocks_loc / nces;

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

  //templates domain
  int nb_templates_loc = 0;
  for(i=0;i<nces;i++)
    nb_templates_loc += ndet * (npoly*nsweeps[i] + n_sss_bins);
  double *out1 = (double *) malloc(nb_templates_loc * sizeof(double));
  double *out2 = (double *) malloc(nb_templates_loc * sizeof(double));

  //time domain
  Ah = (double *) malloc(m * sizeof(double));
  double *Tvec= (double *) malloc(m * sizeof(double));

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

  // to test that we get zero residual for template modes
  // int cumul_offset = 0;
  // for(ces_id=0;ces_id<nces;ces_id++){
  //   for(i=0;i<(npoly*nsweeps[ces_id]+n_sss_bins)*ndet;i++)
  //     out2[cumul_offset+i] = 1+0.01*(i%n_sss_bins);
  //   cumul_offset += (npoly*nsweeps[ces_id]+n_sss_bins)*ndet;
  // }
  // TVecProd(X, nces, (npoly+1) * ndet, sampling_freq, sweeptstamps, az_binned, out2, noise);
  // for(i=0;i<10;i++)
  //   printf("b[%d] = %f\n",i,b[i]);
  MatVecProd(A, x, _g, 0);		//
  // for(i=0; i<10; i++){//
  //     printf("MatVecProd: _g[%d] = %f\n",i,_g[i]);
  // }

  for(i=0; i<m; i++)//
    _g[i] = b[i] + noise[i] - _g[i];
		//
  // for(i=0; i<10; i++){
  //   printf("b-_g:_g[%d] = %f\n",i,_g[i]);
  // }
  // if(rank==0){
  //   char det0filename[256];
  //   char detnfilename[256];
  //   FILE *fp0,*fp1;
  //   sprintf(det0filename,"%s/det0_data_init.dat",outpath);
  //   sprintf(detnfilename,"%s/det%d_data_init.dat",outpath,Nm1.nb_blocks_loc-1);
  //   fp0=fopen(det0filename, "wb");
  //   fp1=fopen(detnfilename, "wb");
  //   fwrite(_g, sizeof(double), Nm1.tpltzblocks[0].n, fp0);
  //   fwrite(_g+Nm1.tpltzblocks[Nm1.nb_blocks_loc-1].idv, sizeof(double), Nm1.tpltzblocks[Nm1.nb_blocks_loc-1].n, fp1);
  //   fclose(fp0);
  //   fclose(fp1);
  // }

  st2=MPI_Wtime();
  int t_id = 0; //time sample index in local data
  for(i=0;i<nb_blocks_loc;i++){
    for(j=0;j<Nm1.tpltzblocks[i].n;j++){
      _g[t_id+j] = Nm1.tpltzblocks[i].T_block[0] * _g[t_id+j];
    }
    t_id += Nm1.tpltzblocks[i].n;
  }
  // stbmmProd(Nm1, _g);		// _g = Nm1 (b-Ax)
  // for(i=0; i<10; i++){//
  //     printf("Nm1*_g: _g[%d] = %f\n",i,_g[i]);
  // }
  t2 = MPI_Wtime();
  if(rank==0)
    printf("Nm1 * v, t=%lf\n",t2-st2);


  st2=MPI_Wtime();
  TrTVecProd(X, nces, (npoly+1) * ndet, sampling_freq, sweeptstamps, az_binned, _g, out1);
  // for(i=78*3; i<78*3+10; i++){//
  //     printf("TrT*_g: out1[%d] = %f\n",i,out1[i]);
  // }
  t2 = MPI_Wtime();
  if(rank==0)
    printf("Tt * v, t=%lf\n",t2-st2);



  st2 = MPI_Wtime();
  for(ces_id=0;ces_id<nces;ces_id++){
    for(i=0;i<ndet;i++){
      for(j=0;j<(npoly*nsweeps[ces_id])+n_sss_bins;j++){
        out2[(ces_id*ndet+i)*(npoly*nsweeps[ces_id]+n_sss_bins) + j] = 0;
        for(k=0;k<(npoly*nsweeps[ces_id])+n_sss_bins;k++){
          out2[(ces_id*ndet+i)*(npoly*nsweeps[ces_id]+n_sss_bins) + j] += B[(ces_id*ndet+i)*(npoly*nsweeps[ces_id]+n_sss_bins)*(npoly*nsweeps[ces_id]+n_sss_bins) + j*(npoly*nsweeps[ces_id]+n_sss_bins) + k] * out1[(ces_id*ndet+i)*(npoly*nsweeps[ces_id]+n_sss_bins) + k];
        }
        // if(out2[i*(npoly*nsweeps) + j]>1+1e-15 || out2[i*(npoly*nsweeps) + j] <1-1e-15)
        //   printf("(TtMT)^-1 * out1: out2[%d] = %f\n",i*(npoly*nsweeps) + j,out2[i*(npoly*nsweeps) + j]);
      }
    }
  }
  t2 = MPI_Wtime();
  if(rank==0)
    printf("(Tt N^-1 T)^-1 * v, t=%lf\n",t2-st2);

  // for(i=0; i<10; i++){//
  //     printf("(TtMT)^-1 * out1: out2[%d] = %f\n",i,out2[i]);
  // }
  st2 = MPI_Wtime();
  TVecProd(X, nces, (npoly+1) * ndet, sampling_freq, sweeptstamps, az_binned, out2, Tvec);
  // for(i=0; i<10; i++){//
  //     printf("(T*out2: Tvec[%d] = %f\n",i,Tvec[i]);
  // }
  t2 = MPI_Wtime();
  if(rank==0)
    printf("T * v, t=%lf\n",t2-st2);

  st2 = MPI_Wtime();
  t_id = 0; //time sample index in local data
  for(i=0;i<nb_blocks_loc;i++){
    for(j=0;j<Nm1.tpltzblocks[i].n;j++){
      Tvec[t_id+j] = Nm1.tpltzblocks[i].T_block[0] * Tvec[t_id+j];
    }
    t_id += Nm1.tpltzblocks[i].n;
  }
  // stbmmProd(Nm1, Tvec);
  // for(i=0; i<10; i++){//
  //     printf("(Nm1*Tvec: Tvec[%d] = %f\n",i,Tvec[i]);
  // }
  t2 = MPI_Wtime();
  if(rank==0)
    printf("N^-1 * v, t=%lf\n",t2-st2);

  st2 = MPI_Wtime();
  for(i=0;i<m;i++){
    _g[i] = _g[i] - Tvec[i];
    // _g[i] = b[i] + noise[i] - Tvec[i];
    // if(_g[i]>1e-15) printf("_g -= Tvec: _g[%d] = %f\n",i,_g[i]);
  }
  // for(i=0; i<10; i++){//
  //     printf("_g -= Tvec: _g[%d] = %f\n",i,_g[i]);
  // }
  t2 = MPI_Wtime();
  if(rank==0)
    printf("_g-=Tvec, t=%lf\n",t2-st2);
  // if(rank==0){
  //   char fdet0filename_filtered[256];
  //   char fdetnfilename_filtered[256];
  //   char fdet0filename_fit[256];
  //   char fdetnfilename_fit[256];
  //   FILE *fp2,*fp3,*fp4,*fp5;
  //   sprintf(fdet0filename_filtered,"%s/det0_data_filtered.dat",outpath);
  //   sprintf(fdetnfilename_filtered,"%s/det%d_data_filtered.dat",outpath,Nm1.nb_blocks_loc-1);
  //   sprintf(fdet0filename_fit,"%s/det0_data_fit.dat",outpath);
  //   sprintf(fdetnfilename_fit,"%s/det%d_data_fit.dat",outpath,Nm1.nb_blocks_loc-1);
  //   fp2=fopen(fdet0filename_filtered, "wb");
  //   fp3=fopen(fdetnfilename_filtered, "wb");
  //   fp4=fopen(fdet0filename_fit, "wb");
  //   fp5=fopen(fdetnfilename_fit, "wb");
  //   fwrite(_g, sizeof(double), Nm1.tpltzblocks[0].n, fp2);
  //   fwrite(_g+Nm1.tpltzblocks[Nm1.nb_blocks_loc-1].idv, sizeof(double), Nm1.tpltzblocks[Nm1.nb_blocks_loc-1].n, fp3);
  //   fwrite(Tvec, sizeof(double), Nm1.tpltzblocks[0].n, fp4);
  //   fwrite(Tvec+Nm1.tpltzblocks[Nm1.nb_blocks_loc-1].idv, sizeof(double), Nm1.tpltzblocks[Nm1.nb_blocks_loc-1].n, fp5);
  //   fclose(fp2);
  //   fclose(fp3);
  //   fclose(fp4);
  //   fclose(fp5);
  // }


  TrMatVecProd(A, _g, g, 0);		//  g = At _g
  // for(i=0; i<10; i++){			//
  //     printf("At*_g: index = %d, g[%d] = %.18f\n", A->lindices[i], i, g[i]);
  // }

//COMMENT
  MatVecProd(&BJ, g, Cg, 0);
  // for(j=0; j<n; j++)                    //
  //   Cg[j]=g[j]; 			//  Cg = C g  with C = Id
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
  k=0;
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
  if(k<=1){
  for(k=1; k<K ; k++){


    MatVecProd(A, h, Ah, 0);		// Ah = A h
    // for(i=0; i<8; i++){//
    //     printf("MatVecProd: Ah[%d] = %f\n",i,Ah[i]);
    // }
    t_id = 0; //time sample index in local data
    for(i=0;i<nb_blocks_loc;i++){
      for(j=0;j<Nm1.tpltzblocks[i].n;j++){
        Nm1Ah[t_id+j] = Nm1.tpltzblocks[i].T_block[0] * Nm1Ah[t_id+j];
      }
      t_id += Nm1.tpltzblocks[i].n;
    }
    // stbmmProd(Nm1, Nm1Ah);		// Nm1Ah = Nm1 Ah   (Nm1Ah == Ah)
    // for(i=0; i<8; i++){//
    //     printf("Nm1Ah: Nm1Ah[%d] = %f\n",i,Nm1Ah[i]);
    // }

    TrTVecProd(X, nces, (npoly+1) * ndet, sampling_freq, sweeptstamps, az_binned, Nm1Ah, out1);

    for(ces_id=0;ces_id<nces;ces_id++){
      for(i=0;i<ndet;i++){
        for(j=0;j<(npoly*nsweeps[ces_id]+n_sss_bins);j++){
          out2[(ces_id*ndet+i)*(npoly*nsweeps[ces_id]+n_sss_bins) + j] = 0;
          for(l=0;l<(npoly*nsweeps[ces_id]+n_sss_bins);l++){
            out2[(ces_id*ndet+i)*(npoly*nsweeps[ces_id]+n_sss_bins) + j] += B[(ces_id*ndet+i)*(npoly*nsweeps[ces_id]+n_sss_bins)*(npoly*nsweeps[ces_id]+n_sss_bins) + j*(npoly*nsweeps[ces_id]+n_sss_bins) + l] * out1[(ces_id*ndet+i)*(npoly*nsweeps[ces_id]+n_sss_bins) + l];
          }
        }
      }
    }

    TVecProd(X, nces, (npoly+1) * ndet, sampling_freq, sweeptstamps, az_binned, out2, Tvec);

    t_id = 0; //time sample index in local data
    for(i=0;i<nb_blocks_loc;i++){
      for(j=0;j<Nm1.tpltzblocks[i].n;j++){
        Tvec[t_id+j] = Nm1.tpltzblocks[i].T_block[0] * Tvec[t_id+j];
      }
      t_id += Nm1.tpltzblocks[i].n;
    }
    // stbmmProd(Nm1, Tvec);

    for(i=0;i<m;i++)
      Nm1Ah[i] = Nm1Ah[i] - Tvec[i];

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
    //   Cg[j]=g[j];                       //  Cg = C g  with C = Id
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
