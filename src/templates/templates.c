/**
 * @file  templates.c
 * @version dev
 * @authors Hamza El Bouhargani
 * @date  October 2019
 * @credit  ANR-B3DCMB
 * @description Defining algebra routines for the filtering operations
 */

/*****************************************************************************/
/*                                  INCLUDE                                  */
/*****************************************************************************/
/* Templates Classes */
#include "templates.h"
#include <math.h>

/******************************************************************************/
/*                          Defining Algebra routines                         */
/******************************************************************************/
/* Projecting templates amplitudes in time domain */
int TVecProd(TemplateClass *X, int m, double sampling_freq, int *sweeptstams,
  double *tau, double *out){

  int i,j,k;
  for(k=0;k<m;k++){
    if((X+k)->flag_w == "OTF"){
      (X+k)->bins = (int *) calloc((X+k)->nsamples * (X+k)->nmult, sizeof(int));
      (X+k)->wghts = (double *) calloc((X+k)->nsamples * (X+k)->numlt, sizeof(double));
      expandpolydata(X+k, (X+k)->bins, (X+k)->wghts, sweeptstamps, sampling_freq, X->order);
    }
    for(i=0;i<(X+k)->nsamples;i++){
      for(j=0;j<(X+k)->nmult;j++){
        out[(X+k)->tinit + i] += (X+k)->wghts[i * (X+k)->nmult + j] * tau[(X+k)->bins[i * (X+k)->nmult + j]];
      }
    }
    if((X+k)->flag_w == "OTF"){
      free((X+k)->bins);
      free((X+k)->wghts);
    }
  }
}

/* Projecting time domain in templates space */
int TrTVecProd(TemplateClass *X, int m, double sampling_freq, int *sweeptstamps,
  double *d, double *out){

  int i,j,k;
  for(k=0;k<m;k++){
    if((X+k)->flag_w == "OTF"){
      (X+k)->bins = (int *) calloc((X+k)->nsamples * (X+k)->nmult, sizeof(int));
      (X+k)->wghts = (double *) calloc((X+k)->nsamples * (X+k)->numlt, sizeof(double));
      expandpolydata(X+k, (X+k)->bins, (X+k)->wghts, sweeptstamps, sampling_freq, X->order);
    }
    for(i=0;i<(X+k)->nsamples;i++){
      for(j=0;j<(X+k)->nmult;j++){
        out[(X+k)->bins[i * (X+k)->nmult + j]] += (X+k)->wghts[i * (X+k)->nmult + j] * d[(X+k)->tinit + i];
      }
    }
    if((X+k)->flag_w == "OTF"){
      free((X+k)->bins);
      free((X+k)->wghts);
    }
  }
}

/* Building Kernel Blocks */
int BuildKernel(TemplateClass *X, int n, double *B){
  // n : number of template classes in one kernel Block (1 det data in 1 CES)
  int l,k;
  for(k=0;k<n;k++){
    if((X+k)->flag_w == "OTF"){
      (X+k)->bins = (int *) calloc((X+k)->nsamples * (X+k)->nmult, sizeof(int));
      (X+k)->wghts = (double *) calloc((X+k)->nsamples * (X+k)->numlt, sizeof(double));
      expandpolydata(X+k, (X+k)->bins, (X+k)->wghts, sweeptstamps, sampling_freq, X->order);
    }
    for(l=0;l<n;l++){
      if((X+l)->flag_w == "OTF"){
        (X+l)->bins = (int *) calloc((X+k)->nsamples * (X+l)->nmult, sizeof(int));
        (X+l)->wghts = (double *) calloc((X+k)->nsamples * (X+l)->numlt, sizeof(double));
        expandpolydata(X+l, (X+l)->bins, (X+l)->wghts, sweeptstamps, sampling_freq, X->order);
      }
      for(i=0;i<(X+k)->nsamples * (X+k)->nmult;i++){
        for(j=0;j<(X+l)->nsamples * (X+l)->nmult;j++){
          B[((X+k)->bins[i]-(X+k)->nbinMin) * (X+l)->nbins + ((X+l)->bins[j]-(X+l)->nbinMin)] += (X+k)->wghts[i] * (X+l)->wghts[j];
        }
      }
      if((X+l)->flag_w == "OTF"){
        free((X+l)->bins);
        free((X+l)->wghts);
      }
    }
    if((X+k)->flag_w == "OTF"){
      free((X+k)->bins);
      free((X+k)->wghts);
    }
  }
}

/* Inverting the kernel Blocks */
/* Description:
  Invert kernel block by minimizing || B*X - I ||_F , where ||.||_F is the Frobenius norm
  The least square problem is solved using SVD, and gives the Moore-Penrose inverse by design
  The singular values (here also eigenvalues) of the kernel are also given in the vector s
  The rank of the Kernel Block is given in 'rank'
*/
/* Parameters:
  - Input:
    @B: the input kernel block
    @Binv: the pseudo-inverse of the kernel block
  - Output:
    @rank: returns the (numerical) rank of the kernel block
*/
int InvKernel(double *B, int n, double *Binv){
  int i,info,rank;
  double epsilon = -1; // machine epsilon: threshold for regularization
  double *s; // vector of singular values (from SVD)
  // Allocate memory and initialize Binv as unit matrix
  s = (double *) calloc(n, sizeof(double));
  for(i=0;i<n;i++)  Binv[i+n*i] = 1.0;
  // Call LAPACK routine to SVD solve the least square problem (Cf. Description)
  info = LAPACKE_dgelss(LAPACK_ROW_MAJOR, n, n, n, B, n, Binv, n, s, epsilon, &rank);
  if (info!=0){
    printf("Failure of the pseudo-inverse computation, LAPACK exit index = %d\n",info);
    exit(1);
  }
  return rank;
}
/* Define Legendre Polynomials */
double P0(double x){
    return 1;
}

double P1(double x){
    return x;
}
//The following is a general functoin that returns the value of the Legendre Polynomial for any given x and n=0,1,2,3,...
double Pn(double x, int n){
    if(n==0){
        return P0(x);
    }else if(n==1){
        return P1(x);
    }else{
        return (double)((2*n-1)*x*Pn(x,n-1)-(n-1)*Pn(x,n-2))/n;
    }
}
// Orthogonal Legendre in [a,b] interval
double Legendre(double x, double a, double b, int n){
  return Pn((2*x-a-b)/(b-a),n);
}

/******************************************************************************/
/*        Utility routines for building templates classes objects             */
/******************************************************************************/
/* Build local templates classes list */
/* Parameters:
  - Input:
    @X: The templates classes list
    @ndet: local number of detectors assigned to the process
    @npoly: maximum ordered considered for polynomial templates
    N.B: Will introduce nground, nhwp and ncustom when we implement support for
    other templates classes types. Flags are ignored for now.
*/
int Tlist_init(TemplateClass *X, int ndet, int *detnsamples, int *detnsweeps,
  int *sweeptstamps, double sampling_freq, int npoly){

  int i,j;
  int tinit,tlast,nbinMin;
  // looping over local list of detectors
  tinit = 0;
  tlast = detnsamples[0] -1;
  nbinMin = 0;
  nbinMax = detnsweeps[0]-1;
  for(i=0;i<ndet;i++){
    // looping over polynomial templates for each detector
    for(j=0;j<npoly;j++){
      // looping over polynomial orders
      // For each order build the corresponding polynomial template class
      Polyinit((X + (i*npoly) + j), tinit, tlast, nbinMin, nbinMax, sweeptstamps, sampling_freq, j);
      nbinMin += detnsweeps[i];
      nbinMax += detnsweeps[i];
    }
    tinit += detnsamples[i];
    tlast += detnsamples[i];
  }
}
/* Build polynomial template class instance */
int Polyinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  int *sweeptstamps, double sampling_freq, int order){

  char flag[256];
  char flag_w[256];
  int *bins;
  double *wghts;
  sprintf(flag,"NULL");
  sprinf(flag_w,"OTF");
  TCinit(X, tinit, tlast, nbinMin, nbinMax, 1, flag, flag, flag ,flag_w, flag);
  X->order = order;
  if(strcmp(X->flag_w,"S")==0){ // Stored weights
    bins = (int *) calloc(X->nsamples * X->nmult,sizeof(int));
    wghts = (double *) calloc(X->nsamples * X->nmult, sizeof(double));
    expandpolydata(X, bins, wghts, sweeptstamps, sampling_freq, order);
  }
}
/* Build Templates Classes objects */
int TCinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  int nmult, char *flag_det, char *flag_CES,
  char* flag_dataset, char *flag_w, char *ID){

  X->tinit = tinit;
  X->tlast = tlast;
  X->nsamples = tlast - tinit + 1;
  X->nbinMin = nbinMin;
  X->nbinMax = nbinMax;
  X->nbins = nbinMax - nbinMin + 1;
  X->nmult = nmult;
  X->flag_det = flag_det;
  X->flag_CES = flag_CES;
  X->flag_dataset = flag_dataset;
  X->flag_w = flag_w;
  X->ID = ID;
}
/* Build bins and weights profile for a polynomial template nmult = 1 for now */
int expandpolydata(TemplateClass *X ,int *bins, double *wghts,
  int *sweepstamps, double sampling_freq, int order){

  int i, j;
  for(i=0;i<X->nbins;i++){
    for(j=sweeptstamps[i];j<sweeptstamps[i+1];j++){
      bins[j] = X->nbinMin + i;
      wghts[j] = Legendre(j*sampling_freq, sweepstamps[i]*sampling_freq, sweepstamps[i+1]*sampling_freq, order);
    }
  }
}
