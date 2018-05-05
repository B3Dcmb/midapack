#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "toeplitz.h"
#include <cblas.h>
#include <time.h>

// Perform unit testing of the low level routine
// Test number as first argument

//extern int BLOCK_SIZE;
extern int NFFT;

int main (int argc, char *argv[])
{
  int test = -1;
  if (argc>1)
    test = atoi(argv[1]);
  int n, m, nfft, lambda;
  
  if(test==0 || test < 0){
    n = 10;      // row dimension
    lambda = 3;  // toeplitz band width
    m = 1;       // column dimension
    nfft = 1;
    printf("test 0 : n = %d, m = %d, nfft=%d and lambda = %d\n", n, m, nfft, lambda);
    test_toeplitz(n, m, lambda, nfft,0,0);
      }


  if(test==1 || test < 0){
    n = 10;      // row dimension
    lambda = 3;  // toeplitz band width
    m = 2;       // column dimension
    nfft = 1;
    printf("test 0 : n = %d, m = %d, nfft=%d and lambda = %d\n", n, m, nfft, lambda);
    test_toeplitz(n, m, lambda, nfft,0,0);
      }

  if(test==2 || test < 0){
    n = 10;      // row dimension
    lambda = 3;  // toeplitz band width
    m = 2;       // column dimension
    nfft = 2;
    printf("test 0 : n = %d, m = %d, nfft=%d and lambda = %d\n", n, m, nfft, lambda);
    test_toeplitz(n, m, lambda, nfft,0,0);
      }


  if(test==3 || test < 0){
    n = 10;      // row dimension
    lambda = 3;  // toeplitz band width
    m = 2;       // column dimension
    nfft = 3;
    printf("test 0 : n = %d, m = %d, nfft=%d and lambda = %d\n", n, m, nfft, lambda);
    test_toeplitz(n, m, lambda, nfft,0,0);
      }



  if(test==4 || test < 0){
    n = 21;      // row dimension
    lambda = 3;  // toeplitz band width
    m = 1;       // column dimension
    nfft = 1;
    printf("test 0 : n = %d, m = %d, nfft=%d and lambda = %d\n", n, m, nfft, lambda);
    test_toeplitz(n, m, lambda, nfft,0,0);
      }
}

// Compare middle level routine output to Cblas
int test_toeplitz(int n, int m, int lambda, int nfft, int blocsize, int vector)
{
  time_t start,end;
  double dif;

  double *T, *T2;  // input toeplitz
  double *V;       // input matrix
  double *TV, *TV2;// output matrix
  
  int i,j; // some index
  
  // alloc arrays
  T  = (double *) calloc(lambda,sizeof(double));
  T2 = (double *) calloc(n*n   ,sizeof(double));
  V  = (double *) calloc(n*m   ,sizeof(double));  
  TV = (double *) calloc(n*m   ,sizeof(double));  
  TV2= (double *) calloc(n*m   ,sizeof(double));  
  
  srand (time (NULL)); // init seed
  
  // init arrays
  for(i=0;i<lambda;i++) // Toeplitz matrix (band only)
    T[i]= 0.0;  //rand()/((double) RAND_MAX);

  T[0]=10; T[1]=2; T[2]=3;

  for (j=0;j<n;j++){ // Full Toeplitz matrix needed for cblas computation
    for(i=0;i<n;i++)
      T2[j*n+i] = 0;
    for (i=0;i<lambda;i++){
      if (j-i>=0)
  	T2[j*n+j-i] = T[i];
      if (j+i<n)
  	T2[j*n+j+i] = T[i];    }}

    for (i=0;i<(n*m);i++)
      V[i] = i+1;//%12;// rand()/((double) RAND_MAX);//1.e-200;//0.0;//rand()/((double) RAND_MAX);

//  V[0]=0;V[1]=0;V[2]=0;
//  V[9]=0;V[10]=0;V[11]=0;

  printf ("n=%d , m=%d, lambda=%d\n", n, m, lambda);
  printf("\n");

  for (i=0;i<(lambda);i++)
    printf("T[%d]=%f\n", i, T[i]);

  printf("\n");

  for (i=0;i<(n*m);i++)
    printf("V[%d]=%f\n", i, V[i]);

  printf("\n");
  
  // perform product using cblas routine
  cblas_dgemm (CblasColMajor, CblasNoTrans,CblasNoTrans, n, m, n, 1, T2, n, V, n, 1, TV2, n);

  // perform product using midas DA lib
  fftw_complex *V_fft, *T_fft;
  double *V_rfft;
  fftw_plan plan_f, plan_b;


  //BLOCK_SIZE = 9;
  NFFT=nfft;

  int blocsize0, nfft0;
  tpltz_init(n, lambda , &nfft0, &blocsize0, &T_fft, T, &V_fft, &V_rfft, &plan_f, &plan_b);

    printf("blocsize0=%d\n", blocsize0);
    printf("nfft0=%d\n", nfft0);


  printf("\n");
//  stmm_core(&V, n, m, T_fft, blocsize0, lambda, V_fft, V_rfft, nfft0, plan_f, plan_b, 0);
  stmm_simple_core(&V, n, m, T, blocsize0, lambda, nfft0, 0);


  printf("======================\n");

  for (i=0;i<(lambda);i++)
    printf("T[%d]=%f\n", i, T[i]);

  printf("\n");

  for (i=0;i<(n*m);i++)
    printf("V[%d]=%f\n", i, V[i]);

  printf("\n");


  tpltz_cleanup(&T_fft, &V_fft, &V_rfft,&plan_f, &plan_b );

  
  //compute difference
  double res =0;
  int nb_false = 0;
  for(j=0;j<m;j++)
    for(i=0;i<n;i++)
      res += fabs(V[j*n+i]-TV2[j*n+i])/TV2[j*n+i];
  res = res/(n*m);
  for (i=0;i<n*m;i++)
    if(fabs(V[i]-TV2[i])>1e-8) {
    //printf("TV[%d,0]=%e \t TV2[%d,0]=%e \t diff=%e\n", i, V[i], i, TV2[i], V[i]-TV2[i]);
    nb_false += 1;
    } 


  if (fabs(res)<1e-8)
    printf("Success ! (difference is about %e)\n", fabs(res));
  else {
    printf("Test failed ... (difference is about %e)\n", fabs(res));
    printf("Number of false values : %d\n",nb_false);}


//  }
//  else {
//    printf ("nothing to do\n");
//  }

  free(T);
  free(T2);
  free(V);
  free(TV);
  free(TV2);
  
  
  return 0;
}


