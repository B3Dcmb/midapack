#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "toeplitz.h"
#include <cblas.h>
#include <time.h>

// Perform unit testing of the low level routine
// Test number as first argument
// -1 (default) : all test will be run. 
// 7 : n = 500, m = 1, nfft = 1, lambda = 150 and bloc_size = n; 
// 8 : n = 500, m = 1, nfft = 1, lambda = 150 and bloc_size = 2**9;


// -1 (default) : all test will be run.
// 0 : n = 5000, m = 10, nfft = 5 and lambda = 150;
// 1 : n = 5000, m = 10, nfft = 1 and lambda = 150;
// 2 : n = lambda, m = 1, nfft = 5 and lambda = 5000;
// 3 : n = 2*lambda, m = 10, nfft = 5 and lambda = 150;
// 4 : n = 3*lambda, m = 10, nfft = 5 and lambda = 150;
// 5 : n = 5000, m = 10, nfft = 5 and lambda = 150 (use an optimal blocsize != 0);
// 6 : n = lambda, m = 1, nfft = 1 and lambda = 150 (use vector product routine);
// 7 : n = 5000, m = 1, nfft = 1 and lambda = 150 (use vector product routine);

//extern int BLOCK_SIZE;
extern int NFFT;
int INDICE;
FILE *FILEOUT;

int main (int argc, char *argv[])
{
  int test = -1;
  if (argc>1)
    test = atoi(argv[1]);
  int n, m, nfft, lambda;

  if(test==0 || test<0) {
    lambda = 3;  // toeplitz band width
    n = 15;   // row dimension
    m = 1;       // column dimension
    nfft = 1;
    printf("test 0 : n = %d, m = %d and lambda = %d\n", n, m, lambda);
    test_toeplitz(n, m, lambda, nfft, 0, 0);
  }

  if(test==1 || test<0) {
    lambda = 3;  // toeplitz band width
    n = 15;   // row dimension
    m = 2;       // column dimension
    nfft = 1;
    printf("test 0 : n = %d, m = %d and lambda = %d\n", n, m, lambda);
    test_toeplitz(n, m, lambda, nfft, 0, 0);
  }


  if(test==2 || test<0) {
    lambda = 3;  // toeplitz band width
    n = 10;   // row dimension
    m = 1;       // column dimension
    nfft = 2;
    printf("test 0 : n = %d, m = %d and lambda = %d\n", n, m, lambda);
    test_toeplitz(n, m, lambda, nfft, 0, 0);
  }


  if(test==3 || test<0) {
    lambda = 3;  // toeplitz band width
    n = 10;   // row dimension
    m = 2;       // column dimension
    nfft = 2;
    printf("test 0 : n = %d, m = %d and lambda = %d\n", n, m, lambda);
    test_toeplitz(n, m, lambda, nfft, 0, 0);
  }


  if(test==4 || test<0) {
    lambda = 3;  // toeplitz band width
    n = 10;   // row dimension
    m = 1;       // column dimension
    nfft = 3;
    printf("test 0 : n = %d, m = %d and lambda = %d\n", n, m, lambda);
    test_toeplitz(n, m, lambda, nfft, 0, 0);
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


  //fich:
  FILE* file;
  char filename [1024];
  sprintf(filename,"output%d.txt",INDICE);
  file = fopen(filename, "w");

//  FILEOUT = file;
  FILEOUT = stdout; 

 
  // alloc arrays
  time (&start);
  T  = (double *) calloc(lambda,sizeof(double));
  T2 = (double *) calloc(n*n   ,sizeof(double));
  V  = (double *) calloc(n*m   ,sizeof(double));  
  TV = (double *) calloc(n*m   ,sizeof(double));  
  TV2= (double *) calloc(n*m   ,sizeof(double));  
  
  srand (time (NULL)); // init seed
  
  // init arrays
  for(i=0;i<lambda;i++) // Toeplitz matrix (band only)
    T[i]=rand()/((double) RAND_MAX);

  T[0]=10; T[1]=2; T[2]=3;

  for (j=0;j<n;j++){ // Full Toeplitz matrix needed for cblas computation
    for(i=0;i<n;i++)
      T2[j*n+i] = 0;
    for (i=0;i<lambda;i++){
      if (j-i>=0)
  	T2[j*n+j-i] = T[i];
      if (j+i<n)
  	T2[j*n+j+i] = T[i];    }}
  for(j=0;j<m;j++) // input matrix
    for(i=0;i<n;i++)
      V[j*n+i] = (j*n+i) +1;//rand()/((double) RAND_MAX);//1.e-200;//0.0;//rand()/((double) RAND_MAX);
  time (&end);
  dif = difftime (end,start);
//  printf ("matrix initialization took:\t %.2lf s\n", dif );
  

  
  // perform product using cblas routine
  time (&start);
  cblas_dgemm (CblasColMajor, CblasNoTrans,CblasNoTrans, n, m, n, 1, T2, n, V, n, 1, TV2, n);
  time (&end);
  dif = difftime (end,start);
//  printf ("c blas routine took:\t \t %.2lf s\n", dif );

  // perform product using midas DA lib
  fftw_complex *V_fft, *T_fft;
  double *V_rfft;
  fftw_plan plan_f, plan_b;

  //BLOCK_SIZE = 0;
  NFFT=nfft;

  int blocsize0, nfft0;
  tpltz_init(n, lambda , &nfft0, &blocsize0, &T_fft, T, &V_fft, &V_rfft, &plan_f, &plan_b);


    printf("blocsize0=%d\n", blocsize0);
    printf("nfft0=%d\n", nfft0);

  time (&start);

  printf("matrix computation...\n");

  int id0=0;
  int l=n*m;

  int ngap;
  int *id0gap;
  int *lgap;

  ngap=2;
  id0gap = (int *) calloc(ngap, sizeof(int));
  lgap = (int *) calloc(ngap, sizeof(int));


  id0gap[0]=3;
  lgap[0]=5;
  id0gap[1]=10;
  lgap[1]=2;

  int k;
  //reset gaps

  for (j=0; j<m; j++)
  for (i=0 ; i<ngap; i++)
    for (k=0 ; k<lgap[i]; k++)
      V[id0gap[i]+k+j*n]=0.0;





//  get_gappy_blocks_params(n, m, nrow, T, nb_blocs_local, nb_blocs_all, lambda, idv, id0gap, lgap, ngap);




  time (&end);


/*
  tpltz_cleanup(&T_fft, &V_fft, &V_rfft,&plan_f, &plan_b );
  dif = difftime (end,start);
  printf ("toeplitz matrix product took:\t %.2lf s\n", dif );
  

  
  // compute difference
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

*/


//  }
//  else {
//    printf ("nothing to do\n");
//  }

  free(T);
  free(T2);
  free(V);
  free(TV);
  free(TV2);
  
  fclose(file);
 
  return 0;
}




