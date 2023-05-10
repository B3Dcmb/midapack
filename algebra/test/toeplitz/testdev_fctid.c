#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "toeplitz.h"
#include <time.h>


char CHAR_RW='\0';  //global variable for write mode

int main (int argc, char *argv[])
{
  int test=0;  //index of the choosen test
  int nb_test=1;  //number of example tests

  if (argc>1 && atoi(argv[1])>0) {
    test=atoi(argv[1]);
  }



  int i;
  int n, m, id0, l, lambda, blocksize, nfft;
  double *Vrshp;
  int lrshp, nrshp, mrshp;
  int i_out;
  int rfirst=id0%n;

  int i_out0;

 int nbcol;
 int *nocol;

  FILE *file;
  file = stdout;



  if (test==1) {
  fprintf(file, "fctid_mat2vect\n");
  fprintf(file, "---------------\n");

  id0=3;
  n=10;
  m=2;
  l=16;
  lambda=3;
  rfirst=id0%n;

  nrshp=l+(lambda-1)*(m-1);//(n+lambda-1)*(m-1) + n;
  mrshp=1;
  lrshp=l+(lambda-1)*(m-1);


  for(i=0;i<lrshp;i++) {
    i_out = fctid_mat2vect(i, rfirst, n, lambda);
    fprintf(file, "i=%d , i_out=%d\n", i, i_out);
  }

}//end test1
//---------------------------------------------------

else if (test==2) {
  fprintf(file, "fctid_concatcol\n");
  fprintf(file, "---------------\n");

  id0=3;
  n=10;
  m=3;
  l=27;
  lambda=3;

  nbcol=2;

  nrshp=n;
  mrshp=nbcol;
  lrshp=nrshp*mrshp;

 nocol = (int *) calloc(nbcol, sizeof(double)); 

 nocol[0]=0;
 nocol[1]=1;

  for(i=0;i<lrshp;i++) {
    i_out = fctid_concatcol(i, id0, n, l, lambda, nocol, nbcol);
    fprintf(file, "i=%d , i_out=%d\n", i, i_out);
  }


}//end test2


//---------------------------------------------------

else if (test==3) {
  fprintf(file, "fctid_mat2vect --> fctid_concatcol\n");
  fprintf(file, "----------------------------------\n");

  id0=0;
  n=10;
  m=3;
  l=30;
  lambda=3;

  nbcol=2;

  nrshp=(n+lambda-1)*(nbcol-1) + n;
  mrshp=1;
  lrshp=nrshp*mrshp;

  fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);

 nocol = (int *) calloc(nbcol, sizeof(double));

 nocol[0]=1;
 nocol[1]=2;

  int i_out2;

  for(i=0;i<lrshp;i++) {
//    i_out0=fctid_concatcol(i, id0, n, lambda, nocol, nbcol);
    i_out = fctid_mat2vect(i , rfirst, n, lambda);
    i_out2=fctid_concatcol(i_out, id0, n, lambda, nocol, nbcol);
    fprintf(file, "i=%d , i_out=%d, i_out2=%d\n", i, i_out, i_out2);
  }


}//end test3



//---------------------------------------------------
else if (test==4) {
  fprintf(file, "fctid_vect2nfftblock\n");
  fprintf(file, "---------------\n");

  id0=0;
  n=10;
  m=1;
  l=10;
  lambda=3;

  int j;
  int v1_size=9;
  int nfft=2;
  int distcorrmin=lambda-1;
  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;

  nrshp=fft_size;
  mrshp=nfft;
  lrshp=nrshp*mrshp;

  fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);
  fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);

  for(j=0;j<mrshp;j++) {
  for(i=0;i<nrshp;i++) {
    i_out = fctid_vect2nfftblock( i+j*nrshp, v1_size, fft_size, nfft, lambda);
    fprintf(file, "i=%d , i_out=%d\n", i+j*nrshp, i_out);
  }
  fprintf(file, "--\n");
  }


}//end test4

//---------------------------------------------------
else if (test==5) {
  fprintf(file, "fctid_vect2nfftblock-->vect2mat\n");
  fprintf(file, "---------------\n");

  id0=0;
  n=3;
  m=2;
  l=6;
  lambda=3;
  
  int j;
  int distcorrmin=lambda-1;
  int m_eff=m;
  int v1_size = l + (m_eff-1)*distcorrmin;
  int nfft=2;
  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;
  int rfirst=id0%n;

  nrshp=fft_size;
  mrshp=nfft;
  lrshp=nrshp*mrshp;

  fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);
  fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);
  fprintf(file, "n=%d , m=%d\n", n, m);

  int i_out2;

  for(j=0;j<mrshp;j++) {
  for(i=0;i<nrshp;i++) {
    i_out = fctid_vect2nfftblock( i+j*nrshp, v1_size, fft_size, nfft, lambda);
    i_out2 = fctid_mat2vect(i_out , rfirst, n, lambda);
    fprintf(file, "i=%d \t i_out=%d \t i_out2=%d\n", i+j*nrshp, i_out, i_out2);
  }
  fprintf(file, "--\n");
  }



}//end test5



//---------------------------------------------------
else if (test==6) {

  nbcol=2;
  nocol = (int *) calloc(nbcol, sizeof(double));
  nocol[0]=0;
  nocol[1]=2;

  int res_is_needconcat;
  res_is_needconcat=is_needconcat(nocol, nbcol);

  fprintf(file, "res_is_needconcat=%d\n", res_is_needconcat);

  nbcol=2;
  nocol = (int *) calloc(nbcol, sizeof(double));
  nocol[0]=2;
  nocol[1]=3;

  res_is_needconcat=is_needconcat(nocol, nbcol);

  fprintf(file, "res_is_needconcat=%d\n", res_is_needconcat);


  nbcol=3;
  nocol = (int *) calloc(nbcol, sizeof(double));
  nocol[0]=2;
  nocol[1]=3;
  nocol[2]=4;

  res_is_needconcat=is_needconcat(nocol, nbcol);
  fprintf(file, "res_is_needconcat=%d\n", res_is_needconcat);


  nbcol=3;
  nocol = (int *) calloc(nbcol, sizeof(double));
  nocol[0]=2;
  nocol[1]=3;
  nocol[2]=5;

  res_is_needconcat=is_needconcat(nocol, nbcol);
  fprintf(file, "res_is_needconcat=%d\n", res_is_needconcat);


  srand (time (NULL));  //init seed
  for(i=0;i<20;i++) {
  int nombre_aleatoire = rand_a_b(0,5);// b borne non incluse
  printf("%d ",nombre_aleatoire);
  }
  printf("\n");



}//end test6


//---------------------------------------------------
else if (test==7) {
  fprintf(file, "fctid_mat2vect\n");
  fprintf(file, "---------------\n");

  id0=3;
  n=10;
  m=2;
  l=16;
  lambda=3;
  rfirst=id0%n;

  nrshp=l+(lambda-1)*(m-1);//(n+lambda-1)*(m-1) + n;
  mrshp=1;
  lrshp=l+(lambda-1)*(m-1);


  i=-1;
    i_out = fctid_mat2vect(i, rfirst, n, lambda);
    fprintf(file, "i=%d , i_out=%d\n", i, i_out);



}//end test7



else if (test==8) {
  fprintf(file, "fctid_mat2vect_inv\n");
  fprintf(file, "---------------\n");

  id0=1;
  n=10;
  m=2;
  l=17;
  lambda=3;
  rfirst=id0%n;

  nrshp=l+(lambda-1)*(m-1);//(n+lambda-1)*(m-1) + n;
  mrshp=1;
  lrshp=l+(lambda-1)*(m-1);


  for(i=0;i<l;i++) {
    i_out = fctid_mat2vect_inv(i, id0, n, lambda);
    fprintf(file, "i=%d , i_out=%d\n", i, i_out);
  }

}//end test8

//---------------------------------------------------

else if (test==9) {
  fprintf(file, "fctid_concatcol\n");
  fprintf(file, "---------------\n");

  id0=3;
  n=10;
  m=3;
  l=27;
  lambda=3;

  nbcol=2;

  nrshp=n;
  mrshp=nbcol;
  lrshp=nrshp*mrshp;

 nocol = (int *) calloc(nbcol, sizeof(double));

 nocol[0]=0;
 nocol[1]=2;

 int *nocol_inv;
  nocol_inv = (int *) calloc(m, sizeof(double));

  for(i=0;i<m;i++)
    nocol_inv[i]=-1;
  for(i=0;i<nbcol;i++)
    nocol_inv[nocol[i]]=i;

  for(i=0;i<m;i++)
    printf("nocol_inv[%d]=%d \t", i, nocol_inv[i]);

  printf("\n");


  for(i=0;i<l+1;i++) {
    i_out = fctid_concatcol_inv(i, id0, n, m, l, lambda, nocol_inv, nbcol);
    fprintf(file, "i=%d , i_out=%d\n", i, i_out);
  }


}//end test9


//---------------------------------------------------

else if (test==10) {
  fprintf(file, "fctid_concatcol\n");
  fprintf(file, "---------------\n");

  id0=0;
  n=10;
  m=3;
  l=30;
  lambda=3;

  nbcol=2;

  nrshp=n;
  mrshp=nbcol;
  lrshp=nrshp*mrshp;

 nocol = (int *) calloc(nbcol, sizeof(double));

 nocol[0]=0;
 nocol[1]=2;

 int *nocol_inv;
  nocol_inv = (int *) calloc(m, sizeof(double));

  for(i=0;i<m;i++)
    nocol_inv[i]=-1;
  for(i=0;i<nbcol;i++)
    nocol_inv[nocol[i]]=i;

  for(i=0;i<m;i++)
    printf("nocol_inv[%d]=%d \t", i, nocol_inv[i]);

  printf("\n");

  int i_out1;

  int idf = id0+l-1;
  int lconc = n*nbcol - (nocol[0]==0)*(id0%n) - (nocol[nbcol]==m)*(n-idf%n);

  for(i=0;i<lconc+1;i++) {
    i_out1 = fctid_concatcol(i, id0, n, m, l, lconc, lambda, nocol, nbcol);
    i_out = fctid_concatcol_inv(i_out1, id0, n, m, l, lconc, lambda, nocol_inv, nbcol);
    fprintf(file, "i=%d , i_out1=%d, i_out=%d\n", i, i_out1, i_out);
  }


}//end test10


//---------------------------------------------------
else if (test==11) {
  fprintf(file, "fctid_vect2nfftblock\n");
  fprintf(file, "---------------\n");

  id0=0;
  n=10;
  m=1;
  l=10;
  lambda=3;

  int j;
  int nfft=2;
  int distcorrmin=lambda-1;

  int idf = id0+l-1;
//  int lconc = n*nbcol - (nocol[0]==0)*(id0%n) - (nocol[nbcol]==m)*(n-idf%n);
  int v1_size=l;//lconc+(distcorrmin)*(nbcol-1);
  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;

  nrshp=fft_size;
  mrshp=nfft;
  lrshp=nrshp*mrshp;

  fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);
  fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);

  for(i=0;i<v1_size;i++) {
    i_out = fctid_vect2nfftblock_inv( i, v1_size, fft_size, nfft, lambda);
    fprintf(file, "i=%d , i_out=%d\n", i+j*nrshp, i_out);
  }


}//end test11


return 0;

}


int rand_a_b(int a, int b){
    return rand()%(b-a) +a;
}



