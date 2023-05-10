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

  FILE *file;
  file = stdout;


//---------------------------------------------------
 if (test==1) {
  fprintf(file, "fctid_vect2nfftblock-->vect2mat\n");
  fprintf(file, "---------------\n");

  id0=0;
  n=3;
  m=3;
  l=9;
  lambda=3;

//  int distcorrmin=lambda-1;
//  int m_eff=m;
//  int v1_size = l + (m_eff-1)*distcorrmin;
  int nfft=2;
//  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;
//  int rfirst=id0%n;

  int nbcol=2;
  int *nocol;
  nocol = (int *) calloc(nbcol, sizeof(double));
  nocol[0]=0;
  nocol[1]=2;

//  nrshp=fft_size;
//  mrshp=nfft;
//  lrshp=nrshp*mrshp;

 // fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);
 // fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);
  fprintf(file, "n=%d , m=%d\n", n, m);

  double *Vin;
  Vin = (double *) calloc(l, sizeof(double));

  for(i=0;i<l;i++)
    Vin[i]=i+1;

  int lconc=6;

  int distcorrmin= lambda-1;
  int v1_size=lconc+(distcorrmin)*(nbcol-1);
  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;

  int nrshp=fft_size;
  int mrshp=nfft;
  int lrshp=nrshp*mrshp;

  double *Vrshp;
  Vrshp = (double *) calloc(lrshp, sizeof(double));

  int flag_outformat=2;

  build_reshape(Vin, nocol, nbcol, lconc, n, m, id0, l, lambda, nfft, Vrshp, nrshp, mrshp, lrshp, flag_outformat);


  fprintf(file, "---------------\n");
  for(i=0;i<lrshp;i++) {
  if (i%nrshp==0)
    fprintf(file, "--\n");
  fprintf(file, "Vrshp[%d]=%f\n", i, Vrshp[i]);
  }

}//end test1



//---------------------------------------------------
 if (test==2) {
  fprintf(file, "fctid_vect2nfftblock-->vect2mat\n");
  fprintf(file, "---------------\n");

  id0=0;
  n=3;
  m=3;
  l=9;
  lambda=3;

//  int distcorrmin=lambda-1;
//  int m_eff=m;
//  int v1_size = l + (m_eff-1)*distcorrmin;
  int nfft=2;
//  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;
//  int rfirst=id0%n;

  int nbcol=2;
  int *nocol;
  nocol = (int *) calloc(nbcol, sizeof(double));
  nocol[0]=0;
  nocol[1]=2;

//  nrshp=fft_size;
//  mrshp=nfft;
//  lrshp=nrshp*mrshp;

 // fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);
 // fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);
  fprintf(file, "n=%d , m=%d\n", n, m);

  double *Vout;
  Vout = (double *) calloc(l, sizeof(double));

  int lconc=6;

  int distcorrmin= lambda-1;
  int v1_size=lconc+(distcorrmin)*(nbcol-1);
  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;

  int nrshp=fft_size;
  int mrshp=nfft;
  int lrshp=nrshp*mrshp;

  fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);
  fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);

  double *Vrshp;
  Vrshp = (double *) calloc(lrshp, sizeof(double));

  for(i=0;i<lrshp;i++)
    Vrshp[i]=i+1;

  int flag_format_rshp=2;

  extract_result(Vout, nocol, nbcol, lconc, n, m, id0, l, lambda, nfft, Vrshp, nrshp, mrshp, lrshp, flag_format_rshp);


  fprintf(file, "---------------\n");
  for(i=0;i<l;i++) {
  if (i%n==0)
    fprintf(file, "--\n");
  fprintf(file, "Vout[%d]=%f\n", i, Vout[i]);
  }

}//end test2



//---------------------------------------------------
 if (test==3) {
  fprintf(file, "fctid_vect2nfftblock-->vect2mat\n");
  fprintf(file, "---------------\n");

  id0=1;
  n=3;
  m=3;
  l=8;
  lambda=3;

  int nfft=2;
  int nbcol=2;
  int *nocol;
  nocol = (int *) calloc(nbcol, sizeof(double));
  nocol[0]=0;
  nocol[1]=2;

//  nrshp=fft_size;
//  mrshp=nfft;
//  lrshp=nrshp*mrshp;

 // fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);
 // fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);
  fprintf(file, "n=%d , m=%d\n", n, m);

  double *Vout;
  Vout = (double *) calloc(l, sizeof(double));

//  int lconc=6;
  int idf = id0+l-1;
  int lconc = n*nbcol - (nocol[0]==0)*(id0%n) - (nocol[nbcol]==m)*(n-idf%n);

  fprintf(file, "lconc=%d\n", lconc);

  int distcorrmin= lambda-1;
  int v1_size=lconc+(distcorrmin)*(nbcol-1);
  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;

  int nrshp=fft_size;
  int mrshp=nfft;
  int lrshp=nrshp*mrshp;

  fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);
  fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);

  double *Vrshp;
  Vrshp = (double *) calloc(lrshp, sizeof(double));

  for(i=0;i<lrshp;i++)
    Vrshp[i]=i+2;

  int flag_format_rshp=2;

  extract_result(Vout, nocol, nbcol, lconc, n, m, id0, l, lambda, nfft, Vrshp, nrshp, mrshp, lrshp, flag_format_rshp);


  fprintf(file, "---------------\n");
  for(i=0;i<l;i++) {
  if ((i+id0)%n==0)
    fprintf(file, "--\n");
  fprintf(file, "Vout[%d]=%f\n", i, Vout[i]);
  }

}//end test3


//---------------------------------------------------
 if (test==4) {
  fprintf(file, "fctid_vect2nfftblock-->vect2mat\n");
  fprintf(file, "---------------\n");

  id0=1;
  n=3;
  m=3;
  l=8;
  lambda=3;

//  int distcorrmin=lambda-1;
//  int m_eff=m;
//  int v1_size = l + (m_eff-1)*distcorrmin;
  int nfft=2;
//  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;
//  int rfirst=id0%n;

  int nbcol=3;
  int *nocol;
  nocol = (int *) calloc(nbcol, sizeof(double));
  nocol[0]=0;
  nocol[1]=1;
  nocol[2]=2;

//  nrshp=fft_size;
//  mrshp=nfft;
//  lrshp=nrshp*mrshp;

 // fprintf(file, "nrshp=%d , mrshp=%d\n", nrshp, mrshp);
 // fprintf(file, "fft_size=%d , v1_size=%d\n", fft_size, v1_size);
  fprintf(file, "n=%d , m=%d\n", n, m);

  double *Vin;
  Vin = (double *) calloc(l, sizeof(double));

  for(i=0;i<l;i++)
    Vin[i]=i+1;

//  int lconc=6;
  int idf = id0+l-1;
  int lconc = n*nbcol - (nocol[0]==0)*(id0%n) - (nocol[nbcol]==m)*(n-idf%n);

  fprintf(file, "lconc=%d\n", lconc);

  int distcorrmin= lambda-1;
  int v1_size=lconc+(distcorrmin)*(nbcol-1);
  int fft_size = ceil(1.0*v1_size/nfft)+2*distcorrmin;

  int nrshp=fft_size;
  int mrshp=nfft;
  int lrshp=nrshp*mrshp;

  double *Vrshp;
  Vrshp = (double *) calloc(lrshp, sizeof(double));

  double *Vout;
  Vout = (double *) calloc(l, sizeof(double));


  int flag_outformat=2;

  build_reshape(Vin, nocol, nbcol, lconc, n, m, id0, l, lambda, nfft, Vrshp, nrshp, mrshp, lrshp, flag_outformat);


  fprintf(file, "---------------\n");
  for(i=0;i<lrshp;i++) {
  if (i%nrshp==0)
    fprintf(file, "--\n");
  fprintf(file, "Vrshp[%d]=%f\n", i, Vrshp[i]);
  }


  int flag_format_rshp=2;
  extract_result(Vout, nocol, nbcol, lconc, n, m, id0, l, lambda, nfft, Vrshp, nrshp, mrshp, lrshp, flag_format_rshp);


  fprintf(file, "---------------\n");
  for(i=0;i<l;i++) {
  if ((i+id0)%n==0)
    fprintf(file, "--\n");
  fprintf(file, "Vout[%d]=%f\n", i, Vout[i]);
  }


}//end test4





return 0;

}



