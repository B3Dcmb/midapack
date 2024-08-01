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
#include "templates.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// choose header based on compilation option
#ifdef HAVE_MKL
#include <mkl.h>
#else
#include <lapacke.h>
#endif
/******************************************************************************/
/*                          Defining Algebra routines                         */
/******************************************************************************/
/* Projecting templates amplitudes in time domain */
int TVecProd(TemplateClass *X, int nces, int m, double sampling_freq, int *scan_size,
  int **sweeptstamps, int **az_binned, double ***hwp_mod, int nhwp, double delta_t,
  int store_hwp, double *tau, double *out){

  int i,j,k,ces_id;
  hwpss_w hwpss_wghts;
  //initialization of output vector
  for(i=X->tinit;i<=(X+nces*m-1)->tlast;i++)
    out[i] = 0;
  //operation
  for(ces_id=0;ces_id<nces;ces_id++){
    if(store_hwp == 0)
      build_hwpss_w(&hwpss_wghts, hwp_mod[ces_id], scan_size[ces_id], nhwp, ces_id);
    for(k=0;k<m;k++){
      if(strcmp((X+(ces_id*m+k))->flag_w,"OTF") == 0){
        (X+(ces_id*m+k))->bins = (int *) calloc((X+(ces_id*m+k))->nsamples * (X+(ces_id*m+k))->nmult, sizeof(int));
        if(strcmp((X+(ces_id*m+k))->ID,"POLY") == 0){
          (X+(ces_id*m+k))->wghts = (double *) calloc((X+(ces_id*m+k))->nsamples * (X+(ces_id*m+k))->nmult, sizeof(double));
          expandpolydata(X+(ces_id*m+k), (X+(ces_id*m+k))->bins, (X+(ces_id*m+k))->wghts, sweeptstamps[ces_id], sampling_freq, (X+(ces_id*m+k))->order);
        }
        else if(strcmp((X+(ces_id*m+k))->ID,"SSS") == 0){
          expandSSSdata(X+(ces_id*m+k), (X+(ces_id*m+k))->bins, az_binned[ces_id]);
        }
        else if(strcmp((X+(ces_id*m+k))->ID,"HWPSS") == 0){
          (X+(ces_id*m+k+1))->bins = (int *) calloc((X+(ces_id*m+k+1))->nsamples * (X+(ces_id*m+k+1))->nmult, sizeof(int));
          expandHWPSSdata(X+(ces_id*m+k), (X+(ces_id*m+k))->bins, &(X+(ces_id*m+k))->wghts,
            (X+(ces_id*m+k+1))->bins, &(X+(ces_id*m+k+1))->wghts, delta_t, sampling_freq,
            (X+(ces_id*m+k))->order, hwpss_wghts);
        }
      }
      //Case1: Polynomial template
      if(strcmp((X+(ces_id*m+k))->ID,"POLY") == 0){
        for(i=0;i<(X+(ces_id*m+k))->nsamples;i++){
          for(j=0;j<(X+(ces_id*m+k))->nmult;j++){
            out[(X+(ces_id*m+k))->tinit + i] += (X+(ces_id*m+k))->wghts[i * (X+(ces_id*m+k))->nmult + j] * tau[(X+(ces_id*m+k))->bins[i * (X+(ces_id*m+k))->nmult + j]];
          }
        }
      }
      //Case2: SSS template
      else if(strcmp((X+(ces_id*m+k))->ID,"SSS") == 0){
        for(i=0;i<(X+(ces_id*m+k))->nsamples;i++){
          for(j=0;j<(X+(ces_id*m+k))->nmult;j++){
            out[(X+(ces_id*m+k))->tinit + i] +=  tau[(X+(ces_id*m+k))->bins[i * (X+(ces_id*m+k))->nmult + j]];
          }
        }
      }
      //Case3: HWPSS template
      else if(strcmp((X+(ces_id*m+k))->ID,"HWPSS") == 0){
        for(i=0;i<(X+(ces_id*m+k))->nsamples;i++){
          for(j=0;j<(X+(ces_id*m+k))->nmult;j++){
            out[(X+(ces_id*m+k))->tinit + i] += (X+(ces_id*m+k))->wghts[i * (X+(ces_id*m+k))->nmult + j] * tau[(X+(ces_id*m+k))->bins[i * (X+(ces_id*m+k))->nmult + j]]; // cosine terms
            out[(X+(ces_id*m+k+1))->tinit + i] += (X+(ces_id*m+k+1))->wghts[i * (X+(ces_id*m+k+1))->nmult + j] * tau[(X+(ces_id*m+k+1))->bins[i * (X+(ces_id*m+k+1))->nmult + j]]; // sine terms
          }
        }
      }
      // Else, something wrong is going on ...
      else{
        printf("Run time error while projecting templates in time domain: invalid template type, please verify templates ID\n");
        exit(1);
      }

      if(strcmp((X+(ces_id*m+k))->flag_w,"OTF") == 0){
        free((X+(ces_id*m+k))->bins);
        if(strcmp((X+(ces_id*m+k))->ID,"POLY") == 0)
          free((X+(ces_id*m+k))->wghts);
        else if(strcmp((X+(ces_id*m+k))->ID,"HWPSS") == 0){
          free((X+(ces_id*m+k+1))->bins);
          k++;
        }
      }
      else{
        if(strcmp((X+(ces_id*m+k))->ID,"HWPSS") == 0)
          k++;
      }
    }
    if(store_hwp == 0)
      free_hwpss_w(&hwpss_wghts, nhwp);
  }
  return 0;
}

/* Projecting time domain in templates space */
int TrTVecProd(TemplateClass *X, int nces, int m, double sampling_freq, int *scan_size,
  int **sweeptstamps, int **az_binned, double ***hwp_mod, int nhwp, double delta_t,
  int store_hwp, double *d, double *out){

  int i,j,k,ces_id;
  hwpss_w hwpss_wghts;
  // initialize output vector
  for(i=X->nbinMin;i<=(X+nces*m-1)->nbinMax;i++)
    out[i] = 0;
  // operation
  for(ces_id=0;ces_id<nces;ces_id++){
    if(store_hwp == 0)
      build_hwpss_w(&hwpss_wghts, hwp_mod[ces_id], scan_size[ces_id], nhwp, ces_id);
    for(k=0;k<m;k++){
      if(strcmp((X+(ces_id*m+k))->flag_w,"OTF") == 0){
        (X+(ces_id*m+k))->bins = (int *) calloc((X+(ces_id*m+k))->nsamples * (X+(ces_id*m+k))->nmult, sizeof(int));
        if(strcmp((X+(ces_id*m+k))->ID,"POLY") == 0){
          (X+(ces_id*m+k))->wghts = (double *) calloc((X+(ces_id*m+k))->nsamples * (X+(ces_id*m+k))->nmult, sizeof(double));
          expandpolydata(X+(ces_id*m+k), (X+(ces_id*m+k))->bins, (X+(ces_id*m+k))->wghts, sweeptstamps[ces_id], sampling_freq, (X+(ces_id*m+k))->order);
        }
        else if(strcmp((X+(ces_id*m+k))->ID,"SSS") == 0){
          expandSSSdata(X+(ces_id*m+k), (X+(ces_id*m+k))->bins, az_binned[ces_id]);
        }
        else if(strcmp((X+(ces_id*m+k))->ID,"HWPSS") == 0){
          (X+(ces_id*m+k+1))->bins = (int *) calloc((X+(ces_id*m+k+1))->nsamples * (X+(ces_id*m+k+1))->nmult, sizeof(int));
          expandHWPSSdata(X+(ces_id*m+k), (X+(ces_id*m+k))->bins, &(X+(ces_id*m+k))->wghts,
            (X+(ces_id*m+k+1))->bins, &(X+(ces_id*m+k+1))->wghts, delta_t, sampling_freq,
            (X+(ces_id*m+k))->order, hwpss_wghts);
        }
      }
      //Case1: Polynomial template
      if(strcmp((X+(ces_id*m+k))->ID,"POLY") == 0){
        for(i=0;i<(X+(ces_id*m+k))->nsamples;i++){
          for(j=0;j<(X+(ces_id*m+k))->nmult;j++){
            out[(X+(ces_id*m+k))->bins[i * (X+(ces_id*m+k))->nmult + j]] += (X+(ces_id*m+k))->wghts[i * (X+(ces_id*m+k))->nmult + j] * d[(X+(ces_id*m+k))->tinit + i];
          }
        }
      }
      //Case2: SSS template
      else if(strcmp((X+(ces_id*m+k))->ID,"SSS") == 0){
        for(i=0;i<(X+(ces_id*m+k))->nsamples;i++){
          for(j=0;j<(X+(ces_id*m+k))->nmult;j++){
            out[(X+(ces_id*m+k))->bins[i * (X+(ces_id*m+k))->nmult + j]] +=  d[(X+(ces_id*m+k))->tinit + i];
          }
        }
      }
      //Case3: HWPSS template
      else if(strcmp((X+(ces_id*m+k))->ID,"HWPSS") == 0){
        for(i=0;i<(X+(ces_id*m+k))->nsamples;i++){
          for(j=0;j<(X+(ces_id*m+k))->nmult;j++){
            out[(X+(ces_id*m+k))->bins[i * (X+(ces_id*m+k))->nmult + j]] += (X+(ces_id*m+k))->wghts[i * (X+(ces_id*m+k))->nmult + j] * d[(X+(ces_id*m+k))->tinit + i]; // cosine terms
            out[(X+(ces_id*m+k+1))->bins[i * (X+(ces_id*m+k+1))->nmult + j]] += (X+(ces_id*m+k+1))->wghts[i * (X+(ces_id*m+k+1))->nmult + j] * d[(X+(ces_id*m+k+1))->tinit + i]; // sine terms
          }
        }
      }
      // Else, something wrong is going on ...
      else{
        printf("Run time error while projecting templates in time domain: invalid template type, please verify templates ID\n");
        exit(1);
      }

      if(strcmp((X+(ces_id*m+k))->flag_w,"OTF") == 0){
        free((X+(ces_id*m+k))->bins);
        if(strcmp((X+(ces_id*m+k))->ID,"POLY") == 0)
          free((X+(ces_id*m+k))->wghts);
        else if(strcmp((X+(ces_id*m+k))->ID,"HWPSS") == 0){
          free((X+(ces_id*m+k+1))->bins);
          k++;
        }
      }
      else{
        if(strcmp((X+(ces_id*m+k))->ID,"HWPSS") == 0)
          k++;
      }
    }
    if(store_hwp == 0)
      free_hwpss_w(&hwpss_wghts, nhwp);
  }
  return 0;
}

/* Building Kernel Blocks */
int BuildKernel(TemplateClass *X, int n, double *B, double w, int *sweeptstamps,
  int *az_binned, hwpss_w hwpss_wghts, double delta_t, double sampling_freq){
  // n : number of template classes in one kernel Block (1 det data in 1 CES)
  int i,j,l,k;
  int m1,m2;
  double sum=0;
  for(k=0;k<n;k++){
    if(strcmp((X+k)->flag_w,"OTF") == 0){
      // printf("I'm here\n");
      // fflush(stdout);
      (X+k)->bins = (int *) calloc((X+k)->nsamples * (X+k)->nmult, sizeof(int));
      if(strcmp((X+k)->ID,"POLY") == 0){
        (X+k)->wghts = (double *) calloc((X+k)->nsamples * (X+k)->nmult, sizeof(double));
        expandpolydata(X+k, (X+k)->bins, (X+k)->wghts, sweeptstamps, sampling_freq, (X+k)->order);
      }
      else if(strcmp((X+k)->ID,"SSS") == 0)
        expandSSSdata(X+k, (X+k)->bins, az_binned);
      else if(strcmp((X+k)->ID,"HWPSS") == 0){
        (X+k+1)->bins = (int *) calloc((X+k+1)->nsamples * (X+k+1)->nmult, sizeof(int));
        expandHWPSSdata(X+k, (X+k)->bins, &(X+k)->wghts, (X+k+1)->bins, &(X+k+1)->wghts,
          delta_t, sampling_freq, (X+k)->order, hwpss_wghts);
      }
    }
    for(l=0;l<n;l++){
      if((strcmp((X+l)->flag_w,"OTF") == 0) && (l!=k)){
        (X+l)->bins = (int *) calloc((X+k)->nsamples * (X+l)->nmult, sizeof(int));
        if(strcmp((X+l)->ID,"POLY") == 0){
          (X+l)->wghts = (double *) calloc((X+k)->nsamples * (X+l)->nmult, sizeof(double));
          expandpolydata(X+l, (X+l)->bins, (X+l)->wghts, sweeptstamps, sampling_freq, (X+l)->order);
        }
        else if(strcmp((X+l)->ID,"SSS") == 0)
          expandSSSdata(X+l, (X+l)->bins, az_binned);
        else if(strcmp((X+l)->ID,"HWPSS") == 0){
          (X+l+1)->bins = (int *) calloc((X+l+1)->nsamples * (X+l+1)->nmult, sizeof(int));
          expandHWPSSdata(X+l, (X+l)->bins, &(X+l)->wghts, (X+l+1)->bins, &(X+l+1)->wghts,
            delta_t, sampling_freq, (X+l)->order, hwpss_wghts);
        }
      }


      // if((X+l)->order == 1){
      //   for(i=0;i<10;i++){
      //     printf("bins[%d] = %d\n",i,(X+l)->bins[i]);
      //     printf("wghts[%d] = %f\n",i,(X+l)->wghts[i]);
      //   }
      // }
      // if((X+l)->order == 1){
      //   for(i=0;i<sweeptstamps[1];i++){
      //     printf("wghts[%d] = %6.2f\n",i,(X+l)->wghts[i]);
      //     sum += (X+l)->wghts[i];
      //   }
      //   printf("sum wgths = %6.2f\n",sum);
      // }

      // Case1: poly/HWPSS x poly/HWPSS kernel
      if(((strcmp((X+k)->ID,"POLY") == 0) || (strcmp((X+k)->ID,"HWPSS") == 0)) && ((strcmp((X+l)->ID,"POLY") == 0) || (strcmp((X+l)->ID,"HWPSS") == 0))){
        if((strcmp((X+k)->ID,"HWPSS") == 0) && (strcmp((X+l)->ID,"HWPSS") == 0)){
          j=0;
          for(i=0;i<(X+k)->nsamples;i++){
            for(m1=0;m1<(X+k)->nmult;m1++){
              for(m2=0;m2<(X+l)->nmult;m2++){
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k)->wghts[i*(X+k)->nmult + m1] * (X+l)->wghts[i*(X+l)->nmult + m2];
                B[((X+k+1)->bins[i*(X+k+1)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k+1)->wghts[i*(X+k+1)->nmult + m1] * (X+l)->wghts[i*(X+l)->nmult + m2];
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l+1)->bins[i*(X+l+1)->nmult + m2]-X->nbinMin)] += w * (X+k)->wghts[i*(X+k)->nmult + m1] * (X+l+1)->wghts[i*(X+l+1)->nmult + m2];
                B[((X+k+1)->bins[i*(X+k+1)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l+1)->bins[i*(X+l+1)->nmult + m2]-X->nbinMin)] += w * (X+k+1)->wghts[i*(X+k+1)->nmult + m1] * (X+l+1)->wghts[i*(X+l+1)->nmult + m2];
              }
            }
          }
        }
        else if((strcmp((X+k)->ID,"HWPSS") == 0) && (strcmp((X+l)->ID,"HWPSS") != 0)){
          for(i=0;i<(X+k)->nsamples;i++){
            for(m1=0;m1<(X+k)->nmult;m1++){
              for(m2=0;m2<(X+l)->nmult;m2++){
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k)->wghts[i*(X+k)->nmult + m1] * (X+l)->wghts[i*(X+l)->nmult + m2];
                B[((X+k+1)->bins[i*(X+k+1)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k+1)->wghts[i*(X+k+1)->nmult + m1] * (X+l)->wghts[i*(X+l)->nmult + m2];
              }
            }
          }
        }
        else if((strcmp((X+k)->ID,"HWPSS") != 0) && (strcmp((X+l)->ID,"HWPSS") == 0)){
          for(i=0;i<(X+k)->nsamples;i++){
            for(m1=0;m1<(X+k)->nmult;m1++){
              for(m2=0;m2<(X+l)->nmult;m2++){
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k)->wghts[i*(X+k)->nmult + m1] * (X+l)->wghts[i*(X+l)->nmult + m2];
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l+1)->bins[i*(X+l+1)->nmult + m2]-X->nbinMin)] += w * (X+k)->wghts[i*(X+k)->nmult + m1] * (X+l+1)->wghts[i*(X+l+1)->nmult + m2];
              }
            }
          }
        }
        else{
          for(i=0;i<(X+k)->nsamples;i++){
            for(m1=0;m1<(X+k)->nmult;m1++){
              for(m2=0;m2<(X+l)->nmult;m2++){
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k)->wghts[i*(X+k)->nmult + m1] * (X+l)->wghts[i*(X+l)->nmult + m2];
              }
            }
          }
        }
      }
      // Case2: SSS x poly/HWPSS kernel
      else if((strcmp((X+k)->ID,"SSS") == 0) && ((strcmp((X+l)->ID,"POLY") == 0) || (strcmp((X+l)->ID,"HWPSS") == 0))){
        if(strcmp((X+l)->ID,"HWPSS") == 0){
          for(i=0;i<(X+k)->nsamples;i++){
            for(m1=0;m1<(X+k)->nmult;m1++){
              for(m2=0;m2<(X+l)->nmult;m2++){
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+l)->wghts[i*(X+l)->nmult + m2];
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l+1)->bins[i*(X+l+1)->nmult + m2]-X->nbinMin)] += w * (X+l+1)->wghts[i*(X+l+1)->nmult + m2];
              }
            }
          }
        }
        else{
          for(i=0;i<(X+k)->nsamples;i++){
            for(m1=0;m1<(X+k)->nmult;m1++){
              for(m2=0;m2<(X+l)->nmult;m2++){
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+l)->wghts[i*(X+l)->nmult + m2];
              }
            }
          }
        }
      }
      // Case3: poly/HWPSS x SSS kernel
      else if(((strcmp((X+k)->ID,"POLY") == 0) || (strcmp((X+k)->ID,"HWPSS") == 0)) && (strcmp((X+l)->ID,"SSS") == 0)){
        if(strcmp((X+k)->ID,"HWPSS") == 0){
          for(i=0;i<(X+k)->nsamples;i++){
            for(m1=0;m1<(X+k)->nmult;m1++){
              for(m2=0;m2<(X+l)->nmult;m2++){
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k)->wghts[i*(X+k)->nmult + m1];
                B[((X+k+1)->bins[i*(X+k+1)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k+1)->wghts[i*(X+k+1)->nmult + m1];
              }
            }
          }
        }
        else{
          for(i=0;i<(X+k)->nsamples;i++){
            for(m1=0;m1<(X+k)->nmult;m1++){
              for(m2=0;m2<(X+l)->nmult;m2++){
                B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w * (X+k)->wghts[i*(X+k)->nmult + m1];
              }
            }
          }
        }
      }
      // Case4: SSS x SSS kernel
      else if((strcmp((X+k)->ID,"SSS") == 0) && (strcmp((X+l)->ID,"SSS") == 0)){
        for(i=0;i<(X+k)->nsamples;i++){
          for(m1=0;m1<(X+k)->nmult;m1++){
            for(m2=0;m2<(X+l)->nmult;m2++){
              B[((X+k)->bins[i*(X+k)->nmult + m1]-X->nbinMin) * ((X+n-1)->nbinMax-X->nbinMin+1) + ((X+l)->bins[i*(X+l)->nmult + m2]-X->nbinMin)] += w ;
            }
          }
        }
      }
      // Else, something wrong is going on ...
      else{
        printf("Run time error while computing templates kernel: invalid templates combination, please verify templates ID\n");
        exit(1);
      }
      // for(i=0;i<(X+k)->nsamples * (X+k)->nmult;i++){
      //   for(j=0;j<(X+l)->nsamples * (X+l)->nmult;j++){
      //     // VERY TIME CONSUMING ON SIMPLE TEST JUMPED FROM 0.02s to ~5s !
      //     // if((i<5) && (j<5))
      //     // printf("\nindex of B = %d\n",((X+k)->bins[i]-(X+k)->nbinMin) * (X+l)->nbins + ((X+l)->bins[j]-(X+l)->nbinMin));
      //     B[((X+k)->bins[i]-(X+k)->nbinMin) * (X+l)->nbins + ((X+l)->bins[j]-(X+l)->nbinMin)] += w * (X+k)->wghts[i] * (X+l)->wghts[j];
      //   }
      // }
      if((strcmp((X+l)->flag_w, "OTF") == 0) && (l!=k)){
        free((X+l)->bins);
        if(strcmp((X+l)->ID,"POLY") == 0)
          free((X+l)->wghts);
        else if(strcmp((X+l)->ID,"HWPSS") == 0){
          free((X+l+1)->bins);
          l++;
        }
      }
      else{
        if(strcmp((X+l)->ID,"HWPSS") == 0)
          l++;
      }
    }
    if(strcmp((X+k)->flag_w, "OTF") == 0){
      free((X+k)->bins);
      if(strcmp((X+k)->ID,"POLY") == 0)
        free((X+k)->wghts);
      else if(strcmp((X+k)->ID,"HWPSS") == 0){
        free((X+k+1)->bins);
        k++;
      }
    }
    else{
      if(strcmp((X+k)->ID,"HWPSS") == 0)
        k++;
    }
  }
  return 0;
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
  double epsilon = 1e-12; // threshold for regularization
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
  free(s);
  return rank;
}

void transpose_nn(double *A, int n)
{
    int i, j;
    double temp;

    for(i = 0; i < n-1 ; i++)
        for (j = i+1; j < n; j++)
        {
            temp = A[i*n+j];
            A[i*n+j] = A[j*n+i];
            A[j*n+i] = temp;
        }

}
/* Test inverse routine */
int inverse_svd(int m, int n, int lda,  double *a)
{
    // A = UDVt
    // A = USVt
    int i, j, k;
    int rank;

    // Setup a buffer to hold the singular values:
    int nsv = m < n ? m : n;
    double *s = malloc(nsv * sizeof(double)); // = D

    // Setup buffers to hold the matrices U and Vt:
    double *u = malloc(m*m * sizeof(double)); // = U
    double *vt = malloc(n*n * sizeof(double));

    // Workspace and status variables:
    /*
    double workSize;//
    double *work = &workSize;//
    int lwork = -1;//
    int *iwork = malloc(8 * nsv * sizeof(int));//
    */
    //double *superb = malloc((nsv-1) * sizeof(double));



    int info = 0;

    transpose_nn(a, m);

    info = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'A', m, n, a, lda, s, u, m, vt, n);
    //info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'A', 'A', m, n, a, lda, s, u, m, vt, n, superb);
    /*
    // Call dgesdd_ with lwork = -1 to query optimal workspace size:
    dgesvd_("A", "A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, &info);
    //dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);
 // info = LAPACKE_dgesdd_work(LAPACK_COL_MAJOR, 'A', m, n, a, lda, s, u, m, vt, n, work, lwork, iwork);
     //info = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'A', m, n, a, lda, s, u, m, vt, n);
    if (info > 0) {
        printf( "The algorithm computing SVD failed to converge (1).\n" );
        exit(1);
    }
    // Optimal workspace size is returned in work[0].
    lwork = workSize;
    work = malloc(lwork * sizeof(double));
    // Call dgesdd_ to do the actual computation:
    dgesvd_("A", "A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, &info);
    //dgesdd_("A", &m, &n, a, &lda, s, u, &m, vt, &n, work, &lwork, iwork, &info);

    */

    if (info > 0) {
        printf( "The algorithm computing SVD failed to converge (1).\n" );
        exit(1);
    }
    if (info < 0)
    { printf( "General error .\n" );
        exit(1);
    }

    // Cleanup workspace:
    // free(work);
    // free(iwork);

    // Computing S-1
    for (k = 0; k < nsv; ++k) {
        if (fabs(s[k]) < 1.0e-15) s[k] = 0.0;
        else {
          s[k] = 1.0 / s[k];
          rank++;
        }
    }

    // do something useful with U, S, Vt ...
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            a[i * m + j] = 0.0;
            for (k = 0; k < n; k++) {
                a[i * m + j] += vt[i * m + k] * s[k] * u[k * m + j];
            }
        }
    }

    // and then clean them up too:
    free(s);
    free(u);
    free(vt);
return rank;
}

/* Define Legendre Polynomials */
double P0(double x){
    return 1;
}

double P1(double x){
    return x;
}

//The following is a general function that returns the value of the Legendre Polynomial for any given x and n=0,1,2,3,...
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

/* Build Legendre polynomial weights through recursive relations */
// void build_poly_w(poly_w *poly_wghts, int *sweeptstamps,int size, int order, int ces_id){
//
//   int i,j;
//   //Initialization of offset
//   poly_wgths->poly = (double **) malloc((order+1) * sizeof(double*));
//   poly_wghts->ces_id = ces_id;
//   poly_wgths->poly[0] = (double *) calloc(size, sizeof(double));
//   for(j=0;j<size;j++){
//     poly_wgths->poly[0][j] = 1;
//   }
//   if(order == 0)
//     return;
//
//   //Linear order
//   for(i=0;i<X->nbins;i++){
//     for(j=sweeptstamps[i];j<sweeptstamps[i+1];j++){
//       bins[j] = X->nbinMin + i;
//       // printf("expandpoly: bins[%d]=%d\n",j,bins[j]);
//       //N.B: the relation is rescaled so sampling_freq shouldn't really matter
//       wghts[j] = Legendre((double)j/sampling_freq, (double)sweeptstamps[i]/sampling_freq, (double)(sweeptstamps[i+1]-1)/sampling_freq, order);
//       // printf("expandpoly: wghts[%d]=%f\n",j,wghts[j]);
//
//     }
//   }
//
// }

/* Build HWPSS harmonics weights through recursive computation of trig functions */
void build_hwpss_w(hwpss_w *hwpss_wghts, double **hwp_mod, int size, int order,
  int ces_id){

  int i,j;
  if(order == 0)
    return;
  //Initialization of fundamental harmonic & ces_id
  hwpss_wghts->hwpcos = (double **) malloc(order * sizeof(double*));
  hwpss_wghts->hwpsin = (double **) malloc(order * sizeof(double*));
  hwpss_wghts->ces_id = ces_id;
  hwpss_wghts->hwpcos[0] = (double *) calloc(size, sizeof(double));
  hwpss_wghts->hwpsin[0] = (double *) calloc(size, sizeof(double));
  memcpy(hwpss_wghts->hwpcos[0],hwp_mod[0],size*sizeof(double));
  memcpy(hwpss_wghts->hwpsin[0],hwp_mod[1],size*sizeof(double));
  if(order == 1)
    return;

  // Recursive computation
  for(i=2;i<=order;i++){
    hwpss_wghts->hwpcos[i-1] = (double *) calloc(size, sizeof(double));
    hwpss_wghts->hwpsin[i-1] = (double *) calloc(size, sizeof(double));
    for(j=0;j<size;j++){// recursive trig relation to compute the next harmonic
      hwpss_wghts->hwpcos[i-1][j] = hwpss_wghts->hwpcos[i-2][j] * hwpss_wghts->hwpcos[0][j] - hwpss_wghts->hwpsin[i-2][j] * hwpss_wghts->hwpsin[0][j];
      hwpss_wghts->hwpsin[i-1][j] = hwpss_wghts->hwpsin[i-2][j] * hwpss_wghts->hwpcos[0][j] + hwpss_wghts->hwpcos[i-2][j] * hwpss_wghts->hwpsin[0][j];
    }
  }
  return;
}

/* Free the CES HWPSS data */
void free_hwpss_w(hwpss_w *hwpss_wghts, int order){

  int i;
  for(i=0;i<order;i++){
    free(hwpss_wghts->hwpcos[i]);
    free(hwpss_wghts->hwpsin[i]);
  }
  if(order != 0){
    free(hwpss_wghts->hwpcos);
    free(hwpss_wghts->hwpsin);
  }
}

/******************************************************************************/
/*        Utility routines for building templates classes objects             */
/******************************************************************************/
/* Build local templates classes list */
/* Parameters:
  - Input:
    @X: The templates classes list
    @ndet: local number of detectors assigned to the process
    @npoly: maximum order considered for polynomial templates
    N.B: Will introduce nground, nhwp and ncustom when we implement support for
    other templates classes types. Flags are ignored for now.
*/
int Tlist_init(TemplateClass *X, int ndet, int nces, int *block_nsamples, int **detnsweeps,
  int *scan_size, int **sweeptstamps, int n_sss_bins, int **az_binned, double sampling_freq,
  int npoly, int ground, int nhwp, double delta_t, int store_hwp, double ***hwp_mod){

  int i,j,k;
  int tinit,tlast,nbinMin, nbinMax;
  int n_class = npoly + ground + 2*nhwp;
  hwpss_w hwpss_wghts;
  // looping over local list of detectors
  tinit = 0;
  tlast = -1;
  nbinMin = 0;
  nbinMax = -1;
  for(i=0;i<nces;i++){
    if(store_hwp != 0)
      build_hwpss_w(&hwpss_wghts, hwp_mod[i], scan_size[i], nhwp, i);
    for(j=0;j<ndet;j++){
      tlast += block_nsamples[i*ndet+j];
      // looping over polynomial templates for each detector
      for(k=0;k<npoly;k++){
        // looping over polynomial orders
        nbinMax += detnsweeps[i][j];
        // For each order build the corresponding polynomial template class
        Polyinit(X + ((i*ndet+j)*n_class + k), tinit, tlast, nbinMin, nbinMax, sweeptstamps[i], sampling_freq, k);
        nbinMin += detnsweeps[i][j];
      }
      if(ground == 1){
        nbinMax += n_sss_bins;
        SSSinit((X + ((i*ndet+j)*n_class) + npoly), tinit, tlast, nbinMin, nbinMax, az_binned[i]);
        nbinMin += n_sss_bins;
      }
      // looping over hwp harmonics orders
      for(k=0;k<nhwp;k++){
        nbinMax += (int)ceil((tlast-tinit+1)/(delta_t*sampling_freq));
        HWPSSinit((X + ((i*ndet+j)*n_class) + npoly+ ground + 2*k), tinit, tlast, nbinMin, nbinMax, delta_t, sampling_freq, k+1, hwpss_wghts);
        nbinMax += (int)ceil((tlast-tinit+1)/(delta_t*sampling_freq));
        nbinMin += 2*((int)ceil((tlast-tinit+1)/(delta_t*sampling_freq)));
      }
      tinit += block_nsamples[i*ndet+j];
    }
  }
  return 0;
}

/* Build polynomial template class instance */
int Polyinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  int *sweeptstamps, double sampling_freq, int order){

  char *flag = (char *) malloc(10 * sizeof(char));
  char *flag_w = (char *) malloc(10 * sizeof(char));
  int *bins;
  double *wghts;
  sprintf(flag,"POLY");
  sprintf(flag_w,"OTF");
  TCinit(X, tinit, tlast, nbinMin, nbinMax, 1, flag, flag, flag ,flag_w, flag);
  X->order = order;
  if(strcmp(X->flag_w,"S")==0){ // Stored weights
    X->bins = (int *) calloc(X->nsamples * X->nmult,sizeof(int));
    X->wghts = (double *) calloc(X->nsamples * X->nmult, sizeof(double));
    expandpolydata(X, X->bins, X->wghts, sweeptstamps, sampling_freq, order);
  }
  return 0;
}

/* Build SSS template class instance */
int SSSinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  int *az_binned){

  char *flag = (char *) malloc(10 * sizeof(char));
  char *flag_w = (char *) malloc(10 * sizeof(char));
  int *bins;
  sprintf(flag,"SSS");
  sprintf(flag_w,"OTF");
  TCinit(X, tinit, tlast, nbinMin, nbinMax, 1, flag, flag, flag ,flag_w, flag);
  if(strcmp(X->flag_w,"S")==0){ // Stored weights
    X->bins = (int *) calloc(X->nsamples * X->nmult,sizeof(int));
    expandSSSdata(X, X->bins, az_binned);
  }
  return 0;
}

/* Build HWP template class instance */
int HWPSSinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  double delta_t, double sampling_freq, int order, hwpss_w hwpss_wghts){

  char *flag = (char *) malloc(10 * sizeof(char));
  char *flag_w = (char *) malloc(10 * sizeof(char));
  sprintf(flag,"HWPSS");
  sprintf(flag_w,"OTF");
  TCinit(X, tinit, tlast, nbinMin, nbinMax, 1, flag, flag, flag ,flag_w, flag);
  TCinit(X+1, tinit, tlast, nbinMax+1, nbinMax+(nbinMax-nbinMin+1), 1, flag, flag, flag ,flag_w, flag);
  X->order = order;
  (X+1)->order = order;
  if(strcmp(X->flag_w,"S")==0 && strcmp((X+1)->flag_w,"S")==0){ // Stored weights
    X->bins = (int *) calloc(X->nsamples * X->nmult,sizeof(int));
    (X+1)->bins = (int *) calloc((X+1)->nsamples * (X+1)->nmult,sizeof(int));
    expandHWPSSdata(X, X->bins, &X->wghts, (X+1)->bins, &(X+1)->wghts, delta_t, sampling_freq, order, hwpss_wghts);
  }
  return 0;
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

  return 0;
}

/* Build bins and weights profile for a polynomial template nmult = 1 for now */
int expandpolydata(TemplateClass *X ,int *bins, double *wghts,
  int *sweeptstamps, double sampling_freq, int order){

  int i, j;
  for(i=0;i<X->nbins;i++){
    for(j=sweeptstamps[i];j<sweeptstamps[i+1];j++){
      bins[j] = X->nbinMin + i;
      // printf("expandpoly: bins[%d]=%d\n",j,bins[j]);
      //N.B: the relation is rescaled so sampling_freq shouldn't really matter
      wghts[j] = Legendre((double)j/sampling_freq, (double)sweeptstamps[i]/sampling_freq, (double)(sweeptstamps[i+1]-1)/sampling_freq, order);
      // printf("expandpoly: wghts[%d]=%f\n",j,wghts[j]);

    }
  }
  return 0;
}

/* Build bins profile for a SSS template nmult = 1 for now */
int expandSSSdata(TemplateClass *X, int *bins, int *az_binned){
  int i;
  for(i=0;i<X->nsamples * X->nmult;i++){
    bins[i] = X->nbinMin + az_binned[i];
  }
  return 0;
}

/* Build bins and weights profile for a HWPSS template nmult = 1 for now */
int expandHWPSSdata(TemplateClass *X, int *bins0, double **wghts0, int *bins1,
  double **wghts1, double delta_t, double sampling_freq, int order,
  hwpss_w hwpss_wghts){
  int i;
  for(i=0;i<X->nsamples * X->nmult;i++){
    bins0[i] = X->nbinMin + (int)floor((double)i/(sampling_freq*delta_t));
    bins1[i] = (X+1)->nbinMin + (int)floor((double)i/(sampling_freq*delta_t));
  }
  // for(i=8300;i<8400;i++){
  //   if(X->order == 1){
  //     printf("X->bins[%d]=%d\n",i,bins0[i]);
  //     printf("X+1->bins[%d]=%d\n",i,bins1[i]);
  //   }
  // }
  *wghts0 = hwpss_wghts.hwpcos[order-1];
  *wghts1 = hwpss_wghts.hwpsin[order-1];
  return 0;
}

/* Bin boresight azimuth array */
int** bin_az(double **az, double *az_min, double *az_max, int *ces_length, int sss, int n_sss_bins, int nces)
{
  int i,j;
  int **az_binned = (int **) malloc(nces * sizeof(int*));

  // Build binned boresight azimuth array
  if(sss){
    for(i=0; i<nces;i++){
      az_binned[i] = (int *) calloc(ces_length[i], sizeof(int));
      for(j=0;j<ces_length[i];j++){
        az_binned[i][j] = miin(n_sss_bins-1,(int)floor((az[i][j]-az_min[i])/((az_max[i]-az_min[i])/n_sss_bins)));
      }
    }
  }

  //free memory
  // for(i=0;i<nces;i++){
  //   free(az[i]);
  // }
  // free(az);
  // free(az_min);
  // free(az_max);
  // free(ces_length);

  return az_binned;
}

/*Function to find minimum of x and y*/
int miin(int x, int y){
  return y ^ ((x ^ y) & -(x < y));
}
