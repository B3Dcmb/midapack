/**
 * @file  templates.h
 * @version dev
 * @authors Hamza El Bouhargani
 * @date  October 2019
 * @credit  ANR-B3DCMB
 * @description Header file for defining templates classes structure + utilities
 */

#ifndef TEMPLATES_H
#define TEMPLATES_H

/******************************************************************************/
/*                             Defining Structures                            */
/******************************************************************************/
/* Template Classes */
typedef struct TemplateClass {
int tinit;          // initial sample (local index)
int tlast;         // final (included, local index)
int nsamples;          // tlast - tinit + 1 (local samples)
int nbinMin;       // minimum bin number (in local indexes)
int nbinMax;  //   maximum bin number (in local indexes)
int nbins;   // number of bins: npixMax - npixMin + 1
int nmult;   // multiplicity: number of bins per sample
int *bins;   //   time ordered list of bin numbers (nmult bins for each sample) -> size: nmult x (tlast - tinit + 1)
double *wghts;  // time ordered list of weights (nmult weight for each sample) -> size: nmult x (tlast - tinit + 1)
// These data fields are not yet settled and are subject to change in future updates
char *flag_det;  // Detector flag
char *flag_CES;  // Constant Elevation Scan flag
char *flag_dataset;  // Data set flag
char *flag_w;  // "S": weights stored in memory, "C": weights computed on the fly
char *ID;  // Unique ID defining the Templates class type (poly of order n, ground, HWP)
int order; // temporary
} TemplateClass;

/******************************************************************************/
/*                      Prototypes of Algebra routines                        */
/******************************************************************************/
/* Projecting templates amplitudes in time domain */
int TVecProd(TemplateClass *X, int m, double *tau, double *out);

/* Projecting time domain in templates space */
int TrTVecProd(TemplateClass *X, int m, double *d, double *out);

/* Building Kernel Blocks */
int BuildKernel(TemplateClass *X, int n, double *B);

/* Inverting the Kernel Blocks */
int InvKernel(double *B, int n, double *Binv);

/******************************************************************************/
/*        Utility routines for building templates classes objects             */
/******************************************************************************/
/* Build Template Class object */
int TCinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  int nmult, int *bins, double *wghts, char *flag_det, char *flag_CES,
  char* flag_dataset, char *flag_w, char *ID);

#endif /* TEMPLATES_H */