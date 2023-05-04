/**
 * @file  templates.h
 * @version dev
 * @authors Hamza El Bouhargani
 * @date  October 2019
 * @credit  ANR-B3DCMB
 * @description Header file for defining templates classes structure + utilities
 */

/******************************************************************************/
/*                             Defining Structures                            */
/******************************************************************************/
/* Metadata used by template operators */
typedef struct TemplateMetadata {
  int nces; // local number of CESs
  double sampling_rate;
  int *scan_size; // Array of local CES lengths (in samples)
  int **sweeptstamps; // 2D array of time stamps corresponding to scan direction changes (one row for each CES)
  int **az_binned; // 2D array of boresight azimuth bins (one row for each CES)
  double ***hwp_mod // 3D Array of cosines and sines of the HWP angle (1d: cos/sin, 2d: CESs, 3d: time)
  int nhwp; // number of HWPSS templates
  double delta_t; //Baseline length for fitting constant amplitude HWPSS templates
}
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

/* Template list */
typedef struct Templates {
  TemplateClass *TemplatesList; // Pointer to list of template classes
  TemplateMetadata mdata; // Contains scan information useful for building and applyin template filters
  int nclass; // number of template classes
  int nb_templates_loc; // local number of templates
  int store_hwp; // Flag to store all HWPSS orders once before single template expansion  
} Templates;

typedef struct hwpss_w {
  int ces_id; // CES local index
  double **hwpcos; // array of cosines harmonics in the HWPSS: hwpcos[order-1][t]
  double **hwpsin; // array of sines harmonics in the HWPSS: hwpsin[order-1][t]
} hwpss_w;

// typedef struct poly_w {
//   int ces_id; // CES local index
//   double **poly; // array of Legenre polynomials: poly[order][t]
// } poly_w;

/******************************************************************************/
/*                      Prototypes of Algebra routines                        */
/******************************************************************************/
/* Projecting templates amplitudes in time domain */
int TVecProd(TemplateClass *X, int nces, int m, double sampling_freq, int *scan_size,
  int **sweeptstamps, int **az_binned, double ***hwp_mod, int nhwp, double delta_t,
  int store_hwp, double *tau, double *out);

/* Projecting time domain in templates space */
int TrTVecProd(TemplateClass *X, int nces, int m, double sampling_freq, int *scan_size,
  int **sweeptstamps, int **az_binned, double ***hwp_mod, int nhwp, double delta_t,
  int store_hwp, double *d, double *out);

/* Building Kernel Blocks */
int BuildKernel(TemplateClass *X, int n, double *B, double w, int *sweeptstamps,
  int *az_binned, hwpss_w hwpss_wghts, double delta_t, double sampling_freq);

/* Inverting the Kernel Blocks */
int InvKernel(double *B, int n, double *Binv);

void transpose_nn(double *A, int n);
int inverse_svd(int m, int n, int lda,  double *a);

double P0(double x);

double P1(double x);

double Pn(double x, int n);

double Legendre(double x, double a, double b, int n);

void build_hwpss_w(hwpss_w *hwpss_wghts, double **hwp_mod, int size, int order,
  int ces_id);

void free_hwpss_w(hwpss_w *hwpss_wghts, int order);
/******************************************************************************/
/*        Utility routines for building templates classes objects             */
/******************************************************************************/
int Tlist_init(TemplateClass *X, int ndet, int nces, int *block_nsamples, int **detnsweeps,
  int *scan_size, int **sweeptstamps, int n_sss_bins, int **az_binned, double sampling_freq,
  int npoly, int ground, int nhwp, double delta_t, int store_hwp, double ***hwp_mod);

int Polyinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  int *sweeptstamps, double sampling_freq, int order);

int SSSinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  int *az_binned);

int HWPSSinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
    double delta_t, double sampling_freq, int order, hwpss_w hwpss_wghts);

/* Build Template Class object */
int TCinit(TemplateClass *X, int tinit, int tlast, int nbinMin, int nbinMax,
  int nmult, char *flag_det, char *flag_CES, char* flag_dataset, char *flag_w,
  char *ID);

int expandpolydata(TemplateClass *X ,int *bins, double *wghts,
  int *sweeptstamps, double sampling_freq, int order);

int expandSSSdata(TemplateClass *X, int *bins, int *az_binned);

int expandHWPSSdata(TemplateClass *X, int *bins0, double **wghts0, int *bins1,
  double **wghts1, double delta_t, double sampling_freq, int order,
  hwpss_w hwpss_wghts);

int** bin_az(double **az, double *az_min, double *az_max, int *ces_length,
  int sss, int n_sss_bins, int nces);

int miin(int x, int y);
