
#ifndef   	TOEPLITZ_H_
#define   	TOEPLITZ_H_

#ifdef MPI
#include <mpi.h>
#endif

#include <fftw3.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//=========================================================================
//Basic functions definition
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
 
//=========================================================================
//Fixed parameters

/// Verbose mode
/** Print some informative messages during the computation.
 */
#ifndef VERBOSE
#define VERBOSE 0 
#endif

/// MPI user defined tag
/** Tag used for MPI communications.
 */
#ifndef MPI_USER_TAG
#define MPI_USER_TAG 123
#endif


//=========================================================================
// User routines definition

//Sequential routines (11)
int tpltz_init(int n, int lambda, int *nfft, int *blocksize, fftw_complex **T_fft, double *T, fftw_complex **V_fft, double **V_rfft, fftw_plan *plan_f, fftw_plan *plan_b);

int tpltz_cleanup(fftw_complex **T_fft, fftw_complex **V_fft, double **V_rfft,fftw_plan *plan_f, fftw_plan *plan_b);

int stmm_core(double **V, int n, int m, fftw_complex *T_fft, int blocksize, int lambda, fftw_complex *V_fft, double *V_rfft, int nfft, fftw_plan plan_f, fftw_plan plan_b, int flag_offset);

int stmm(double **V, int n, int m, int id0, int l, fftw_complex *T_fft, int lambda, fftw_complex *V_fft, double *V_rfft, fftw_plan plan_f, fftw_plan plan_b, int blocksize, int nfft);

int gstmm(double **V, int n, int m, int id0, int l, fftw_complex *T_fft, int lambda, fftw_complex *V_fft, double *V_rfft, fftw_plan plan_f, fftw_plan plan_b, int blocksize, int nfft, int *id0gap, int *lgap, int ngap);

//int stbmm(double **V, int *n, int m, int nrow, double *T, int nb_blocks_local, int nb_blocks_all, int *lambda, int *idv, int idp, int local_V_size);

//int gstbmm(double **V, int *n, int m, int nrow, double *T, int nb_blocks_local, int nb_blocks_all, int *lambda, int *idv, int id0p, int local_V_size, int *id0gap, int *lgap, int ngap);

int gap_masking(double **V, int l, int *id0gap, int *lgap, int ngap);

int gap_filling(double **V, int l, int *id0gap, int *lgap, int ngap);

int reset_gaps(double **V, int id0,int local_V_size, int m, int nrow, int *id0gap, int *lgap, int ngap);


//Mpi routines (12)
#ifdef MPI
int mpi_stmm(double **V, int n, int m, int id0, int l, double *T, int lambda, MPI_Comm comm);

int mpi_stbmm(double **V, int *n, int m, int nrow, double *T, int nb_blocks_local, int nb_blocks_all, int *lambda, int *idv, int idp, int local_V_size, MPI_Comm comm);

int mpi_gstbmm(double **V, int *n, int m, int nrow, double *T, int nb_blocks_local, int nb_blocks_all, int *lambda, int *idv, int id0p, int local_V_size, int *id0gap, int *lgap, int ngap, MPI_Comm comm);

#endif


//=========================================================================
// User routines definition

//Low level routines (21)
int optimal_blocksize(int n, int lambda, int bs_flag);

int fftw_init_omp_threads();

int rhs_init_fftw(int *nfft, int fft_size, fftw_complex **V_fft, double **V_rfft, fftw_plan *plan_f, fftw_plan *plan_b, int fftw_flag);

int circ_init_fftw(double *T, int fft_size, int lambda, fftw_complex **T_fft);

int scmm_direct(int fft_size, fftw_complex *C_fft, int ncol, double *V_rfft, double **CV, fftw_complex *V_fft, fftw_plan plan_f_V, fftw_plan plan_b_CV);

int scmm_basic(double **V, int blocksize, int m, fftw_complex *C_fft, int lambda, double **CV, fftw_complex *V_fft, double *V_rfft, int nfft, fftw_plan plan_f_V, fftw_plan plan_b_CV);

int stmm_reshape(double **V, int n, int m, int id0, int l, fftw_complex *T_fft, int lambda, fftw_complex *V_fft, double *V_rfft, fftw_plan plan_f, fftw_plan plan_b, int blocksize, int nfft);

int build_gappy_blocks(int *n, int m, int nrow, double *T, int nb_blocks_local, int nb_blocks_all, int *lambda, int *idv, int *id0gap, int *lgap, int ngap, int *nb_blocks_gappy_final, double *Tgappy, int *idvgappy, int *ngappy, int *lambdagappy, int flag_param_distmin_fixed);


//Internal routines (22)
int print_error_message(int error_number, char const *file, int line);

int copy_block(int ninrow, int nincol, double *Vin, int noutrow, int noutcol, double *Vout, int inrow, int incol, int nblockrow, int nblockcol, int outrow, int outcol, double norm, int set_zero_flag);

int vect2nfftblock(double *V1, int v1_size, double *V2, int fft_size, int nfft, int lambda);

int nfftblock2vect(double *V2, int fft_size, int nfft, int lambda, double *V1, int v1_size);

int gap_reduce(double **V, int id0, int l, int lambda, int *id0gap, int *lgap, int ngap, int *newl, int id0out);

int get_overlapping_blocks_params(int nbloc, int *idv, int *n, int local_V_size, int nrow, int idp, int *idpnew, int *local_V_size_new, int *nnew, int *ifirstBlock, int *ilastBlock);


//=========================================================================
#endif 	    /* !TOEPLITZ_H_ */


