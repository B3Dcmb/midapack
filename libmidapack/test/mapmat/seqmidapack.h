/*
** File: seqmidapack.h, version 1.1b, July 2012  
** This file is the main part of the Toeplitz algebra module 
**
** Project:  Midapack library, ANR MIDAS'09 
** Purpose:  Provide algebra tools suitable for Cosmic Microwave Background (CMB)
**           data analysis.
**
** Authors:  Pierre Cargemel, Frederic Dauvergne, Maude Le Jeune, Antoine Rogier, 
**           Radek Stompor (APC, Paris)
**
**
** 
** Copyright (c) 2010-2012 APC CNRS Université Paris Diderot
** 
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation; either version 3 of the License, or
** (at your option) any later version.
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
** GNU Lesser General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, see http://www.gnu.org/licenses/lgpl.html
**
** For more information about ANR MIDAS'09 project see
** http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
** 
** ACKNOWLEDGMENT: This work has been supported in part by the French National 
** Research Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
**
**
*/


//=========================================================================
//=======================  Toeplitz algebra module ========================
//=========================================================================

//Basic functions definition
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))


//=========================================================================
//Fixed parameters


/// MPI user defined tag
/** Tag used for MPI communications.
 */
#ifndef MPI_USER_TAG
#define MPI_USER_TAG 123
#endif

//Define this parameter to use fftw multithreading
//This is not fully tested
//#ifndef fftw_MULTITHREADING
//#define fftw_MULTITHREADING 
//#endif



/// Number of FFT's performed at the same time (default value). 
/** FFT's can be performed simultaneously using advanced fftw plans.
*/
#ifndef NFFT_DEFAULT
#define NFFT_DEFAULT 1
#endif


/// fftw plan allocation flag (default value). 
/** fftw plan allocation flag can be one of (from fastest to lowest):
    ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE. Default is MEASURE.
    Fastest plans lead to sub optimal FFT computation. 
*/
#ifndef FFTW_FLAG_AUTO
#define FFTW_FLAG_AUTO FFTW_ESTIMATE
#endif


//Parameters to define the computational strategy
#ifndef FLAG_STGY
#define FLAG_STGY

#define FLAG_BS 0  //0:auto  1:fixed  2:zero  3:3lambda  4:4lambda  5:formula2 
#define FLAG_NFFT 0  //0:auto  1:fixed  2:numthreads  3:fftwthreads 
#define FLAG_FFTW FFTW_FLAG_AUTO  //ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE. Default is MEASURE
#define FLAG_NO_RSHP 0  //0:auto  1:yes  1:no
#define FLAG_NOFFT 0 //0:auto  1:yes  1:no
#define FLAG_BLOCKINGCOMM 0  //0:auto  1:noblocking  2:blocking 
#define FIXED_NFFT 0  //fixed init value for nfft
#define FIXED_BS 0   //fixed init value for blockside
#define FLAG_VERBOSE 0
#define FLAG_SKIP_BUILD_GAPPY_BLOCKS 0
#define FLAG_PARAM_DISTMIN_FIXED 0

#endif

//=========================================================================
//Global parameters

extern int VERBOSE;
extern int PRINT_RANK;

//=========================================================================
//Strutures definition

typedef struct Block {
     int idv;
     double *T_block;  //pointer of the Toeplitz data
     int lambda;
     int n;
} Block;


typedef struct Flag {
       int flag_bs;  //bs used formula
       int flag_nfft;
       int flag_fftw;
       int flag_no_rshp;  //with or without
       int flag_nofft;
       int flag_blockingcomm;
       int fixed_nfft;  //init value for nfft
       int fixed_bs;
       int flag_verbose;
       int flag_skip_build_gappy_blocks;
       int flag_param_distmin_fixed;
} Flag;

typedef struct Gap {
  int id0;
  int l;
} Gap;


//=========================================================================
// User routines definition (API)

//Sequential routines (group 11)
int tpltz_init(int n, int lambda, int *nfft, int *blocksize, fftw_complex **T_fft, double *T, fftw_complex **V_fft, double **V_rfft, fftw_plan *plan_f, fftw_plan *plan_b, Flag flag_stgy);

int tpltz_cleanup(fftw_complex **T_fft, fftw_complex **V_fft, double **V_rfft,fftw_plan *plan_f, fftw_plan *plan_b);

int stmm_core(double **V, int n, int m, double *T, fftw_complex *T_fft, int blocksize, int lambda, fftw_complex *V_fft, double *V_rfft, int nfft, fftw_plan plan_f, fftw_plan plan_b, int flag_offset, int flag_nofft);

int stmm_main(double **V, int n, int m, int id0, int l, double *T, fftw_complex *T_fft, int lambda, fftw_complex *V_fft, double *V_rfft, fftw_plan plan_f, fftw_plan plan_b, int blocksize, int nfft, Flag flag_stgy);

int stmm(double **V, int n, int m, double *T, int lambda, Flag flag_stgy);

int stbmm(double **V, int nrow, int m, int m_rowwise, Block *tpltzblocks, int nb_blocks_local, int nb_blocks_all, int idp, int local_V_size, Flag flag_stgy);

int gstbmm(double **V, int nrow, int m, int m_rowwise, Block *tpltzblocks, int nb_blocks_local, int nb_blocks_all, int id0p, int local_V_size, int *id0gap, int *lgap, int ngap,Flag flag_stgy);

int reset_gaps(double **V, int id0,int local_V_size, int m, int nrow, int m_rowwise, int *id0gap, int *lgap, int ngap);



//=========================================================================
// User routines definition

//Low level routines (group 21)
int flag_stgy_init_auto(Flag *flag_stgy);

int flag_stgy_init_zeros(Flag *flag_stgy);

int flag_stgy_init_defined(Flag *flag_stgy);

int print_flag_stgy_init(Flag flag_stgy);

int define_blocksize(int n, int lambda, int bs_flag, int fixed_bs);

int define_nfft(int n_thread, int flag_nfft, int fixed_nfft);

int scmm_direct(int fft_size, int nfft, fftw_complex *C_fft, int ncol, double *V_rfft, double **CV, fftw_complex *V_fft, fftw_plan plan_f_V, fftw_plan plan_b_CV);

int scmm_basic(double **V, int blocksize, int m, fftw_complex *C_fft, double **CV, fftw_complex *V_fft, double *V_rfft, int nfft, fftw_plan plan_f_V, fftw_plan plan_b_CV);

int stmm_simple_basic(double **V, int n, int m, double *T, int lambda, double **TV);

int build_gappy_blocks(int nrow, int m, Block *tpltzblocks, int nb_blocks_local, int nb_blocks_all, int *id0gap, int *lgap, int ngap, Block *tpltzblocks_gappy, int *nb_blocks_gappy_final, int flag_param_distmin_fixed);



//=========================================================================
//============================  Mapmat module =============================
//=========================================================================

//Parameters

#define NONE 0
#define RING 1
#define BUTTERFLY 2
#define NONBLOCKING 3
#define NOEMPTY 4 
#define SEQ 0 
#define OMP 1
#define GPU 2


//=========================================================================
//Strutures definition
 
typedef struct { 
  int		flag;			
  int		m;			 
  int		nnz;                     
  int		*indices;		
  double	*values;		
  int		lcount;
  int		*lindices;		  
}Mat;


//=========================================================================
// User routines definition (API)

int MatInit(Mat *A, int m, int nnz, int *indices, double *values);

int MatAllocate(Mat *A, int m, int nnz);

int MatFree(Mat *A);						

int MatSetIndices(Mat *A, int count, int offset, int stride, int *indices);

int MatSetValues(Mat *A, int count, int offset, int stride, double *values);

int MatLocalShape(Mat *A, int sflag, int pflag);

int MatVecProd(Mat *A, double *x, double *y, int pflag);

int TrMatVecProd(Mat *A, double *y, double* x, int pflag);

int MatLoad(Mat *A, char *filename);

int MatSave(Mat *A, char* filename);

int MatInfo(Mat *A);							



