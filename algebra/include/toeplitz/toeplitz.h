/**
@file: toeplitz.h version 1.2b, November 2012
@brief Header file with main definitions and declarations for the Toeplitz
algebra module
@author  Frederic Dauvergne
**
** Project:  Midapack library, ANR MIDAS'09 - Toeplitz Algebra module
** Purpose:  Provide Toeplitz algebra tools suitable for Cosmic Microwave
Background (CMB)
**           data analysis.
**
***************************************************************************
@note Copyright (c) 2010-2012 APC CNRS Université Paris Diderot
@note
@note This program is free software; you can redistribute it and/or modify it
under the terms
@note of the GNU Lesser General Public License as published by the Free Software
Foundation;
@note either version 3 of the License, or (at your option) any later version.
This program is
@note distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even
@note the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU
@note Lesser General Public License for more details.
@note
@note You should have received a copy of the GNU Lesser General Public License
along with this
@note program; if not, see http://www.gnu.org/licenses/lgpl.html
@note
@note For more information about ANR MIDAS'09 project see :
@note http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
@note
@note ACKNOWLEDGMENT: This work has been supported in part by the French
National Research
@note Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
***************************************************************************
** Log: toeplitz*.c
**
** Revision 1.0b  2012/05/07  Frederic Dauvergne (APC)
** Official release 1.0beta. The first installement of the library is the
Toeplitz algebra
** module.
**
** Revision 1.1b  2012/07/-  Frederic Dauvergne (APC)
** - mpi_stbmm allows now rowi-wise order per process datas and no-blocking
communications.
** - OMP improvment for optimal cpu time.
** - bug fixed for OMP in the stmm_basic routine.
** - distcorrmin is used to communicate only lambda-1 datas when it is needed.
** - new reshaping routines using transformation functions in stmm. Thus, only
one copy
**   at most is needed.
** - tpltz_init improvement using define_nfft and define_blocksize routines.
** - add Block struture to define each Toeplitz block.
** - add Flag structure and preprocessing parameters to define the computational
strategy.
**   All the flag parameters are then available directly from the API.
**
** Revision 1.2b  2012/11/30  Frederic Dauvergne (APC)
** - extend the mpi product routine to rowwise order data distribution. This is
now allowing
** tree kinds of distribution.
** - add int64 for some variables to extend the global volume of data you can
use.
** - Openmp improvments.
** - Add toeplitz_wizard.c, which contains a set of easy to use routines with
defined structures.
**
***************************************************************************
**
*/

#ifndef TOEPLITZ_H
#define TOEPLITZ_H

#ifdef W_MPI
#include <mpi.h>
#endif

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

//=========================================================================
// Fixed parameters

/// MPI user defined tag
/** Tag used for MPI communications.
 */
#ifndef MPI_USER_TAG
#define MPI_USER_TAG 123
#endif

// Define this parameter to use fftw multithreading
// This is not fully tested
// #ifndef fftw_MULTITHREADING
// #define fftw_MULTITHREADING
// #endif

/// Number of FFT's performed at the same time (default value).
/** FFT's can be performed simultaneously using advanced fftw plans.
 */
#ifndef NFFT_DEFAULT
#define NFFT_DEFAULT 1 /*1*/
#endif

/// fftw plan allocation flag (default value).
/** fftw plan allocation flag can be one of (from fastest to lowest):
    ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE. Default is MEASURE.
    Fastest plans lead to sub optimal FFT computation.
*/
#ifndef FFTW_FLAG_AUTO
#define FFTW_FLAG_AUTO FFTW_ESTIMATE
#endif

// Parameters to define the computational strategy
#ifndef FLAG_STGY
#define FLAG_STGY

#define FLAG_BS 0   // 0:auto  1:fixed  2:zero  3:3lambda  4:4lambda  5:formula2
#define FLAG_NFFT 0 // 0:auto  1:fixed  2:numthreads  3:fftwthreads
#define FLAG_FFTW                                                              \
    FFTW_FLAG_AUTO // ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE. Default is MEASURE
#define FLAG_NO_RSHP 0      // 0:auto  1:yes  1:no
#define FLAG_NOFFT 0        // 0:auto  1:yes  1:no
#define FLAG_BLOCKINGCOMM 0 // 0:auto  1:noblocking  2:blocking
#define FIXED_NFFT 0        // fixed init value for nfft
#define FIXED_BS 0          // fixed init value for blockside
#define FLAG_VERBOSE 0
#define FLAG_SKIP_BUILD_GAPPY_BLOCKS 0
#define FLAG_PARAM_DISTMIN_FIXED 0
#define FLAG_PRECOMPUTE_LVL                                                    \
    0 // 0: no precompute  1: precompute plans  2: precomputes Toeplitz and
      // plans

#endif

//=========================================================================
// Global parameters

extern int VERBOSE;
extern int PRINT_RANK;

//=========================================================================
// Strutures definition

typedef struct tpltz_block_t {
    int64_t idv;
    double *T_block; // pointer of the Toeplitz data
    int lambda;
    int n;
    int first_gap; // index of the first relevant local gap
    int last_gap;  // index of the last relevant local gap
    /* For precomputed fftw
         int bs;
         int nfft;
         fftw_complex *T_fft;
         fftw_complex *V_fft;
         double *V_rfft;
         fftw_plan plan_f;
         fftw_plan plan_b;
         fftw_plan plan_f_T;
    */
} Block;

typedef struct tpltz_flag_t {
    int flag_bs; // bs used formula
    int flag_nfft;
    int flag_fftw;
    int flag_no_rshp; // with or without
    int flag_nofft;
    int flag_blockingcomm;
    int fixed_nfft; // init value for nfft
    int fixed_bs;   // long long int
    int flag_verbose;
    int flag_skip_build_gappy_blocks;
    int flag_param_distmin_fixed;
    int flag_precompute_lvl;
} Flag;

typedef struct tpltz_gap_t {
    int64_t *id0gap;
    int *lgap;
    int ngap;
} Gap;

typedef struct tpltz_mat_t {
    int64_t nrow; // n total
    int m_cw;     // V column number in the linear row-wise order (vect row-wise
                  // order)
    int m_rw; // V column number in the uniform row-wise order (matrix row-wise
              // order)
    Block *tpltzblocks;
    int nb_blocks_loc;
    int nb_blocks_tot;
    int64_t idp;
    int local_V_size;
    Flag flag_stgy;
#ifdef W_MPI
    MPI_Comm comm;
#endif
} Tpltz;

//=========================================================================
// Groups definition for documentation

/** @defgroup toeplitz TOEPLITZ module
 *  Toeplitz matrix algebra module
 */

/** @defgroup group1 user interface (API)
 *  These routines provide main functionality of the Toeplitz algebra library.
 *  They are divided in two groups:
 *  - shared-memory: multithreaded (openMP/sequential) routines
 *  - distributed-memory (MPI) routines
 *  @ingroup toeplitz
 */

/** @defgroup wizard wizard routines
 *  These are easy to use drivers for the Toeplitz routines.
 * @ingroup group1
 */

/** @defgroup group11 multithreaded/sequential routines
 *  These are shared-memory routines.
 *  @ingroup group1
 */

/** @defgroup group12 distributed memory (MPI) routines
 *  These are distributed-memory routines.
 *  @ingroup group1
 */

/** @defgroup group2 internal routines
 *  These are auxiliary, internal routines, not intended to be used by no-expert
 * user. They are divided in two groups:
 *  - low level routines
 *  - internal routines
 *  @ingroup toeplitz
 */

/** @defgroup group21 low-level routines
 *  These are low-level routines.
 *  @ingroup group2
 */

/** @defgroup group22 lower internal routines
 *  These are lower internal routines.
 *  @ingroup group2
 */

//=========================================================================
// User routines definition (API)

// Sequential routines (group 11)
int tpltz_init(int n, int lambda, int *nfft, int *blocksize,
               fftw_complex **T_fft, double *T, fftw_complex **V_fft,
               double **V_rfft, fftw_plan *plan_f, fftw_plan *plan_b,
               Flag flag_stgy);

int tpltz_cleanup(fftw_complex **T_fft, fftw_complex **V_fft, double **V_rfft,
                  fftw_plan *plan_f, fftw_plan *plan_b);

int stmm_core(double **V, int n, int m, double *T, fftw_complex *T_fft,
              int blocksize, int lambda, fftw_complex *V_fft, double *V_rfft,
              int nfft, fftw_plan plan_f, fftw_plan plan_b, int flag_offset,
              int flag_nofft);

int stmm_main(double **V, int n, int m, int id0, int l, double *T,
              fftw_complex *T_fft, int lambda, fftw_complex *V_fft,
              double *V_rfft, fftw_plan plan_f, fftw_plan plan_b, int blocksize,
              int nfft, Flag flag_stgy);

int stmm(double **V, int n, int m, double *T, int lambda, Flag flag_stgy);

int stbmm(double **V, int nrow, int m_cw, int m_rw, Block *tpltzblocks,
          int nb_blocks, int64_t idp, int local_V_size, Flag flag_stgy);

int gstbmm(double **V, int nrow, int m_cw, int m_rw, Block *tpltzblocks,
           int nb_blocks, int64_t idp, int local_V_size, int64_t *id0gap,
           int *lgap, int ngap, Flag flag_stgy);

int reset_gaps(double **V, int64_t id0, int local_V_size, int m, int64_t nrow,
               int m_rowwise, const int64_t *id0gap, const int *lgap, int ngap);

// Mpi routines (group 12)
#ifdef W_MPI

int mpi_stmm(double **V, int n, int m, int id0, int l, double *T, int lambda,
             Flag flag_stgy, MPI_Comm comm);

int mpi_stbmm(double **V, int64_t nrow, int m, int m_rowwise,
              Block *tpltzblocks, int nb_blocks_local, int nb_blocks_all,
              int64_t idp, int local_V_size, Flag flag_stgy, MPI_Comm comm);

int mpi_gstbmm(double **V, int nrow, int m, int m_rowwise, Block *tpltzblocks,
               int nb_blocks_local, int nb_blocks_all, int id0p,
               int local_V_size, int64_t *id0gap, int *lgap, int ngap,
               Flag flag_stgy, MPI_Comm comm);

#endif

//=========================================================================
// User routines definition

// Low level routines (group 21)
int flag_stgy_init_auto(Flag *flag_stgy);

int flag_stgy_init_zeros(Flag *flag_stgy);

int flag_stgy_init_defined(Flag *flag_stgy);

int print_flag_stgy_init(Flag flag_stgy);

int define_blocksize(int n, int lambda, int bs_flag, int fixed_bs);

int define_nfft(int n_thread, int flag_nfft, int fixed_nfft);

// int fftw_init_omp_threads(int fftw_n_thread);

int rhs_init_fftw(const int *nfft, int fft_size, fftw_complex **V_fft,
                  double **V_rfft, fftw_plan *plan_f, fftw_plan *plan_b,
                  int fftw_flag);

int circ_init_fftw(const double *T, int fft_size, int lambda,
                   fftw_complex **T_fft);

int scmm_direct(int fft_size, int nfft, fftw_complex *C_fft, int ncol,
                double *V_rfft, double **CV, fftw_complex *V_fft,
                fftw_plan plan_f_V, fftw_plan plan_b_CV);

int scmm_basic(double **V, int blocksize, int m, fftw_complex *C_fft,
               double **CV, fftw_complex *V_fft, double *V_rfft, int nfft,
               fftw_plan plan_f_V, fftw_plan plan_b_CV);

int stmm_simple_basic(double **V, int n, int m, double *T, int lambda,
                      double **TV);

int build_gappy_blocks(int nrow, int m, Block *tpltzblocks, int nb_blocks_local,
                       int nb_blocks_all, const int64_t *id0gap,
                       const int *lgap, int ngap, Block *tpltzblocks_gappy,
                       int *nb_blocks_gappy_final,
                       int flag_param_distmin_fixed);

// Internal routines (group 22)
int print_error_message(int error_number, char const *file, int line);

int copy_block(int ninrow, int nincol, double *Vin, int noutrow, int noutcol,
               double *Vout, int inrow, int incol, int nblockrow, int nblockcol,
               int outrow, int outcol, double norm, int set_zero_flag);

int vect2nfftblock(double *V1, int v1_size, double *V2, int fft_size, int nfft,
                   int lambda);

int nfftblock2vect(double *V2, int fft_size, int nfft, int lambda, double *V1,
                   int v1_size);

int get_overlapping_blocks_params(int nbloc, Block *tpltzblocks,
                                  int local_V_size, int64_t nrow, int64_t idp,
                                  int64_t *idpnew, int *local_V_size_new,
                                  int *nnew, int *ifirstBlock, int *ilastBlock);

// Wizard routines
int stbmmProd(Tpltz *Nm1, double *V);

int gstbmmProd(Tpltz *Nm1, double *V, Gap *Gaps);

// Toeplitz rshp routines
int fctid_mat2vect(int i, int id0, int n, int lambda);

int fctid_mat2vect_inv(int i, int id0, int n, int lambda);

int fctid_concatcol(int i, int id0, int n, int m, int l, int lconc, int lambda,
                    int *nocol, int nbcol);

int fctid_concatcol_inv(int i, int id0, int n, int m, int l, int lconc,
                        int lambda, int *nocol_inv, int nbcol);

int fctid_vect2nfftblock(int i, int v1_size, int fft_size, int nfft,
                         int lambda);

int is_needconcat(int *nocol, int nbcol);

int fctid_vect2nfftblock_inv(int i, int v1_size, int fft_size, int nfft,
                             int lambda);

int define_rshp_size(int flag_format_rshp, int fft_size, int nfft, int v1_size,
                     int vedge_size, int *nrshp, int *mrshp, int *lrshp);

int build_nocol_inv(int *nocol, int nbcol, int m);

int build_reshape(double *Vin, int *nocol, int nbcol, int lconc, int n, int m,
                  int id0, int l, int lambda, int nfft, double *Vrshp,
                  int nrshp, int mrshp, int lrshp, int flag_format_rshp);

int extract_result(double *Vout, int *nocol, int nbcol, int lconc, int n, int m,
                   int id0, int l, int lambda, int nfft, double *Vrshp,
                   int nrshp, int mrshp, int lrshp, int flag_format_rshp);

#ifdef __cplusplus
}
#endif

//=========================================================================
#endif // TOEPLITZ_H
