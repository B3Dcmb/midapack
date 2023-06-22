/**
@file toeplitz.c version 1.2b, November 2012
@brief Contains the main part of the sequential routines for Toeplitz algebra
@author  Frederic Dauvergne, Maude Le Jeune, Antoine Rogier, Radek Stompor
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


#include "toeplitz.h"

#ifdef _OPENMP

#include <omp.h>

#endif

#define max(a, b)                                                              \
    ({                                                                         \
        __typeof__(a) _a = (a);                                                \
        __typeof__(b) _b = (b);                                                \
        _a > _b ? _a : _b;                                                     \
    })

#define min(a, b)                                                              \
    ({                                                                         \
        __typeof__(a) _a = (a);                                                \
        __typeof__(b) _b = (b);                                                \
        _a < _b ? _a : _b;                                                     \
    })

// r1.2 - Frederic Dauvergne (APC)
// This file contains the main part of the Toeplitz algebra module. This include
// the elementary product routines (using FFT) and initialization routines.
// This also contains the mpi version of the Toeplitz matrix product with global
// row-wise order distribution of the data.
//
// todo:
//- add in stmm non blocking communication as it is done for the stbmm routine
//- scmm_direct dont need nfft parameter


//=========================================================================
// Global parameters

/// Verbose mode
/** Prints some informative messages during the computation.
 */
int VERBOSE;
int VERBOSE_FIRSTINIT = 1;

// Parameter just to know the rank for printing when VERBOSE mode is on
int PRINT_RANK = -1;

//=========================================================================

/// Prints error message corresponding to an error number.
/** @ingroup group22
   \param error_number error number
   \param file file name
   \param line line number
*/
int print_error_message(int error_number, char const *file, int line) {
    char *str_mess;
    str_mess = (char *) malloc(100 * sizeof(char));
    if (error_number == 1)
        sprintf(str_mess,
                "Error on line %d of %s. Toeplitz band width > vector size\n",
                line, file);
    if (error_number == 2)
        sprintf(str_mess, "Error on line %d of %s. Bad allocation.\n", line,
                file);
    if (error_number == 3)
        sprintf(str_mess,
                "Error on line %d of %s. Error at fftw multithread "
                "initialization.\n",
                line, file);
    if (error_number == 7)
        sprintf(str_mess, "Error on line %d of %s.\n", line, file);
    fprintf(stderr, "%s", str_mess);
    printf("%s", str_mess);
    return error_number;
}


//=========================================================================

/// Defines an optimal size of the block used in the sliding windows algorithm.
/** @ingroup group21
    The optimal block size is computed as the minimum power of two above
   3*lambda, i.e. the smallest value equal to 2^x, where x is an integer, and
   above 3*lambda. If bs_flag is set to one, a different formula is used to
   compute the optimal block size (see MADmap: A MASSIVELY PARALLEL MAXIMUM
   LIKELIHOOD COSMIC MICROWAVE BACKGROUND MAP-MAKER, C. M. Cantalupo, J. D.
   Borrill, A. H. Jaffe, T. S. Kisner, and R. Stompor, The Astrophysical Journal
   Supplement Series, 187:212–227, 2010 March). To avoid using block size much
   bigger than the matrix, the block size is set to 3*lambda when his previous
   computed size is bigger than the matrix size n. This case append mostly for
   small matrix compared to his bandwith. \param n matrix row dimension \param
   lambda half bandwidth of the Toeplitz matrix \param bs_flag flag to use a
   different formula for optimal block size computation \param fixed_bs fixed
   blocksize value if needed
 */
int define_blocksize(int n, int lambda, int bs_flag, int fixed_bs) {
    int bs;       // computed optimal block size
    int min_bs;   // minimum block size used for the computation
    int min_pow2; // minimum power of two index used for the block size
                  // computation

    // cheating
    //  bs_flag = 5;//1;//5;
    //  fixed_bs = pow(2,15);  //2^14 winner because smaller block than 2^15 (as
    //  same speed)

    if (bs_flag == 1) {
        bs = fixed_bs;

    } else if (bs_flag
               == 2) { // this formula need to be check - seems there is a pb
        min_bs = 2 * lambda; // when bs = 2 lambda. Not enough data left in the
                             // middle
        min_pow2 = (int) ceil(log(min_bs) / log(2));
        bs       = pow(2, min_pow2);
        if (bs > n) // This is to avoid block size much bigger than the matrix.
                    // Append mostly
            bs = min_bs; // when the matrix is small compared to his bandwith

    } else if (bs_flag == 3) {
        min_bs   = 3 * lambda;
        min_pow2 = (int) ceil(log(min_bs) / log(2));
        bs       = pow(2, min_pow2);
        if (bs > n) // This is to avoid block size much bigger than the matrix.
                    // Append mostly
            bs = min_bs; // when the matrix is small compared to his bandwith

    } else if (bs_flag == 4 || bs_flag == 0) {
        min_bs   = 4 * lambda;
        min_pow2 = (int) ceil(log(min_bs) / log(2));
        bs       = pow(2, min_pow2);
        if (bs > n) // This is to avoid block size much bigger than the matrix.
                    // Append mostly
            bs = min_bs; // when the matrix is small compared to his bandwith

    } else if (bs_flag == 5) {
        // Different formula to compute the optimal block size
        bs = 1;
        while (bs < 2 * (lambda - 1) * log(bs + 1) && bs < n) { bs = bs * 2; }

    } else if (bs_flag == 6) { // the same as bs_flag==5 but with constrain on
                               // the minimal size
        // and the number of subblocks.

        min_bs   = 4 * lambda;
        min_pow2 = (int) ceil(log(min_bs) / log(2));

        min_pow2 = max(
                min_pow2,
                pow(2, 14)); // add condition to have a minimum size 2^14 for bs
        // This is based on empirical estimation and can be justified
        // by the speed benchmark of FFTW3 (see the FFTW official website).
        bs = pow(2, min_pow2);

        if (bs > n) // This is to avoid block size much bigger than the matrix.
                    // Append mostly
            bs = min_bs; // when the matrix is small compared to his bandwith

        // test if enough subblock for sliding windows algorithm:
        //     int nbloc_bs = ceil( (1.0*n)/(bs-2*distcorrmin));
        //     if (nbloc_bs<8) //Empirical condition to avoid small number of
        //     subblocks
        //       bs = 0;   //Switch to no sliding windows algorithm

    } else {
        printf("Error. Wrong value for bs_flag. Set to auto mode.\n");
        min_bs   = 4 * lambda;
        min_pow2 = (int) ceil(log(min_bs) / log(2));
        bs       = pow(2, min_pow2);
        if (bs > n) // This is to avoid block size much bigger than the matrix.
                    // Append mostly
            bs = min_bs; // when the matrix is small compared to his bandwith
    }


    if (PRINT_RANK == 0 && VERBOSE > 1)
        printf("Computed optimal blocksize is %d (with lambda = %d)\n", bs,
               lambda);

    return bs;
}


//=========================================================================

/// Defines the number of simultaneous ffts for the Toeplitz matrix product
/// computation.
/** @ingroup group21
    \param n_thread number of omp threads
    \param flag_nfft flag to set the strategy to define nfft
    \param fixed_nfft fixed nfft value if nedeed (used for the case where
   flag_nfft=1)
 */
int define_nfft(int n_thread, int flag_nfft, int fixed_nfft) {
    int nfft;

    if (flag_nfft == 0) nfft = NFFT_DEFAULT;
    else if (flag_nfft == 1)
        nfft = fixed_nfft;
    else if (flag_nfft == 2)
        nfft = n_thread;
    else {
        printf("Error. Wrong value for flag_nfft. Set to auto mode.\n");
        nfft = NFFT_DEFAULT;
    }

    return nfft;
}


//=========================================================================

/// Sets a block size and initializes all fftw arrays and plans needed for the
/// computation.
/** @ingroup group11
    Initializes the fftw arrays and plans is necessary before any computation of
   the Toeplitz matrix matrix product. Use tpltz_cleanup afterwards. \sa
   tpltz_cleanup \param n row size of the matrix used for later product \param
   lambda Toeplitz band width \param nfft maximum number of FFTs you want to
   compute at the same time \param blocksize optimal block size used in the
   sliding window algorithm to compute an optimize value) \param T_fft complex
   array used for FFTs \param T Toeplitz matrix \param V_fft complex array used
   for FFTs \param V_rfft real array used for FFTs \param plan_f fftw plan
   forward (r2c) \param plan_b fftw plan backward (c2r)
*/
int tpltz_init(int n, int lambda, int *nfft, int *blocksize,
               fftw_complex **T_fft, double *T, fftw_complex **V_fft,
               double **V_rfft, fftw_plan *plan_f, fftw_plan *plan_b,
               Flag flag_stgy) {

    // Set the VERBOSE global variable
    VERBOSE = flag_stgy.flag_verbose;


    // initialize block size
    *blocksize =
            define_blocksize(n, lambda, flag_stgy.flag_bs, flag_stgy.fixed_bs);


    // if (bs==0)
//     flag_stgy.flag_bs = 9999 //swich to noslidingwindowsalgo


// #pragma omp parallel
//{  n_thread = omp_get_num_threads(); }

//  if ((NB_OMPTHREADS <= n_thread) && (NB_OMPTHREADS != 0))
//    omp_set_num_threads(NB_OMPTHREADS);
#ifdef _OPENMP
    int n_thread = omp_get_max_threads();
#else
    int n_thread = 1;
#endif

    // initialize nfft
    *nfft = define_nfft(n_thread, flag_stgy.flag_nfft,
                        flag_stgy.fixed_nfft); //*nfft=n_thread;


    if (PRINT_RANK == 0 && VERBOSE > 0 && VERBOSE_FIRSTINIT == 1) {
        printf("Using %d threads\n", n_thread);
        printf("nfft = %d\n", *nfft);
    }

    // initialize fftw plan allocation flag
    int fftw_flag = flag_stgy.flag_fftw; // FFTW_FLAG;

    // initialize fftw for omp threads
    #ifdef fftw_MULTITHREADING
        fftw_init_omp_threads(n_thread);
    #endif

    // initialize fftw array and plan for T (and make it circulant first)
    // t1=MPI_Wtime();
    circ_init_fftw(T, (*blocksize), lambda, T_fft);
    //  t2=  MPI_Wtime();

    //  if (PRINT_RANK==0 && VERBOSE>0)
    //    printf("time circ_init_fftw=%f\n", t2-t1);

    // initialize fftw array and plan for V
    // t1=MPI_Wtime();
    rhs_init_fftw(nfft, (*blocksize), V_fft, V_rfft, plan_f, plan_b, fftw_flag);
    //  t2=  MPI_Wtime();

    //  if (PRINT_RANK==0 && VERBOSE>0)
    //    printf("time rhs_init_fftw=%f\n", t2-t1);

    if (PRINT_RANK == 0 && VERBOSE > 1)
        printf("Initialization finished successfully\n");

    VERBOSE_FIRSTINIT = 0;

    return 0;
}


//=========================================================================
#ifdef fftw_MULTITHREADING
/// Initialize omp threads for fftw plans.
/** @ingroup group21
    Initialize omp threads for fftw plans. The number of threads used for ffts
   (define by the variable n_thread) is read from OMP_NUM_THREAD environment
   variable. fftw multithreaded option is controlled by fftw_MULTITHREADING
   macro.
*/
int fftw_init_omp_threads(int fftw_n_thread) {
    int status;

    // initialize fftw omp threads
    status = fftw_init_threads();
    if (status == 0) return print_error_message(3, __FILE__, __LINE__);

    // set the number of FFTW threads
    fftw_plan_with_nthreads(fftw_n_thread);

    if (PRINT_RANK == 0 && VERBOSE > 1 && VERBOSE_FIRSTINIT == 1)
        printf("Using multithreaded FFTW with %d threads\n", fftw_n_thread);

    return 0;
}
#endif


//=========================================================================

/// Initializes fftw array and plan for the right hand side, general matrix V.
/** @ingroup group21
    Initialize fftw array and plan for the right hand side matrix V.
    \param nfft maximum number of FFTs you want to compute at the same time
    \param fft_size effective FFT size for the general matrix V (usually equal
   to blocksize) \param V_fft complex array used for FFTs \param V_rfft real
   array used for FFTs \param plan_f fftw plan forward (r2c) \param plan_b fftw
   plan backward (c2r) \param fftw_flag fftw plan allocation flag
 */
int rhs_init_fftw(const int *nfft, int fft_size, fftw_complex **V_fft,
                  double **V_rfft, fftw_plan *plan_f, fftw_plan *plan_b,
                  int fftw_flag) {
    // allocate fftw arrays and plans for V
    *V_fft  = (fftw_complex *) fftw_malloc((*nfft) * (fft_size / 2 + 1)
                                           * sizeof(fftw_complex));
    *V_rfft = (double *) fftw_malloc((*nfft) * fft_size * sizeof(double));
    if (*V_fft == 0 || *V_rfft == 0)
        return print_error_message(2, __FILE__, __LINE__);

    *plan_f = fftw_plan_many_dft_r2c(1, &fft_size, (*nfft), *V_rfft, &fft_size,
                                     1, fft_size, *V_fft, NULL, 1,
                                     fft_size / 2 + 1, fftw_flag);
    *plan_b = fftw_plan_many_dft_c2r(1, &fft_size, (*nfft), *V_fft, NULL, 1,
                                     fft_size / 2 + 1, *V_rfft, &fft_size, 1,
                                     fft_size, fftw_flag);


    return 0;
}


//=========================================================================

/// Initializes fftw array and plan for the circulant matrix T_circ obtained
/// from T.
/** @ingroup group21
    Builds the circulant matrix T_circ from T and initilizes its fftw arrays and
   plans. Use tpltz_cleanup afterwards. \sa tpltz_cleanup \param T Toeplitz
   matrix. \param fft_size effective FFT size for the circulant matrix (usually
   equal to blocksize) \param lambda Toeplitz band width. \param T_fft complex
   array used for FFTs.
 */
int circ_init_fftw(const double *T, int fft_size, int lambda,
                   fftw_complex **T_fft) {
    // routine variable
    int i;
    int circ_fftw_flag = FFTW_ESTIMATE;
    // allocation for T_fft
    *T_fft = (fftw_complex *) fftw_malloc((fft_size / 2 + 1)
                                          * sizeof(fftw_complex));
    if (*T_fft == 0) return print_error_message(2, __FILE__, __LINE__);
    double *T_circ = (double *) (*T_fft);

    // inplace fft
    fftw_plan plan_f_T;
    plan_f_T = fftw_plan_dft_r2c_1d(fft_size, T_circ, *T_fft, circ_fftw_flag);

    // make T circulant
#pragma omp parallel for
    for (i = 0; i < fft_size + 2; i++) T_circ[i] = 0.0;

    T_circ[0] = T[0];
    for (i = 1; i < lambda; i++) {
        T_circ[i]            = T[i];
        T_circ[fft_size - i] = T[i];
    }

    fftw_execute(plan_f_T);
    fftw_destroy_plan(plan_f_T);

    return 0;
}


//=========================================================================

/// Cleans fftw workspace used in the Toeplitz matrix matrix product's
/// computation.
/** @ingroup group11
    Destroy fftw plans, free memory and reset fftw workspace. \sa tpltz_init
    \param T_fft complex array used for FFTs
    \param V_fft complex array used for FFTs
    \param V_rfft real array used for FFTs
    \param plan_f fftw plan forward (r2c)
    \param plan_b fftw plan backward (c2r)
*/
int tpltz_cleanup(fftw_complex **T_fft, fftw_complex **V_fft, double **V_rfft,
                  fftw_plan *plan_f, fftw_plan *plan_b) {
    fftw_destroy_plan(*plan_f);
    fftw_destroy_plan(*plan_b);
    fftw_free(*T_fft);
    fftw_free(*V_fft);
    fftw_free(*V_rfft);
#ifdef fftw_MULTITHREADING
    fftw_cleanup_threads();
#endif
    fftw_cleanup();

    return 0;
}


//=========================================================================

/// Copies (and potentially reshapes) a selected block of the input matrix to a
/// specified position of the output matrix.
/** @ingroup group22
    Copy a matrix block of a size nblockrow x nblockcol from the input matrix
   Vin (size ninrow x nincol) starting with the element (inrow, incol) to the
   output matrix Vout (size notrow x noutcol) starting with the element (outrow,
   outcol) after multiplying by norm. If the output matrix is larger than the
   block the extra elements are either left as they were on the input or zeroed
   if zero_flag is set to 1. If the block to be copied is larger than either the
   input or the output matrix an error occurs.
*/
int copy_block(int ninrow, int nincol, double *Vin, int noutrow, int noutcol,
               double *Vout, int inrow, int incol, int nblockrow, int nblockcol,
               int outrow, int outcol, double norm, int set_zero_flag) {
    int i, j, p, offsetIn, offsetOut;

    // do some size checks first
    if ((nblockcol > nincol) || (nblockrow > ninrow) || (nblockcol > noutcol)
        || (nblockrow > noutrow)) {
        printf("Error in routine copy_block. Bad size setup.\n");
        return print_error_message(7, __FILE__, __LINE__);
    }

    if (set_zero_flag) {
#pragma omp parallel for // private(i) num_threads(NB_OMPTHREADS_CPBLOCK)
        for (i = 0; i < noutrow * noutcol;
             i++)        // could use maybe memset but how about threading
            Vout[i] = 0.0;
    }

    offsetIn  = ninrow * incol + inrow;
    offsetOut = noutrow * outcol + outrow;

    // #pragma omp parallel for private(i,j,p)
    // num_threads(NB_OMPTHREADS_CPBLOCK)
    for (i = 0; i < nblockcol * nblockrow; i++) { // copy the block
        j = i / nblockrow;
        p = i % nblockrow;
        Vout[offsetOut + j * noutrow + p] =
                Vin[offsetIn + j * ninrow + p] * norm;
    }

    return 0;
}


//=========================================================================

/// Performs the product of a circulant matrix C_fft by a matrix V_rfft using
/// fftw plans.
/** @ingroup group21
    Performs the product of a circulant matrix C_fft by a matrix V_rfft using
   fftw plans: forward - plan_f_V; and backward - plan_b_CV. C_fft is a Fourier
   (complex representation of the circulant matrix) of length fft_size/2+1;
    V_rfft is a matrix with ncol columns and fft_size rows; V_fft is a workspace
   of fft_size/2+1 complex numbers as required by the backward FFT (plan_b_CV);
   CV is the output matrix of the same size as the input V_rfft one. The FFTs
   transform ncol vectors simultanously. \param fft_size row dimension \param
   nfft number of simultaneous FFTs \param C_fft complex array used for FFTs
    \param ncol column dimension
    \param V_rfft real array used for FFTs
    \param[out] CV product of the circulant matrix C_fft by the matrix V_rfft
    \param V_fft complex array used for FFTs
    \param plan_f_V fftw plan forward (r2c)
    \param plan_b_CV fftw plan backward (c2r)
*/
int scmm_direct(int fft_size, int nfft, fftw_complex *C_fft, int ncol,
                double *V_rfft, double **CV, fftw_complex *V_fft,
                fftw_plan plan_f_V, fftw_plan plan_b_CV) {
    // routine variables
    int sizeT = fft_size / 2 + 1;
    int i, idx;

    // perform forward FFT
    fftw_execute(plan_f_V); // input in V_rfft; output in V_fft

    //  printf("ncol=%d, fft_size=%d, sizeT=%d\n", ncol, fft_size, sizeT);

    // double t1, t2;
    //   t1=MPI_Wtime();

#pragma omp parallel for private(idx) // num_threads(nfft)
    for (i = 0; i < ncol * sizeT; i++) {
        idx         = i % sizeT;
        V_fft[i][0] = C_fft[idx][0] * V_fft[i][0] - C_fft[idx][1] * V_fft[i][1];
        V_fft[i][1] = C_fft[idx][0] * V_fft[i][1] + C_fft[idx][1] * V_fft[i][0];
    }

    //  t2=  MPI_Wtime();
    //  printf("Computation time : %lf s.\n", t2-t1);


    // This is wrong :
    /*
    int icol;
    double t1, t2;
      t1=MPI_Wtime();
    #pragma omp parallel for private(i, idx)
      for(icol=0;icol<ncol;icol++) {
      for(idx=0;idx<sizeT;idx++) {
        i=icol*idx;
        V_fft[i][0] = C_fft[idx][0]*V_fft[i][0]-C_fft[idx][1]*V_fft[i][1];
        V_fft[i][1] = C_fft[idx][0]*V_fft[i][1]+C_fft[idx][1]*V_fft[i][0];
      }}
      t2=  MPI_Wtime();
    */
    //  printf("Computation time : %lf s.\n", t2-t1);


    // perform  backward FFts
    fftw_execute(plan_b_CV); // input in V_fft; output in V_rfft

    return 0;
}


//=========================================================================

/// Performs the product of a circulant matrix by a matrix using FFT's (an
/// INTERNAL routine)
/** @ingroup group21
    This routine multiplies a circulant matrix, represented by C_fft, by a
   general matrix V, and stores the output as a matrix CV. In addition the
   routine requires two workspace objects, V_fft and V_rfft, to be allocated
   prior to a call to it as well as two fftw plans: one forward (plan_f_V), and
   one backward (plan_b_TV). The sizes of the input general matrix V and the
   ouput CV are given by blocksize rows and m columns. They are stored as a
   vector in the column-wise order. The circulant matrix, which is assumed to be
   band-diagonal with a band-width lambda, is represented by a Fourier transform
   with its coefficients stored in a vector C_fft (length blocksize). blocksize
   also defines the size of the FFTs, which will be performed and therefore this
   is the value which has to be used while creating the fftw plans and
   allocating the workspaces. The latter are given as: nfft*(blocksize/2+1) for
    V_fft and nfft*blocksize for V_rfft. The fftw plans should correspond to
   doing the transforms of nfft vectors simultaneously. Typically, the
   parameters of this routine are fixed by a preceding call to Toeplitz_init().
    The parameters are :
    \param V matrix (with the convention V(i,j)=V[i+j*n])
    \param blocksize row dimension of V
    \param m column dimension of V
    \param C_fft complex array used for FFTs (FFT of the Toeplitz matrix)
    \param[out] CV product of the circulant matrix C_fft by the matrix V_rfft
    \param V_fft complex array used for FFTs
    \param V_rfft real array used for FFTs
    \param nfft number of simultaneous FFTs
    \param plan_f_V fftw plan forward (r2c)
    \param plan_b_CV fftw plan backward (c2r)
*/
int scmm_basic(double **V, int blocksize, int m, fftw_complex *C_fft,
               double **CV, fftw_complex *V_fft, double *V_rfft, int nfft,
               fftw_plan plan_f_V, fftw_plan plan_b_CV) {
    // routine variables
    int i, k;                                 // loop index
    int nloop = (int) ceil((1.0 * m) / nfft); // number of subblocks

    // Loop over set of columns
    int ncol = min(nfft, m); // a number of columns to be copied from the data
                             // to working matrix
    // equal the number of simultaneous FFTs


#pragma omp parallel for // num_threads(NB_OMPTHREADS_BASIC)//schedule(dynamic,1)
    for (i = 0; i < blocksize * ncol; i++)
        V_rfft[i] = 0.0; // could use maybe memset but how about threading


    // bug fixed conflit between num_threads and nfft
    // #pragma omp parallel for schedule(dynamic,1) num_threads(8)
    // //num_threads(nfft)
    for (k = 0; k < nloop;
         k++) {             // this is the main loop over the set of columns
        if (k == nloop - 1) // last loop ncol may be smaller than nfft
            ncol = m - (nloop - 1) * nfft;

        // init fftw matrices.
        // extracts a block of ncol full-length columns from the data matrix and
        // embeds in a bigger matrix padding each column with lambda zeros. Note
        // that all columns will be zero padded thanks to the "memset" call
        // above

        copy_block(blocksize, m, (*V), blocksize, ncol, V_rfft, 0, k * nfft,
                   blocksize, ncol, 0, 0, 1.0, 0);
        // note: all nfft vectors are transformed below ALWAYS in a single go
        // (if ncol < nfft) the extra useless work is done.

        scmm_direct(blocksize, nfft, C_fft, ncol, V_rfft, CV, V_fft, plan_f_V,
                    plan_b_CV);
        // note: the parameter CV is not really used

        // extract the relevant part from the result
        copy_block(blocksize, ncol, V_rfft, blocksize, m, (*CV), 0, 0,
                   blocksize, ncol, 0, k * nfft, 1.0 / ((double) blocksize), 0);

    } // end of loop over the column-sets


    return 0;
}


//=========================================================================

/// Performs the stand alone product of a Toeplitz matrix by a matrix using the
/// sliding window algorithm.  (an INTERNAL routine)
/** @ingroup group21
    The product is performed block-by-block with a defined block size or a
   computed optimized block size that reflects a trade off between cost of a
   single FFT of a length block_size and a number of blocks needed to perform
   the mutiplicaton. The latter determines how many spurious values are computed
   extra due to overlaps between the blocks. Use flag_offset=0 for "classic"
   algorithm and flag_offset=1 to put an offset to avoid the first and last
   lambdas terms. Usefull when a reshaping was done before with optimal column
    for a nfft. Better be inside the arguments of the routine.
    The parameters are:
    \param V \b [input] data matrix (with the convention V(i,j)=V[i+j*n]) ;
             \b [out] result of the product TV
    \param n number of rows of V
    \param m number of columns of V
    \param T Toeplitz matrix data composed of the non-zero entries of his first
   row \param T_fft complex array used for FFTs \param blocksize block size used
   in the sliding window algorithm \param lambda Toeplitz band width \param
   V_fft complex array used for FFTs \param V_rfft real array used for FFTs
    \param nfft number of simultaneous FFTs
    \param plan_f fftw plan forward (r2c)
    \param plan_b fftw plan backward (c2r)
    \param flag_offset flag to avoid extra 2*lambda padding to zeros on the
   edges \param flag_nofft flag to do product without using fft
*/
int stmm_core(double **V, int n, int m, double *T, fftw_complex *T_fft,
              int blocksize, int lambda, fftw_complex *V_fft, double *V_rfft,
              int nfft, fftw_plan plan_f, fftw_plan plan_b, int flag_offset,
              int flag_nofft) {

    // double t1,t2;

    // t1=  MPI_Wtime();

    // cheating:
    //  flag_offset = 1;

    // routine variable
    int status;
    int i, j, k, p; // loop index
    int currentsize;
    int distcorrmin = lambda - 1;

    int blocksize_eff =
            blocksize
            - 2 * distcorrmin; // just a good part after removing the overlaps
    int nbloc;                 // number of subblock of slide/overlap algorithm

    if (flag_offset == 1)
        nbloc = ceil((1.0 * (n - 2 * distcorrmin)) / blocksize_eff);
    else
        nbloc = ceil((1.0 * n)
                     / blocksize_eff); // we need n because of reshaping

    // if(PRINT_RANK==0 && VERBOSE>0)
    //   printf("nbloc=%d, n=%d, m=%d, blocksize=%d, blocksize_eff=%d\n", nbloc,
    //   n, m, blocksize, blocksize_eff);

    double *V_bloc, *TV_bloc;
    V_bloc  = (double *) calloc(blocksize * m, sizeof(double));
    TV_bloc = (double *) calloc(blocksize * m, sizeof(double));
    if ((V_bloc == 0) || (TV_bloc == 0))
        return print_error_message(2, __FILE__, __LINE__);

    int offset = 0;
    if (flag_offset == 1) offset = distcorrmin;

    int iV  = 0;      //"-distcorrmin+offset";  //first index in V
    int iTV = offset; // first index in TV

    //"k=0";
    // first subblock separately as it requires some padding. prepare the block
    // of the data vector with the overlaps on both sides
    currentsize = min(blocksize - distcorrmin + offset, n - iV);
    // note: if flag_offset=0, pad first distcorrmin elements with zeros (for
    // the first subblock only)
    //  and if flag_offset=1 there is no padding with zeros.
    copy_block(n, m, *V, blocksize, m, V_bloc, 0, 0, currentsize, m,
               distcorrmin - offset, 0, 1.0, 0);

    // do block computation
    if (flag_nofft == 1)
        status = stmm_simple_basic(&V_bloc, blocksize, m, T, lambda, &TV_bloc);
    else
        status = scmm_basic(&V_bloc, blocksize, m, T_fft, &TV_bloc, V_fft,
                            V_rfft, nfft, plan_f, plan_b);

    if (status != 0) {
        printf("Error in stmm_core.");
        return print_error_message(7, __FILE__, __LINE__);
    }


    // now copy first the new chunk of the data matrix **before** overwriting
    // the input due to overlaps !
    iV = blocksize_eff - distcorrmin + offset;

    if (nbloc > 1) {
        currentsize = min(blocksize, n - iV); // not to overshoot

        int flag_reset =
                (currentsize
                 != blocksize); // with flag_reset=1, always "memset" the block.
        copy_block(n, m, *V, blocksize, m, V_bloc, iV, 0, currentsize, m, 0, 0,
                   1.0, flag_reset);
    }

    // and now store the ouput back in V
    currentsize = min(blocksize_eff, n - iTV); // to trim the extra rows
    copy_block(blocksize, m, TV_bloc, n, m, *V, distcorrmin, 0, currentsize, m,
               iTV, 0, 1.0, 0);


    iTV += blocksize_eff;
    // now continue with all the other subblocks
    for (k = 1; k < nbloc; k++) {

        // do bloc computation
        if (flag_nofft == 1)
            status = stmm_simple_basic(&V_bloc, blocksize, m, T, lambda,
                                       &TV_bloc);
        else
            status = scmm_basic(&V_bloc, blocksize, m, T_fft, &TV_bloc, V_fft,
                                V_rfft, nfft, plan_f, plan_b);

        if (status != 0) break;


        iV += blocksize_eff;
        // copy first the next subblock to process
        if (k != nbloc - 1) {
            currentsize = min(blocksize, n - iV); // not to overshoot

            int flag_resetk =
                    (currentsize != blocksize); // with flag_reset=1, always
                                                // "memset" the block.
            copy_block(n, m, *V, blocksize, m, V_bloc, iV, 0, currentsize, m, 0,
                       0, 1.0, flag_resetk);
        }

        // and then store the output in V
        currentsize = min(blocksize_eff, n - iTV); // not to overshoot
        copy_block(blocksize, m, TV_bloc, n, m, *V, distcorrmin, 0, currentsize,
                   m, iTV, 0, 1.0, 0);
        iTV += blocksize_eff;

    } // end bloc computation


    free(V_bloc);
    free(TV_bloc);


    // t2=  MPI_Wtime();

    // if (PRINT_RANK==0 && VERBOSE>0)
    //   printf("time stmm_core=%f\n", t2-t1);

    return status;
}


//=========================================================================

/// Performs the product of a Toeplitz matrix by a general matrix using the
/// sliding window algorithm with optimize reshaping. (an INTERNAL routine)
/** @ingroup group21
    The input matrix is formatted into an optimized matrix depending on the
   block size and the number of simultaneous ffts (defined with the variable
   nfft). The obtained number of columns represent the number of vectors FFTs of
   which are computed simulatenously. The multiplication is then performed
   block-by-block with the chosen block size using the core routine. The
   parameters are : \param V \b [input] data matrix (with the convention
   V(i,j)=V[i+j*n]) ; \b [out] result of the product TV \param n number of rows
   of V \param m number of columns of V \param id0 first index of V \param l
   length of V \param T Toeplitz matrix data composed of the non-zero entries of
   his first row \param T_fft complex array used for FFTs \param lambda Toeplitz
   band width \param V_fft complex array used for FFTs \param V_rfft real array
   used for FFTs \param plan_f fftw plan forward (r2c) \param plan_b fftw plan
   backward (c2r) \param blocksize block size \param nfft number of simultaneous
   FTTs \param flag_stgy flag strategy for the product computation
*/
int stmm_main(double **V, int n, int m, int id0, int l, double *T,
              fftw_complex *T_fft, int lambda, fftw_complex *V_fft,
              double *V_rfft, fftw_plan plan_f, fftw_plan plan_b, int blocksize,
              int nfft, Flag flag_stgy) {

    // routine variable
    int i, j, k, p;                   // loop index
    int distcorrmin              = lambda - 1;
    int flag_prod_strategy_nofft = 0; // 0: ffts   1: no ffts
    int flag_shortcut_m_eff_eq_1 = 1; // 1;//1;
    int flag_shortcut_nbcol_eq_1 = 1; // 1;//1;
    int flag_nfullcol_in_middle =
            0; // 0; //in the case where m=1 can be good to direct stmm_core too
    int flag_optim_offset_for_nfft = 0;
    int flag_no_rshp               = flag_stgy.flag_no_rshp; // 0;
    int flag_nofft                 = flag_stgy.flag_nofft;   // 1;

    int m_eff = (id0 + l - 1) / n - id0 / n + 1; // number of columns
    int nfullcol;
    int nloop_middle; // change it to number of full column to improve memory

    FILE *file;
    file = stdout;

    if (l < distcorrmin) // test to avoid communications errors
        return print_error_message(1, __FILE__, __LINE__);


    // shortcut for m==1 if flag_shortcut_m_eff_eq_1==1  && nfft==1 ??
    if (m_eff == 1 && flag_shortcut_m_eff_eq_1 == 1 && nfft == 1
        || flag_no_rshp == 1 && id0 == 0 && l == n * m) {

        int flag_offset = 0;

        //  if (flag_prod_strategy_nofft==1)   //need to have T as input to make
        //  it work
        // stmm_simple_core(V, n, m, T, blocksize, lambda, nfft, flag_offset);
        //  else
        int nr = min(l, n);
        stmm_core(V, nr, m_eff, T, T_fft, blocksize, lambda, V_fft, V_rfft,
                  nfft, plan_f, plan_b, flag_offset, flag_nofft);


        return 0;
    } // End shortcut for m==1


    // the middle
    int m_middle;

    // define splitting for the product computation
    nfullcol = max(
            0, (l - (n - id0 % n) % n - (id0 + l) % n)
                       / n); // check how many full columns input data we have

    if (flag_nfullcol_in_middle == 1)
        nloop_middle = ceil(1.0 * (nfullcol) / nfft);
    else
        nloop_middle = (nfullcol) / nfft;

    if (flag_nfullcol_in_middle == 1) m_middle = nfullcol;
    else
        m_middle = nfft * nloop_middle;


    int vmiddle_size = n * m_middle;


    if (PRINT_RANK == 0 && VERBOSE > 2)
        printf("nloop_middle=%d , m_middle=%d\n", nloop_middle, m_middle);


    // compute the middle if needed
    if (nloop_middle > 0) {
        double *Vmiddle;
        int     offset_middle = (n - id0 % n) % n;
        Vmiddle               = (*V) + offset_middle;

        int flag_offset = 0;
        stmm_core(&Vmiddle, n, m_middle, T, T_fft, blocksize, lambda, V_fft,
                  V_rfft, nfft, plan_f, plan_b, flag_offset, flag_nofft);

    } //(nloop_middle>0)


    // edge  (first+last columns + extra column from the euclidian division)
    int v1edge_size = min(l, (n - id0 % n) % n);
    int v2edge_size = max(l - (v1edge_size + vmiddle_size), 0);
    int vedge_size  = v1edge_size + v2edge_size;

    // compute the edges if needed
    if (vedge_size > 0) {

        int m_v1edge, m_v2edge;
        m_v1edge   = (v1edge_size > 0) * 1; // m_v1 = 1 or 0 cannot be more
        m_v2edge   = m - (m_v1edge + m_middle);
        int  nbcol = m_v1edge + m_v2edge;
        int *nocol;
        nocol = (int *) calloc(nbcol, sizeof(double));

        // define the columns for the edge computation
        if (m_v1edge == 1) nocol[0] = 0;
        for (i = (m_v1edge); i < nbcol; i++) nocol[i] = m_middle + i;

        if (PRINT_RANK == 0 && VERBOSE > 2)
            printf("nbcol=%d , m_v1edge=%d , m_v2edge=%d\n", nbcol, m_v1edge,
                   m_v2edge);

        // shorcut for nbcol==1
        if (nbcol == 1 && nfft == 1 && flag_shortcut_nbcol_eq_1 == 1) {
            // this is the case where no reshaping is needed. This is equivalent
            // to flag_format_rshp==0
            double *Vedge;
            int offset_edge = n * nocol[0]; // work because all the previous
                                            // columns are obligatory full
            Vedge           = (*V) + offset_edge;
            int flag_offset = 0;
            stmm_core(&Vedge, vedge_size, nbcol, T, T_fft, blocksize, lambda,
                      V_fft, V_rfft, nfft, plan_f, plan_b, flag_offset,
                      flag_nofft);

        } else { // general case to compute de edges

            double *Vin;
            Vin = (*V);

            // size for the different kinds of reshaping
            int lconc = vedge_size; // another way to compute : lconc = n*nbcol
                                    // - (nocol[0]==0)*(id0%n) -
                                    // (nocol[nbcol-1]==(m-1))*(n-(id0+l)%n);
            int v1_size  = lconc + (distcorrmin) * (nbcol - 1);
            int fft_size = ceil(1.0 * v1_size / nfft) + 2 * distcorrmin;

            int flag_format_rshp = (nfft > 1) * 2 + (nfft == 1 && nbcol > 1) * 1
                                 + (nfft == 1 && nbcol == 1) * 0;
            int nrshp, mrshp, lrshp;

            define_rshp_size(flag_format_rshp, fft_size, nfft, v1_size,
                             vedge_size, &nrshp, &mrshp, &lrshp);

            // allocate Vrshp for computation
            double *Vrshp;
            Vrshp = (double *) calloc(lrshp, sizeof(double));
            double *Vout;
            Vout = (*V);

            if (PRINT_RANK == 0 && VERBOSE > 2) {
                fprintf(file, "nrshp=%d , mrshp=%d , lrshp=%d\n", nrshp, mrshp,
                        lrshp);
                fprintf(file, "flag_format_rshp=%d\n", flag_format_rshp);
            }

            build_reshape(Vin, nocol, nbcol, lconc, n, m, id0, l, lambda, nfft,
                          Vrshp, nrshp, mrshp, lrshp, flag_format_rshp);

            int flag_offset;
            if (flag_format_rshp == 2 && flag_optim_offset_for_nfft == 1)
                flag_offset = 1;
            else
                flag_offset = 0;

            // compute Vrshp
            stmm_core(&Vrshp, nrshp, mrshp, T, T_fft, blocksize, lambda, V_fft,
                      V_rfft, nfft, plan_f, plan_b, flag_offset, flag_nofft);

            extract_result(Vout, nocol, nbcol, lconc, n, m, id0, l, lambda,
                           nfft, Vrshp, nrshp, mrshp, lrshp, flag_format_rshp);


        } // End general case to compute de edges
    }     // End (vedge_size>0)


    return 0;
}


//=========================================================================
#ifdef W_MPI
/// Performs the product of a Toeplitz matrix by a general matrix using MPI. We
/// assume that the matrix has already been scattered.  (a USER routine)
/** @ingroup group12
    The multiplication is performed using FFT applied to circulant matrix in
   order to diagonalized it. The parameters are : \param V \b [input]
   distributed data matrix (with the convention V(i,j)=V[i+j*n]); \b [out]
   result of the product TV \param n number of rows of V \param m number of
   columns of V \param id0 first index of scattered V \param l length of the
   scattered V \param T Toeplitz matrix. \param lambda Toeplitz band width.
    \param flag_stgy flag strategy for the product computation
    \param comm communicator (usually MPI_COMM_WORLD)
*/
int mpi_stmm(double **V, int n, int m, int id0, int l, double *T, int lambda,
             Flag flag_stgy, MPI_Comm comm) {

    // mpi variables
    int        rank; // rank process
    int        size; // number of processes
    MPI_Status status;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);


    // routine variables
    int     i, j, k;          // some index
    int     idf    = id0 + l; // first index of scattered V for rank "rank + 1";
    int     cfirst = id0 / n; // first column index
    int     clast  = idf / n; // last column index
    int     clast_r = (idf - 1) / n;
    int     m_eff   = clast_r - cfirst + 1;
    double *V1, *Lambda;

    // Mpi communication conditions
    // Mpi comm is needed when columns are truncated
    int right   = rank + 1;
    int left    = rank - 1;
    int v1_size = l + 2 * lambda;         // size including comm
    if (rank == 0 || cfirst * n == id0) { // no left comm
        v1_size -= lambda;
        left = MPI_PROC_NULL;
    }
    if (rank == (size - 1) || clast * n == idf) { // no right comm
        v1_size -= lambda;
        right = MPI_PROC_NULL;
    }

    // init data to send
    Lambda = (double *) malloc(2 * lambda * sizeof(double));
    if (Lambda == 0) return print_error_message(2, __FILE__, __LINE__);

    for (i = 0; i < lambda; i++) {
        Lambda[i]          = (*V)[i];
        Lambda[i + lambda] = (*V)[i + l - lambda];
    }

    if (PRINT_RANK == 0 && VERBOSE > 2)
        printf("[rank %d] Left comm with %d | Right comm with %d\n", rank, left,
               right);

    // send and receive data
    MPI_Sendrecv_replace(Lambda, lambda, MPI_DOUBLE, left, MPI_USER_TAG, right,
                         MPI_USER_TAG, comm,
                         &status); // 1st comm
    MPI_Sendrecv_replace((Lambda + lambda), lambda, MPI_DOUBLE, right,
                         MPI_USER_TAG, left, MPI_USER_TAG, comm,
                         &status); // 2nd comm


    if (l < lambda) // After sendrecv to avoid problems of communication for
                    // others processors
        return print_error_message(1, __FILE__, __LINE__);

    // copy received data
    if (left == MPI_PROC_NULL && right == MPI_PROC_NULL) // 0--0 : nothing to do
        V1 = *V;
    else if (left == MPI_PROC_NULL) {                    // 0--1 : realloc
        *V = realloc(*V, v1_size * sizeof(double));
        if (*V == NULL) return print_error_message(2, __FILE__, __LINE__);
        V1 = *V;
    } else // 1--1 or 1--0 : new allocation
        V1 = (double *) malloc(v1_size * sizeof(double));

    if (left != MPI_PROC_NULL) {
        for (i = 0; i < lambda; i++) V1[i] = Lambda[i + lambda];
        id0 -= lambda;
    }
    if (right != MPI_PROC_NULL) {
        for (i = 0; i < lambda; i++) V1[i + v1_size - lambda] = Lambda[i];
    }

    // Copy input matrix V
    int offset = 0;
    if (left != MPI_PROC_NULL) {
        offset = lambda;
#pragma omp parallel for
        for (i = offset; i < l + offset; i++) V1[i] = (*V)[i - offset];
    }

    fftw_complex *V_fft, *T_fft;
    double       *V_rfft;
    fftw_plan     plan_f, plan_b;


    // Compute matrix product
    int nfft, blocksize;

    tpltz_init(v1_size, lambda, &nfft, &blocksize, &T_fft, T, &V_fft, &V_rfft,
               &plan_f, &plan_b, flag_stgy);

    if (PRINT_RANK == 0 && VERBOSE > 1)
        printf("[rank %d] Before middle-level call : blocksize=%d, nfft=%d\n",
               rank, blocksize, nfft);

    stmm_main(&V1, n, m, id0, v1_size, T, T_fft, lambda, V_fft, V_rfft, plan_f,
              plan_b, blocksize, nfft, flag_stgy);


    tpltz_cleanup(&T_fft, &V_fft, &V_rfft, &plan_f, &plan_b);

    // Copy output matrix TV
    offset = 0;
    if (left != MPI_PROC_NULL) offset = lambda;
    if (left == MPI_PROC_NULL && right == MPI_PROC_NULL) // 0--0
        *V = V1;
    else if (left == MPI_PROC_NULL) {                    // 0--1
        V1 = realloc(V1, l * sizeof(double));
        if (V1 == NULL) return print_error_message(2, __FILE__, __LINE__);
        *V = V1;
    } else { // 1--0 or 1--1
#pragma omp parallel for
        for (i = offset; i < l + offset; i++) (*V)[i - offset] = V1[i];
    }

    if (left != MPI_PROC_NULL) free(V1);

    return 0;
}

#endif
