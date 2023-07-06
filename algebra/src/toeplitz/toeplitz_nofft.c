/**
@file toeplitz_nofft.c version 1.1b, July 2012
@brief Contains basic product without using ffts for Toeplitz algebra
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
@note For more information about ANR MIDAS'09 project see :
@note http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
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

#include <math.h>
#include <stdlib.h>

#define max(a, b)               \
    ({                          \
        __typeof__(a) _a = (a); \
        __typeof__(b) _b = (b); \
        _a > _b ? _a : _b;      \
    })

#define min(a, b)               \
    ({                          \
        __typeof__(a) _a = (a);                                                \
        __typeof__(b) _b = (b);                                                \
        _a < _b ? _a : _b;                                                     \
    })

extern int PRINT_RANK;

// r1.1 - Frederic Dauvergne (APC)
// basic product without fft use.
// stmm_simple_core is not used by the API. This is similar to stmm_core by
// using a sliding windows algorithm with differents parameters.


//=========================================================================
/// Perform the product of a Toeplitz matrix by a matrix without using FFT's.
/** @ingroup group21
    This routine multiplies the values directly between them. This exploit the
   fact that the bandwith is small compared to the matrix size. The number of
   operation is then no more than (lambda*2-1) multiplications and
   (lambda*2-1)-1 additions per row.
*/
int stmm_simple_basic(double **V, int n, int m, double *T, int lambda,
                      double **TV) {

    int j_first, j_last;
    int i, j, k, Tid;
    int n_thread;
    int idx;

    int flag_nocomputeedges = 1;
    int offset_edges        = 0;

    int distcorrmin = lambda - 1;

    if (flag_nocomputeedges == 1) offset_edges = distcorrmin;


    for (k = 0; k < m; k++) {

#pragma omp parallel for shared(k, lambda, n) private(i, j, j_first, j_last,   \
                                                              Tid)
        for (i = 0 + offset_edges; i < n - offset_edges; i++) {

            (*TV)[i + k * n] = 0;
            j_first          = max(i - (lambda - 1), 0);
            j_last           = min(i + lambda, n);

            for (j = j_first; j < j_last; j++) {
                Tid = abs(j - i);
                (*TV)[i + k * n] += T[Tid] * (*V)[j + k * n];
            } // End j loop

        }     // End i loop
    }         // End k loop

    return 0;
}


//=========================================================================

/// Perform the stand alone product of a Toeplitz matrix by a matrix using the
/// sliding window algorithm.
/** @ingroup group21
    The product is performed block-by-block with a defined block size or a
   computed optimized blocksize. This routine is not used by th API. \param V \b
   [input] data matrix (with the convention V(i,j)=V[i+j*n]) ; \b [out] result
   of the product TV \param n number of rows of V \param m number of columns of
   V \param T Toeplitz matrix data composed of the non-zero entries of his first
   row \param blocksize block size used in the sliding window algorithm \param
   lambda Toeplitz band width \param nfft number of simultaneous FFTs \param
   flag_offset flag to avoid extra 2*lambda padding to zeros on the edges
*/
int stmm_simple_core(double **V, int n, int m, double *T, int blocksize,
                     int lambda, int nfft, int flag_offset) {

    // routine variable
    int status;
    int i, j, k, p; // loop index
    int currentsize;
    int distcorrmin = lambda - 1;
    int blocksize_eff =
            blocksize
            - 2 * distcorrmin; // just a good part after removing the overlaps
    int nbloc; // a number of subblock of slide/overlap algorithm

    if (flag_offset == 1)
        nbloc = ceil((1.0 * (n - 2 * distcorrmin)) / blocksize_eff);
    else
        nbloc = ceil((1.0 * n) / blocksize_eff);


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
    status = stmm_simple_basic(&V_bloc, blocksize, m, T, lambda, &TV_bloc);

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
        status = stmm_simple_basic(&V_bloc, blocksize, m, T, lambda, &TV_bloc);
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

    return status;
}
