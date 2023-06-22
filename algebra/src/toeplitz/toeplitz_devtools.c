/**
@file toeplitz_devtools.c version 1.1b, July 2012
@brief Contains developpement tools routines for Toeplitz algebra
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
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 3 of the License, or
@note (at your option) any later version. This program is distributed in the
hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of
@note MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
General Public License for more details.
@note
@note You should have received a copy of the GNU Lesser General Public License
along with this program; if not, see http://www.gnu.org/licenses/lgpl.html
@note For more information about ANR MIDAS'09 project see
http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
@note ACKNOWLEDGMENT: This work has been supported in part by the French
National Research Agency (ANR) through COSINUS program (project MIDAS no.
ANR-09-COSI-009).
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
***************************************************************************
**
*/

#include "toeplitz.h"
#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//
// dev tools for cblas and print - fd@apc

int stmm_cblas(int n_loc, int m_loc, double *T2_loc, double *V, double *TV2) {

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n_loc, m_loc, n_loc,
                1, T2_loc, n_loc, (V), n_loc, 1, TV2, n_loc);

    return 0;
}


// Build full Toeplitz matrix needed for cblas computation
int build_full_Toeplitz(int n_loc, double *T_loc, int lambda_loc,
                        double *T2_loc) {

    int i, j;

    for (j = 0; j < n_loc; j++) { // reset all the matrix to zeros
        for (i = 0; i < n_loc; i++) { T2_loc[j * n_loc + i] = 0; }
    }

    for (j = 0; j < n_loc;
         j++) { // Full Toeplitz matrix needed for cblas computation
        for (i = 0; i < lambda_loc; i++) {
            if (j - i >= 0) T2_loc[j * n_loc + j - i] = T_loc[i];
            if (j + i < n_loc) T2_loc[j * n_loc + j + i] = T_loc[i];
        }
    }


    return 0;
}


int print_full_Toeplitz(int n_loc, double *T2_loc) {

    int i, j;

    FILE *file;
    file = stdout;

    for (i = 0; i < n_loc; i++) {
        for (j = 0; j < n_loc; j++) {
            fprintf(file, "%.1f\t", T2_loc[i + j * n_loc]);
        }
        fprintf(file, "\n");
    }


    return 0;
}


int print_full_matrix(int n_loc, int m_loc, double *Mat) {

    int i, j;

    FILE *file;
    file = stdout;

    for (i = 0; i < n_loc; i++) {
        for (j = 0; j < m_loc; j++) {
            fprintf(file, "%.1f\t", Mat[i + j * n_loc]);
        }
        fprintf(file, "\n");
    }


    return 0;
}