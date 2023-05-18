/**
@file toeplitz_params.c version 1.1b, July 2012
@brief Routines to set the flag strategy parameters for Toeplitz algebra
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


// r1.1 - Frederic Dauvergne (APC)
// This is some routines to set the flag strategy parameters

// todo:
//- add some routines to estimate the best choice strategy when automatic
// parameters are choosen.

//=========================================================================

/// Set the flag to automatic paramaters.
/** @ingroup group21
    \param flag_stgy flag strategy for the product computation
 */
int flag_stgy_init_auto(Flag *flag_stgy) {
    static const Flag z = {0};
    *flag_stgy          = z;

    flag_stgy->flag_fftw = FLAG_FFTW;
    //  flag_stgy->flag_verbose=FLAG_VERBOSE;

    return 0;
}


//=========================================================================

/// Set the flag parameters to zeros. This is almost the same as automatic.
/** @ingroup group21
    \param flag_stgy flag strategy for the product computation
 */
int flag_stgy_init_zeros(Flag *flag_stgy) {
    static const Flag z = {0};
    *flag_stgy          = z;

    return 0;
}


//=========================================================================

/// Set the parameters flag to the defined ones.
/** @ingroup group21
    \param flag_stgy flag strategy for the product computation
 */
int flag_stgy_init_defined(Flag *flag_stgy) {

    flag_stgy->flag_bs = FLAG_BS;  //0:auto 1:fixed 2:zero 3:3lambda 4:4lambda 5:formula2
    flag_stgy->flag_nfft = FLAG_NFFT;  //0:auto  1:fixed  2:numthreads  3:fftwthreads
    flag_stgy->flag_fftw = FLAG_FFTW;
    flag_stgy->flag_no_rshp = FLAG_NO_RSHP;  //0:auto  1:yes  1:no
    flag_stgy->flag_nofft = FLAG_NOFFT; //0:auto  1:yes  1:no
    flag_stgy->flag_blockingcomm = FLAG_BLOCKINGCOMM;  //0:auto 1:noblocking 2:blocking 3:blocking_nooptim
    flag_stgy->fixed_nfft = FIXED_NFFT;  //fixed init value for nfft
    flag_stgy->fixed_bs = FIXED_BS;    //fixed init value for blockside
    flag_stgy->flag_verbose = FLAG_VERBOSE;
    flag_stgy->flag_skip_build_gappy_blocks = FLAG_SKIP_BUILD_GAPPY_BLOCKS;
    flag_stgy->flag_param_distmin_fixed = FLAG_PARAM_DISTMIN_FIXED;
    flag_stgy->flag_precompute_lvl = FLAG_PRECOMPUTE_LVL;
    return 0;
}


//=========================================================================

/// Print the flag parameters values.
/** @ingroup group21
    \param flag_stgy flag strategy for the product computation
 */
int print_flag_stgy_init(Flag flag_stgy) {

    FILE *file;
    file = stdout;

    fprintf(file, "flag_bs=%d\n", flag_stgy.flag_bs);
    fprintf(file, "flag_nfft=%d\n", flag_stgy.flag_nfft);
    fprintf(file, "flag_fftw=%d\n", flag_stgy.flag_fftw);
    fprintf(file, "flag_no_rshp=%d\n", flag_stgy.flag_no_rshp);
    fprintf(file, "flag_nofft=%d\n", flag_stgy.flag_nofft);
    fprintf(file, "flag_blockingcomm=%d\n", flag_stgy.flag_blockingcomm);
    fprintf(file, "fixed_nfft=%d\n", flag_stgy.fixed_nfft);
    fprintf(file, "fixed_bs=%d\n", flag_stgy.fixed_bs);
    fprintf(file, "flag_verbose=%d\n", flag_stgy.flag_verbose);
    fprintf(file, "flag_skip_build_gappy_blocks=%d\n",
            flag_stgy.flag_skip_build_gappy_blocks);
    fprintf(file, "flag_param_distmin_fixed=%d\n",
            flag_stgy.flag_param_distmin_fixed);
    fprintf(file, "flag_precompute_lvl=%d\n", flag_stgy.flag_precompute_lvl);

    return 0;
}
