/** 
@file toeplitz_wizard.c
@brief easy-to-use and "all-in-one" wizard routines for the Toeplitz module
@author  Frederic Dauvergne
**  
** Project:  Midapack library, ANR MIDAS'09 - Toeplitz Algebra module
** Purpose:  Provide Toeplitz algebra tools suitable for Cosmic Microwave Background (CMB)
**           data analysis.
**
***************************************************************************
@note Copyright (c) 2010-2012 APC CNRS Université Paris Diderot
@note 
@note This program is free software; you can redistribute it and/or modify it under the terms
@note of the GNU Lesser General Public License as published by the Free Software Foundation; 
@note either version 3 of the License, or (at your option) any later version. This program is
@note distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
@note the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
@note Lesser General Public License for more details.
@note 
@note You should have received a copy of the GNU Lesser General Public License along with this
@note program; if not, see http://www.gnu.org/licenses/lgpl.html
@note For more information about ANR MIDAS'09 project see :
@note http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
@note ACKNOWLEDGMENT: This work has been supported in part by the French National Research 
@note Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
***************************************************************************
** Log: toeplitz*.c
**
** Revision 1.0b  2012/05/07  Frederic Dauvergne (APC)
** Official release 1.0beta. The first installement of the library is the Toeplitz algebra
** module.
**
** Revision 1.1b  2012/07/-  Frederic Dauvergne (APC)
** - mpi_stbmm allows now rowi-wise order per process datas and no-blocking communications.
** - OMP improvment for optimal cpu time.
** - bug fixed for OMP in the stmm_basic routine.
** - distcorrmin is used to communicate only lambda-1 datas when it is needed.
** - new reshaping routines using transformation functions in stmm. Thus, only one copy
**   at most is needed.
** - tpltz_init improvement using define_nfft and define_blocksize routines.
** - add Block struture to define each Toeplitz block.
** - add Flag structure and preprocessing parameters to define the computational strategy.
**   All the flag parameters are then available directly from the API.
**
** Revision 1.2b  2012/11/30  Frederic Dauvergne (APC)
** - extend the mpi product routine to rowwise order data distribution. This is now allowing
** tree kinds of distribution.
** - add int64 for some variables to extend the global volume of data you can use.
** - Openmp improvments.
** - Add toeplitz_wizard.c, which contains a set of easy to use routines with defined structures.
**
***************************************************************************
**
*/


#include "toeplitz.h"

/// Performs the product of a Toeplitz matrix by a general matrix either sequentially or using MPI. The complexity is hidden in the input structure, which needs to be defined by a user.
/** @ingroup wizard
 **/

int stbmmProd( Tpltz Nm1, double *V)
{

#ifdef W_MPI

  mpi_stbmm(&V, Nm1.nrow, Nm1.m_cw, Nm1.m_rw, Nm1.tpltzblocks, Nm1.nb_blocks_loc, Nm1.nb_blocks_tot, Nm1.idp, Nm1.local_V_size, Nm1.flag_stgy, Nm1.comm);

#else

//int stbmm(double **V, int nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks, int64_t idp, int local_V_size, Flag flag_stgy)

  stbmm(&V, Nm1.nrow, Nm1.m_cw, Nm1.m_rw, Nm1.tpltzblocks, Nm1.nb_blocks_loc, Nm1.idp, Nm1.local_V_size, Nm1.flag_stgy);

#endif

  return 0;
}



int gstbmmProd( Tpltz Nm1, double *V, Gap Gaps)
{

#ifdef W_MPI

  mpi_gstbmm(&V, Nm1.nrow, Nm1.m_cw, Nm1.m_rw, Nm1.tpltzblocks, Nm1.nb_blocks_loc, Nm1.nb_blocks_tot, Nm1.idp, Nm1.local_V_size, Gaps.id0gap, Gaps.lgap, Gaps.ngap, Nm1.flag_stgy, Nm1.comm);

#else

//int gstbmm0(double **V, int nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks, int64_t idp, int local_V_size, int *id0gap, int *lgap, int ngap, Flag flag_stgy)

  gstbmm(&V, Nm1.nrow, Nm1.m_cw, Nm1.m_rw, Nm1.tpltzblocks, Nm1.nb_blocks_loc, Nm1.idp, Nm1.local_V_size, Gaps.id0gap, Gaps.lgap, Gaps.ngap, Nm1.flag_stgy);

#endif

  return 0;
}


