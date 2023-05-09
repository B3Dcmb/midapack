/**
@file toeplitz_seq.c version 1.1b, July 2012
@brief Contains sequential/openMP routines for Toeplitz algebra
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
***************************************************************************
**
*/

#include "toeplitz.h"

//r1.1 - Frederic Dauvergne (APC)
//This is the sequential version of the mpi routines for the API.
//The stbmm and gstbmm are the same as the mpi_stbmm and mpi_gstbmm but without
//any communication. stmm is a simplifed call of the sequential product including
//initialization and cleaning.


//=========================================================================
/// Perform the product of a Toeplitz matrix by a general matrix using the sliding window algorithm.
/** @ingroup group11
    This is a simplifed call of the sequential product including initialization and cleaning.
    This use the flag parameters to set the comutational options.
    \param V \b [input] data matrix (with the convention V(i,j)=V[i+j*n]) ;
             \b [out] result of the product TV
    \param n number of rows of V
    \param m number of columns of V
    \param m number of columns of V
    \param T Toeplitz matrix data composed of the non-zero entries of his first row
    \param lambda number of non-zero in the first row of the Toeplitz and size of T
    \param flag_stgy flag strategy for the product computation
*/
int stmm(double **V, int n, int m, double *T, int lambda, Flag flag_stgy)
{

//fftw variables
  fftw_complex *V_fft, *T_fft;
  double *V_rfft;
  fftw_plan plan_f, plan_b;

//product parameters
  int nfft, blocksize;

  FILE *file;
  file = stdout;


  tpltz_init(n, lambda, &nfft, &blocksize, &T_fft, T, &V_fft, &V_rfft, &plan_f, &plan_b, flag_stgy);

  //Toeplitz computation
  if(VERBOSE)
    fprintf(file, "Before stmm_main call : nfft = %d, blocksize = %d\n", nfft, blocksize);
  stmm_main(V, n, m, 0, n*m, T, T_fft, lambda, V_fft, V_rfft, plan_f, plan_b, blocksize, nfft, flag_stgy);

  tpltz_cleanup(&T_fft, &V_fft, &V_rfft, &plan_f, &plan_b);


  return 0;
}


//=========================================================================
/// Performs the multiplication of a symmetric, Toeplitz block-diagonal matrix, T, by an arbitrary matrix, V, distributed over processes in the generalized column-wise way.
/** @ingroup group12
    Each process performs the multiplication sequentially for each diagonal block and based on
    the sliding window algorithm. Prior to that MPI calls are used to exchange data between
    neighboring process. Each of the diagonal blocks is a symmetric, band-diagonal Toeplitz
    matrix, which can be different for each block.
    The parameters are :
    \param V \b [input] distributed data matrix (with the convention V(i,j)=V[i+j*n]) ;
             \b [out] result of the product TV
    \param nrow number of rows of the global data matrix V
    \param m number of columns for the data matrix V in the global rowwise order
    \param m_rowwise number of columns for the data matrix V in the rowwise order per processor
    \param tpltzblocks list of the toeplitz blocks struture with its own parameters
    (idv, n, T_block, lambda) :
    \param nb_blocks number of Toeplitz blocks as stored in T
    \param idp global index of the first element of the local part of V
    \param local_V_size a number of all elements in local V
    \param flag_stgy flag strategy for the product computation
*/
int stbmm(double **V, int nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks, int64_t idp, int local_V_size, Flag flag_stgy)
{


  int nb_blocks_local=nb_blocks;
  int nb_blocks_all=nb_blocks;
 // int idp=0;
 // int local_V_size=nrow;
 MPI_Comm comm=MPI_COMM_NULL;

  mpi_stbmm(V, nrow, m_cw, m_rw, tpltzblocks, nb_blocks, nb_blocks, idp, local_V_size, flag_stgy, comm);



  return 0;
}


//====================================================================

/// Performs the multiplication of a symmetric, Toeplitz block-diagonal matrix with gaps, T, by an arbitrary matrix, V, distributed over processes.
// This matrix V contains defined gaps which represents the useless data for the comutation. The gaps indexes are defined in the global time space as the generized toeplitz matrix,
// meaning the row dimension. Each of his diagonal blocks is a symmetric, band-diagonal Toeplitz matrix, which can be different for each block.
/** @ingroup group11
    We first rebuild the Toeplitz block matrix structure to reduce the computation cost and
    skip the computations of the values on the defined gaps.

    The parameters are :
    \param V \b [input] distributed data matrix (with the convention V(i,j)=V[i+j*n]) ;
             \b [out] result of the product TV
    \param nrow number of rows of the global data matrix V
    \param m number of columns for the data matrix V in the global rowwise order
    \param m_rowwise number of columns for the data matrix V in the rowwise order per processor
    \param tpltzblocks list of the toeplitz blocks struture with its own parameters
    (idv, n, T_block, lambda) :
    - idv is the global row index defining for each Toeplitz block as stored in the vector T ;
    - n size of each Toeplitz block
    - T_block location of each Toeplitz matrix data composed of the non-zero entries of the first
    row ;
    - lambda size of each Toeplitz block data T_block. The bandwith size is then equal to lambda*2-1
    \param nb_blocks_all number of all Toeplitz block on the diagonal of the full Toeplitz matrix
    \param nb_blocks_local number of Toeplitz blocks as stored in T
    \param idp global index of the first element of the local part of V
    \param local_V_size a number of all elements in local V
    \param id0gap index of the first element of each defined gap
    \param lgap length of each defined gaps
    \param ngap number of defined gaps
    \param flag_stgy flag strategy for the product computation
*/
int gstbmm(double **V, int nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks, int64_t idp, int local_V_size, int64_t *id0gap, int *lgap, int ngap, Flag flag_stgy)
{
  int nb_blocks_local=nb_blocks;
  int nb_blocks_all=nb_blocks;
 // int idp=0;
 // int local_V_size=nrow;
 MPI_Comm comm=MPI_COMM_NULL;

  mpi_gstbmm(V, nrow, m_cw, m_rw, tpltzblocks, nb_blocks, nb_blocks, idp, local_V_size, id0gap, lgap, ngap, flag_stgy, comm);



  return 0;
}


int gstbmm0(double **V, int nrow, int m, int m_rowwise, Block *tpltzblocks, int nb_blocks_local, int nb_blocks_all, int id0p, int local_V_size, int64_t *id0gap, int *lgap, int ngap, Flag flag_stgy)
{

  int rank=0;
  int i,j,k;   //some indexes

  int flag_skip_build_gappy_blocks = flag_stgy.flag_skip_build_gappy_blocks;

  FILE *file;
  file = stdout;
  PRINT_RANK=rank ;

//put zeros at the gaps location
  reset_gaps( V, id0p, local_V_size, m, nrow, m_rowwise, id0gap, lgap, ngap);


//allocation for the gappy structure of the diagonal block Toeplitz matrix
  int nb_blocks_gappy;
  int nb_blockgappy_max;
  int Tgappysize_max;

  Block *tpltzblocks_gappy;

//some computation usefull to determine the max size possible for the gappy variables
  int Tsize=0;
  int lambdamax=0;

if (VERBOSE)
  fprintf(file, "[%d] flag_skip_build_gappy_blocks=%d\n", rank, flag_skip_build_gappy_blocks);

  if (flag_skip_build_gappy_blocks==1) {  //no build gappy blocks strategy, just put zeros at gaps location

  //compute the product using only the input Toeplitz blocks structure with zeros at the gaps location
//to remake  stbmm(V, nrow, m, m_rowwise, tpltzblocks, nb_blocks_local, nb_blocks_all, id0p, local_V_size, flag_stgy);

  }
  else { //build gappy blocks strategy

  for(Tsize=i=0;i<nb_blocks_local;i++)
    Tsize += tpltzblocks[i].lambda;

  for(i=0;i<nb_blocks_local;i++) {
    if (tpltzblocks[i].lambda>lambdamax)
      lambdamax = tpltzblocks[i].lambda;
  }

//compute max size possible for the gappy variables
  nb_blockgappy_max = nb_blocks_local+ngap;
  Tgappysize_max = Tsize + lambdamax*ngap;

//allocation of the gappy variables with max size possible
  tpltzblocks_gappy = (Block *) calloc(nb_blockgappy_max,sizeof(Block));


//build gappy Toeplitz block structure considering significant gaps locations, meaning we skip
//the gaps lower than the minimum correlation distance. You can also use the flag_param_distmin_fixed
//parameter which allows you to skip the gap lower than these value. Indeed, sometimes it's
//better to just put somes zeros than to consider two separates blocks.
//ps: This criteria could be dependant of the local lambda in futur impovements.
  int flag_param_distmin_fixed = flag_stgy.flag_param_distmin_fixed;
  build_gappy_blocks(nrow, m, tpltzblocks, nb_blocks_local, nb_blocks_all, id0gap, lgap, ngap, tpltzblocks_gappy, &nb_blocks_gappy, flag_param_distmin_fixed);


if (VERBOSE) {
    fprintf(file, "[%d] nb_blocks_gappy=%d\n", rank, nb_blocks_gappy);
    for(i=0;i<nb_blocks_gappy;i++)
      fprintf(file, "[%d] idvgappy[%d]=%ld ; ngappy[%d]=%d\n", rank, i, tpltzblocks_gappy[i].idv, i, tpltzblocks_gappy[i].n );
}
//ps: we could reallocate the gappy variables to their real size. Not sure it's worth it.

//compute the product using the freshly created gappy Toeplitz blocks structure
//to remake  stbmm(V, nrow, m, m_rowwise, tpltzblocks_gappy, nb_blocks_local, nb_blocks_all, id0p, local_V_size, flag_stgy);

  } //end flag_skip_build_gappy_blocks==1


//put zeros on V for the gaps location again. Indeed, some gaps are just replaced by zeros
//in input, it's created some fakes results we need to clear after the computation.
  reset_gaps( V, id0p, local_V_size, m, nrow, m_rowwise, id0gap, lgap, ngap);


  return 0;
}
