/**
@file toeplitz_utils.c version 1.2b, July 2012
@brief Contains a set of utilitaries routines for Toeplitz algebra
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
#include <time.h>
extern int PRINT_RANK;

//set of utilitaries routines - fd@apc


int defineTpltz( Tpltz *Nm1, int64_t nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks_loc, int nb_blocks_tot, int64_t idp, int local_V_size, Flag flag_stgy, MPI_Comm comm)
{

  Nm1->nrow = nrow; //glob //recup du fichier params apres (en variables globales)
  Nm1->m_cw = m_cw; //glob
  Nm1->m_rw = m_rw; //glob
  Nm1->tpltzblocks = tpltzblocks; //toep
  Nm1->nb_blocks_loc = nb_blocks_loc; //toep
  Nm1->nb_blocks_tot = nb_blocks_tot;  //toep
  Nm1->idp = idp; //comput
  Nm1->local_V_size = local_V_size; //comput
  Nm1->flag_stgy = flag_stgy; //param
  Nm1->comm = comm; //param


  return 0;
}



int defineBlocks_avg(Block *tpltzblocks, double *T, int nb_blocks_loc, int n_block_avg, int lambda_block_avg, int64_t id0 )
{

int i;


  for ( i=0; i<nb_blocks_loc; i++)
    tpltzblocks[i].n = n_block_avg;

  for ( i=0; i<nb_blocks_loc; i++)
    tpltzblocks[i].lambda = lambda_block_avg;

  tpltzblocks[0].idv = (int64_t) (id0/n_block_avg) * n_block_avg ;
  for(i=1;i<nb_blocks_loc;i++)
    tpltzblocks[i].idv = (int64_t) tpltzblocks[i-1].idv + tpltzblocks[i-1].n;

  for( i=0; i<nb_blocks_loc; i++) {
    tpltzblocks[i].T_block = (T);
  }


  return 0;
}


//=============================================


int createRandomT(double *T, int Tsize)
{

  int i;
  srand (time (NULL));  //init seed

  //input matrix definition of T
    for(i=0;i<Tsize;i++)
      T[i]= rand()/((double) RAND_MAX);

  return 0;
}



int createTbasic1(double *T, int Tsize)
{

  int i;
  srand (Tsize);

  //input matrix definition of T
    for(i=0;i<Tsize;i++)
      T[i]= 1.0 + rand()/((double) RAND_MAX);

  return 0;
}



int createTbasic2(double *T, int Tsize)
{

  int i;
  srand (Tsize);

  //input matrix definition of T
    for(i=0;i<Tsize;i++) {
      if (i == 0) {
        T[i]=10.;}
      else if (i == 1) {
        T[i]=2.;}
      else if (i == 2) {
        T[i]=3.;}
      else {
        T[i]=rand()/((double) RAND_MAX);
     }}

  return 0;
}


int createTbasic3(double *T, int Tsize)
{

  int i;
  srand (Tsize);

  //input matrix definition of T
    for(i=0;i<Tsize;i++) {
      if (i == 0) {
        T[i]=2.;}
      else if (i == 1) {
        T[i]=-1.;}
      else if (i == 2) {
        T[i]=0.;}
      else {
        T[i]=0.;//rand()/((double) RAND_MAX);
     }}


  return 0;
}


int createTfrominvtt(double *T, int Tsize)
{

  int i;

//#include "invtt_params.h"

  double *invtt;

  T = invtt;
//  createinvtt(invtt);


  return 0;
}
