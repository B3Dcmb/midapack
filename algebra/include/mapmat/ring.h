/** @file   ring.h
    @brief <b> Declaration of routines for ring-like communication scheme. </b>
    @note Copyright (C) 2010 APC CNRS Université Paris Diderot @n This program
   is free software; you can redistribute it and/or modify it under the terms of
   the GNU General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later version. This
   program is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
   A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
   @n You should have received a copy of the GNU General Public License along
   with this program; if not, see http://www.gnu.org/licenses/gpl.html
    @author Pierre Cargemel
    @date April 2012*/

#ifndef MAPMAT_RING_H
#define MAPMAT_RING_H

#ifdef W_MPI
#include <mpi.h>

int ring_init(int *indices, int count, int **R, int *nR, int **S, int *nS,
              int steps, MPI_Comm comm);

int ring_reduce(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax,
                double *val, double *res_val, int steps, MPI_Comm comm);

int ring_noempty_step_reduce(int **R, int *nR, int nRmax, int **S, int *nS,
                             int nSmax, double *val, double *res_val, int steps,
                             MPI_Comm comm);

int ring_nonblocking_reduce(int **R, int *nR, int **S, int *nS, double *val,
                            double *res_val, int steps, MPI_Comm comm);

int ring_noempty_reduce(int **R, int *nR, int nneR, int **S, int *nS, int nneS,
                        double *val, double *res_val, int steps, MPI_Comm comm);

int alltoallv_reduce(int **R, int *nR, int nRtot, int **S, int *nS, int nStot,
                     double *val, double *res_val, int steps, MPI_Comm comm);
#endif

#endif // MAPMAT_RING_H
