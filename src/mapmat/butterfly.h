/** @file   butterfly.h
    @brief <b> Declaration of routines for butterfly-like communication scheme. </b>
    @note Copyright (C) 2010 APC CNRS Université Paris Diderot @nThis program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. @nYou should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/gpl.html
    @author Pierre Cargemel
    @date April 2012*/


int butterfly_init(int *indices, int count, int **R, int *nR, int **S, int *nS, int **com_indices, int *com_count, int steps, MPI_Comm comm);

int butterfly_reduce(int **R, int *nR, int  nRmax, int **S, int *nS, int nSmax, double *val, int steps, MPI_Comm comm);

int butterfly_blocking_1instr_reduce(int **R, int *nR, int  nRmax, int **S, int *nS, int nSmax, double *val, int steps, MPI_Comm comm);
int butterfly_blocking_2instr_reduce(int **R, int *nR, int  nRmax, int **S, int *nS, int nSmax, double *val, int steps, MPI_Comm comm);

double butterfly_reduce_b(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax, double *val, int b, int steps, MPI_Comm comm);
