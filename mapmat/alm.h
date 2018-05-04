/** @file alm.h
    @brief <b> Declaration of subroutines handling maps (associated sets of indices and values). </b>
    @note Copyright (C) 2010 APC CNRS Université Paris Diderot @n This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. @n You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/gpl.html
    @author Pierre Cargemel
    @date April 2012*/

/** @brief <b>Set some map values into a submap values array </b>
    @param mapval array of values
    @param submapval  array of values 
    @return array of indices*/ 
void m2s(double *mapval, double *submapval, int *subset, int count);

void s2m_sum(double *mapval, double *submapval, int *subset, int count);
void s2m(double *mapval, double *submapval, int *subset, int count);

void cnt_nnz_dot_prod(double *out, double *in, int cnt, int *ind, double *val, int nnz);
void lmatvecprod(int *ind, double *val, int m, int nnz, double *in, double *out);

#if OPENMP
void omp_cnt_nnz_dot_prod(double *out, double *in, int cnt, int *ind, double *val, int nnz);
void omp_lmatvecprod(int *ind, double *val, int m, int nnz, double *in, double *out);
#endif
int m2m(double *vA1, int *A1, int n1, double *vA2, int *A2, int n2);

int m2m_sum(double *vA1, int *A1, int n1, double *vA2, int *A2, int n2);
