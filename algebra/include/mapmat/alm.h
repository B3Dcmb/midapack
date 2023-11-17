/** @file alm.h
    @brief <b> Declaration of subroutines handling maps (associated sets of
   indices and values). </b>
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

#ifndef MAPMAT_ALM_H
#define MAPMAT_ALM_H

void m2s(const double *mapval, double *submapval, const int *subset, int count,
         int mult);

void s2m(double *mapval, const double *submapval, const int *subset, int count,
         int mult);

__attribute__((unused)) int m2m(const double *vA1, const int *A1, int n1,
                                double *vA2, const int *A2, int n2);

__attribute__((unused)) int m2m_sum(const double *vA1, const int *A1, int n1,
                                    double *vA2, const int *A2, int n2);

__attribute__((unused)) int m2m_sum_i(const int *vA1, const int *A1, int n1,
                                      int *vA2, const int *A2, int n2);

#endif // MAPMAT_ALM_H
