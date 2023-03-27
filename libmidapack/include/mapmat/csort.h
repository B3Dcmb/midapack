/** @file csort.h 
    @note Copyright (C) 2010 APC CNRS Université Paris Diderot @n This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. @n You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/gpl.html
    @author Pierre Cargemel
    @date March 2012*/

#ifndef MAPMAT_CSORT_H
#define MAPMAT_CSORT_H

#ifdef __cplusplus
extern "C" {
#endif

int ssort(int *indices, int count, int flag);

#if OPENMP
int omp_psort(int *indices, int count, int flag);
#endif

int sorted(int *indices, int count);

int monotony(int *indices, int count);

#ifdef __cplusplus
}
#endif

#endif //MAPMAT_CSORT_H
