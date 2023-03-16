/** @file als.h 
    @brief <b> Declaration of subroutines handling sets of indices or sets of values. </b>
    @note Copyright (C) 2010 APC CNRS Université Paris Diderot @n This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. @n You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/gpl.html
    @author Pierre Cargemel
    @date April 2012*/

#ifndef MAPMAT_ALS_H
#define MAPMAT_ALS_H

int card(int *A, int nA);

void merge(int *A, int nA, int *B);

int card_or(int *A1, int n1, int *A2, int n2);

int card_and(int *A1, int n1, int *A2, int n2);

int map_and(int *A1, int n1, int *A2, int n2, int *mapA1andA2);

int set_or(int *A1, int n1, int *A2, int n2, int *A1orA2);

int set_and(int *A1, int n1, int *A2, int n2, int *A1andA2);

void subset2map(int *A, int nA, int *subA, int nsubA);

#endif //MAPMAT_ALS_H
