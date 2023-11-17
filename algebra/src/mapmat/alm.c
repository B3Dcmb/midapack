/** @file alm.c
    @brief Implementation of subroutines handling maps, distributions or
   functions. That means, almost all structures describes as sets of indices
   associated to sets of values).
    @note  Copyright (c) 2010-2012 APC CNRS Université Paris Diderot. This
   program is free software; you can redistribute it and/or modify it under the
   terms of the GNU Lesser General Public License as published by the Free
   Software Foundation; either version 3 of the License, or (at your option) any
   later version. This program is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
   General Public License for more details. You should have received a copy of
   the GNU General Public License along with this program; if not, see
   http://www.gnu.org/licenses/lgpl.html
    @note For more information about ANR MIDAS'09 project see
   http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
    @note ACKNOWLEDGMENT: This work has been supported in part by the French
   National  Research Agency (ANR) through COSINUS program (project MIDAS no.
   ANR-09-COSI-009).
    @author Pierre Cargemel
    @date April 2012*/

#include "alm.h"

/** Set some map values into a submap values array
    @param mapval array of values
    @param submapval array of values
    @return array of indices*/
void m2s(const double *mapval, double *submapval, const int *subset, int count,
         int mult) {
    for (int i = 0; i < count; i++) {
        for (int j = 0; j < mult; j++) {
            submapval[mult * i + j] = mapval[mult * subset[i] + j];
        }
    }
}

/** @brief assign submap values the submap values array
    @param mapval array of values
    @param submapval  array of values
    @return array of indices*/
void s2m(double *mapval, const double *submapval, const int *subset, int count,
         int mult) {
    for (int i = 0; i < count; i++) {
        for (int j = 0; j < mult; j++) {
            mapval[mult * subset[i] + j] = submapval[mult * i + j];
        }
    }
}

/** Function m2m for "map to map"
    Extract values from one map (A1, vA1), and for each pixel shared with an
   other map (A2, vA2), assign pixel value in vA1 and to pixel value in vA2.
    @return a number of elements shared between A1 and A2
    @sa m2m_sum
    @ingroup matmap_group22*/
__attribute__((unused)) int m2m(const double *vA1, const int *A1, int n1,
                                double *vA2, const int *A2, int n2) {
    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2) {
        if (A1[i] < A2[j]) {
            i++;
        } else if (A1[i] > A2[j]) {
            j++;
        } else {
            vA2[j] = vA1[i];
            k++;
            i++;
            j++;
        }
    }
    return k;
}

/** Function m2m_sum for "sum map to map"
    Extract values from one map (A1, vA1), and for each pixel shared with an
   other map (A2, vA2), sum pixel value in vA1 to pixel value in vA2.
    @return a number of elements shared between A1 and A2
    @sa m2m
    @ingroup matmap_group22*/
__attribute__((unused)) int m2m_sum(const double *vA1, const int *A1, int n1,
                                    double *vA2, const int *A2, int n2) {
    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2) {
        if (A1[i] < A2[j]) {
            i++;
        } else if (A1[i] > A2[j]) {
            j++;
        } else {
            vA2[j] += vA1[i];
            k++;
            i++;
            j++;
        }
    }
    return k;
}

/** Function m2m_sum_i for "sum map to map" (integer version)
    Extract values from one integer map (A1, vA1), and for each pixel shared
   with an other integer map (A2, vA2), sum pixel value in vA1 to pixel value in
   vA2.
    @return a number of elements shared between A1 and A2
    @sa m2m
    @ingroup matmap_group22*/
__attribute__((unused)) int m2m_sum_i(const int *vA1, const int *A1, int n1,
                                      int *vA2, const int *A2, int n2) {
    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2) {
        if (A1[i] < A2[j]) {
            i++;
        } else if (A1[i] > A2[j]) {
            j++;
        } else {
            vA2[j] += vA1[i];
            k++;
            i++;
            j++;
        }
    }
    return k;
}
