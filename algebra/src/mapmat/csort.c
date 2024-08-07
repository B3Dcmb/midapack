/** @file csort.c
    @brief subroutines for sequential or parallel sort and/or merge sets of
   integer.
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
    @date April 2012 */

#include <stdlib.h>

#include <mapmat/csort.h>

/** Insertion sort :
    complexity is n square; ascending order
    @param indices array of integer
    @param count  number of elements in indices array
    @return void*/
void insertion_sort(int *indices, int count) {
    int i, j;
    int tmp;
    for (i = 0; i < count - 1; i++) {
        tmp = indices[i + 1];
        j = i;
        while (j != -1 && tmp < indices[j]) {
            indices[j + 1] = indices[j];
            indices[j] = tmp;
            j--;
        }
    }
}

/** Quick sort :
    complexity n square (Average is n log(n)).
    Sort in ascending order indices between left index and right index.
    Notice that this algorithm uses recursive calls, and can lead to stack
   overflow.
    @param indices array of integer
    @param first element number
    @param last element number
    @return void*/
void quick_sort(int *indices, int left, int right) {
    int pivot;
    int tmp, key;
    int i, j;
    if (left < right) {
        key = indices[left];
        i = left + 1;
        j = right;
        while (i <= j) {
            while ((i <= right) && (indices[i] <= key))
                i++;
            while ((j > left) && (indices[j] > key))
                j--;
            if (i < j) {
                tmp = indices[i];
                indices[i] = indices[j];
                indices[j] = tmp;
                i++;
                j--;
            }
        }
        tmp = indices[left];
        indices[left] = indices[j];
        indices[j] = tmp;
        pivot = j;
        quick_sort(indices, left, pivot - 1);
        quick_sort(indices, pivot + 1, right);
    }
}

/** Bubble sort :
    complexity n square
    @param indices array of integer
    @param count  number of elements in indices array
    @return void **/
void bubble_sort(int *indices, int count) {
    int i, j, tmp;
    for (i = (count - 1); i > 0; i--) {
        for (j = 1; j <= i; j++) {
            if (indices[j - 1] > indices[j]) {
                tmp = indices[j - 1];
                indices[j - 1] = indices[j];
                indices[j] = tmp;
            }
        }
    }
}

/** Counting sort :
    sort indices in strictly ascending order (monotony).
    If the array initially ows reundants elements, they will be merged.
    Thus the returned number of sorted elements is less than the initial number
   of elements. Assuming elements belong to [min, max], complexity is
   n+(max-min+1). Due to allocation, memory overhead is
   (max-min+1)sizeof(element)
    @param indices array of integer
    @param count  number of elements in indices array
    @return number of sorted elements*/
int counting_sort(int *indices, int count) {
    int *buf;
    int i, j, k;
    int min, max;
    min = indices[0];
    max = indices[0];
    for (i = 1; i < count; i++) {
        if (indices[i] > max) {
            max = indices[i];
        } else {
            if (indices[i] < min)
                min = indices[i];
        }
    }
    buf = (int *)calloc(max - min + 1, sizeof(int));
    for (i = 0; i < count; i++) {
        buf[indices[i] - min] = 1;
    }
    j = 0;
    for (i = 0; i < (max - min + 1); i++) {
        if (buf[i] == 1) {
            indices[j] = min + i;
            j++;
        }
    }
    free(buf);
    return j;
}

/** Shell sort
    @param indices array of integer
    @param count  number of elements in indices array
    @return void*/
void shell_sort(int *a, int n) {
    int j, i, k, m, mid;
    for (m = n / 2; m > 0; m /= 2) {
        for (j = m; j < n; j++) {
            for (i = j - m; i >= 0; i -= m) {
                if (a[i + m] >= a[i])
                    break;
                else {
                    mid = a[i];
                    a[i] = a[i + m];
                    a[i + m] = mid;
                }
            }
        }
    }
}

/** Sort and merge redundant elements of a set of indices using a specified
   method. The indices tab, initially an arbitrary set of integers, becomes a
   monotony set. Available methods :
    - quick sort
    - bubble sort
    - insertion sort
    - counting sort
    - shell sort

    @param indices tab (modified)
    @param count number of indices
    @param flag method
    @return number of sorted elements
    @ingroup matmap_group22*/
int ssort(int *indices, int count, int flag) {
    int i, n;
    int *ptr_i, *ptr_o;
    switch (flag) {
    case 0:
        quick_sort(indices, 0, count - 1);
        break;
    case 1:
        bubble_sort(indices, count);
        break;
    case 2:
        insertion_sort(indices, count);
        break;
    case 3:
        n = counting_sort(indices, count);
        return n;
    case 4:
        shell_sort(indices, count);
        break;
    }
    ptr_i = indices;
    ptr_o = indices;
    n = 1;
    for (i = 0; i < count - 1; i++) {
        ptr_i++;
        if (*ptr_i != *ptr_o) {
            ptr_o++;
            n++;
            *ptr_o = *ptr_i;
        }
    }
    return n;
}

// optimized version is faster than the other implementation but there is a bug
// !!!
#ifdef W_OPENMP

#include <omp.h>
#define PAGE 1024

int omp_psort_opt(int *A, int nA, int flag) {
    int i;
    int *count, *disp;
    int q, r;
    int p, l;
    int tid, nths;
    int *buf;
    int *ptr_i, *ptr_o;
    int n, k, d;

#pragma omp parallel shared(nths)
    { //---fork---just to get the number of threads
        nths = omp_get_num_threads();
    } //---join---

    p = nA / PAGE; // number of full pages
    q = p / nths;  // full pages per thread
    l = p % nths;  // full pages left
    r = nA % PAGE; // number of elements the last page not full

    count = (int *)malloc(nths * sizeof(int));
    disp = (int *)malloc(nths * sizeof(int));

    for (i = 0; i < nths; i++) {
        count[i] = q * PAGE;
        if (i < l)
            count[i] += PAGE;
        if (i == l)
            count[i] += r;
    }

    disp[0] = 0;
    for (i = 0; i < nths - 1; i++) {
        disp[i + 1] = disp[i] + count[i];
    }

#pragma omp parallel private(tid, n, k, d, buf) shared(nths, A, nA, disp, count)
    { //---fork---1st step, sort on local chunk
        tid = omp_get_thread_num();

        buf = (int *)malloc(nA * sizeof(int));
        // buf = (int *) malloc(count[tid]*sizeof(int));
        memcpy(buf, A + disp[tid], count[tid] * sizeof(int));

        n = ssort(buf, count[tid], flag);
        count[tid] = n;

        memcpy(A + disp[tid], buf, n * sizeof(int));

        nths = omp_get_num_threads();

        k = nths;
        d = 1;
        while (k > 1) {
#pragma omp barrier
            if (tid % (2 * d) == 0 && tid + d < nths) {
                set_or(A + disp[tid], count[tid], A + disp[tid + d],
                       count[tid + d], buf);
                count[tid] = n;
                memcpy(A + disp[tid], buf, n * sizeof(int));
                d *= 2;
                k /= 2;
            }
        }
        free(buf);
    } //---join---

    nA = count[0];
    //  printf("\nA :");
    //  for(i=0; i<nA; i++){
    //    printf(" %d",A[i]);
    //  }
    free(count);
    free(disp);
    return nA;
}

/** Sort and merge redundant elements of a set of indices, using openMP.
    The indices tab, initially an arbitrary set of integers, becomes a monotony
   set. Algorithm is diivided in two steps :
    - each thread sorts, in parallel, a subpart of the set using a specified
   method.
    - subsets obtained are merged successively in a binary tree manner.

    Available methods for the fully parallel step :
    - quick sort
    - bubble sort
    - insertion sort
    - counting sort
    - shell sort

    @param indices tab (modified)
    @param count number of elements to sort
    @return flag method
    @ingroup matmap_group22*/
int omp_psort(int *A, int nA, int flag) {
    int i;
    int *count, *disp;
    int q, r;
    int tid, nths;
    int *buf;
    int *ptr_i, *ptr_o;
    int n, k, d;

#pragma omp parallel private(tid) shared(nths)
    { //---fork---just to get the number of threads
        nths = omp_get_num_threads();
    } //---join---

    q = nA / nths;
    r = nA % nths;

    count = (int *)malloc(nths * sizeof(int));
    disp = (int *)malloc(nths * sizeof(int));

    for (i = 0; i < nths; i++) {
        if (i < r) {
            count[i] = q + 1;
        } else {
            count[i] = q;
        }
    }

    disp[0] = 0;
    for (i = 0; i < nths - 1; i++) {
        disp[i + 1] = disp[i] + count[i];
    }

#pragma omp parallel private(tid, n) shared(A, disp, count)
    { //---fork---1st step, sort on local chunk
        tid = omp_get_thread_num();
        n = ssort(A + disp[tid], count[tid], flag);
        count[tid] = n;
    } //---join---

    buf = (int *)malloc(nA * sizeof(int));

#pragma omp parallel private(tid, n, k, d) shared(nths, nA, A, disp, count, buf)
    { //---fork---2nd step, gathering with a binary tree scheme
        tid = omp_get_thread_num();
        nths = omp_get_num_threads();
        k = nths;
        d = 1;
        while (k > 1) {
            if (tid % (2 * d) == 0 && tid + d < nths) {
                //        printf("\nd %d, thread %d, count+ %d, disp %d",d ,
                //        tid, count[tid], disp[tid]);
                n = card_or(A + disp[tid], count[tid], A + disp[tid + d],
                            count[tid + d]);
                set_or(A + disp[tid], count[tid], A + disp[tid + d],
                       count[tid + d], buf + disp[tid]);
                count[tid] = n;
                memcpy(A + disp[tid], buf + disp[tid], n * sizeof(int));
                d *= 2;
                k /= 2;
            }
#pragma omp barrier
        }
    } //---join---

    nA = count[0];
    free(buf);
    free(count);
    free(disp);
    return nA;
}
#endif

///=========================Checker
/// Routines=============================================
int sorted(int *indices, int count) {
    int i = 0;
    while (i < count - 2) {
        if (indices[i] > indices[i + 1]) {
            return 1;
        } else {
            i++;
        }
    }
    return 0;
}

int monotony(int *indices, int count) {
    int i = 0;
    while (i < count - 2) {
        if (indices[i] >= indices[i + 1]) {
            return 1;
        } else {
            i++;
        }
    }
    return 0;
}
