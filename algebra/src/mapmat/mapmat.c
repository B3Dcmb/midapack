/** @file   mapmat.c
    @brief  Matrix routines implementation
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
    @date   November 2011*/

#include "mapmat.h"
#include "alm.h"
#include "cindex.h"
#include "csort.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Create a matrix specifying the number of local rows m,
    the number of non-zero elements per row nnz ,
    indices tab, values tab, ts_flags tab, flag for communication and a
   communicator comm. indices and values tabs must be allocated and contain at
   least m*nnz elements. It represents column indices of the nonzero elements.
   Respectively values tab represents the non-zero values. After call MatInit,
   all precomputation are done and the matrix structure is ready to use. That
   means you can start applying matrix operation. Another way to initialize a
   matrix structure is to apply step by step :
    - MatSetIndices
    - MatSetValues
    - MatLocalShape
    - MatComShape
    @warning do not modify indices tab until you will use the matrix structure.
    @warning [MPI COMM!] with Midapack sequential version, there is no
   communicator argument.
    @param A pointer to a Mat struct
    @param m number of local rows
    @param nnz number of non-zero per row
    @param indices input tab (modified)
    @param values input tab
    @param flag communication flag
    @param comm MPI communicator
    @ingroup matmap_group11
    @sa MatFree */
void MatInit(Mat *A, int m, int nnz, int *indices, double *values,
             uint8_t *flags
#ifdef W_MPI
             ,
             MPI_Comm comm
#endif
) {
    MatSetIndices(A, m, indices);

    MatSetValues(A, m, nnz, values);

    // set flags
    A->flags = flags;

    // compute lindices (local indexation)
    MatLocalShape(A, 3 /* counting sort */);
    A->trash_pix = 0;

#ifdef W_MPI
    // build communication scheme
    MatComShape(A, comm);
#endif
}

/** Set column indices of the nonzero elements.
    indices tab must be allocated and contains at least m elements.
    @param A pointer to a Mat struct
    @param m number of local rows
    @param indices input tab
    @return void
    @ingroup matmap_group11*/
void MatSetIndices(Mat *A, int m, int *indices) {
    A->m = m;             // set number of local rows
    A->indices = indices; // point to indices
}

/** Set values of the nonzero elements.
    values tab must be allocated and contains at least m*nnz values.
    @param A pointer to a Mat struct
    @param m number of local rows
    @param nnz number of non-zero per row
    @param values input tab
    @return void
    @ingroup matmap_group11*/
void MatSetValues(Mat *A, int m, int nnz, double *values) {
    A->m = m;           // set number of local rows
    A->nnz = nnz;       // set number of non-zero values per row
    A->values = values; // point to values
}

/** Free allocated tabs of a matrix structure including local indices tab and
   communication tabs. Do not free indices and values which user is responsible
   for.
    @param A pointer to a Mat struct
    @return void
    @sa MatInit MatLocalShape
    @ingroup matmap_group11 */
void MatFree(Mat *A) {
    free(A->lindices);
#if W_MPI
    free(A->com_indices);
#endif
}

/** Compute a local indices into a dense vector, lindices,
    and reindices indices tab according the this local dense vector.
    For this three steps are performed :
    - sort and merge indices tab,
    - allocate lindices of size lcount and copy the sorted indices
    - reindex indices according the local indices

    @warning lindices is internally allocated ( to free it, use MatFree )
    @sa MatComShape MatFree MatSetIndices
    @ingroup matmap_group11 */
void MatLocalShape(Mat *A, int sflag) {
    // allocate a tmp copy of indices tab to sort
    int *tmp_indices = malloc(sizeof(int) * A->m);
    if (tmp_indices == NULL) {
        fputs("allocation of tmp_indices in MatLocalShape failed", stderr);
        exit(EXIT_FAILURE);
    }

    memcpy(tmp_indices, A->indices, sizeof(int) * A->m);

    // sort tmp_indices
    // A->lcount = omp_psort(tmp_indices, A->m, sflag);
    A->lcount = ssort(tmp_indices, A->m, sflag);

    // allocate internal lindices tab
    A->lindices = malloc(sizeof(int) * A->lcount);
    if (A->lindices == NULL) {
        fputs("allocation of lindices in MatLocalShape failed", stderr);
        exit(EXIT_FAILURE);
    }

    // copy sorted tmp_indices and free
    memcpy(A->lindices, tmp_indices, sizeof(int) * A->lcount);
    free(tmp_indices);

    // re-index A->indices according to lindices
    sindex(A->lindices, A->lcount, A->indices, A->m);
}

#if W_MPI
/** Transform the matrix data structure, identifying columns shared by several
   processors
    @warning [MPI ONLY!] this function does not exist in Midapack sequential
   version
    @sa MatLocalShape MatInit TrMatVecProd
    @ingroup matmap_group11 */
void MatComShape(Mat *A, MPI_Comm comm) {
    // set communicator
    A->comm = comm;

    // determine maximum and minimum pixel indices
    int min, max;
    MPI_Allreduce(&(A->lindices[A->lcount - 1]), &max, 1, MPI_INT, MPI_MAX,
                  A->comm);
    // minimum index
    MPI_Allreduce(&(A->lindices[A->trash_pix]), &min, 1, MPI_INT, MPI_MIN,
                  A->comm);

    A->com_count = max - min + 1;
    A->com_indices = malloc(sizeof(int) * (A->lcount - A->trash_pix));
    if (A->com_indices == NULL) {
        fputs("failed to allocate com_indices in MatComShape", stderr);
        exit(EXIT_FAILURE);
    }

    int i = A->trash_pix;
    int j = 0;
    // same as subsetmap for a contiguous set
    while (j < A->com_count && i < A->lcount) {
        if (min + j < A->lindices[i]) {
            j++;
        } else {
            A->com_indices[i - A->trash_pix] = j;
            i++;
            j++;
        }
    }
}

void greedyreduce(Mat *A, double *x) {
    int size;
    MPI_Comm_size(A->comm, &size);

    // if there is only one process, return immediately
    if (size == 1) {
        return;
    }

    // allocate buffer that will be reduced
    double *lvalues = NULL;
    lvalues = malloc((sizeof *lvalues) * A->nnz * (A->lcount - A->trash_pix));
    if (lvalues == NULL) {
        int rank;
        MPI_Comm_rank(A->comm, &rank);
        fprintf(stderr, "[proc %d] malloc of lvalues failed in greedyreduce",
                rank);
        exit(EXIT_FAILURE);
    }

    // copy local values into result values
    memcpy(lvalues, x, (sizeof *x) * A->nnz * (A->lcount - A->trash_pix));

    double *com_val = calloc(A->nnz * A->com_count, sizeof(double));
    double *out_val = calloc(A->nnz * A->com_count, sizeof(double));
    s2m(com_val, lvalues, A->com_indices, A->lcount - A->trash_pix, A->nnz);
#if 0
    for (i = 0; i < A->com_count; i++) {
            printf("%lf ", com_val[i]);
        }
#endif
    MPI_Allreduce(com_val, out_val, A->com_count * A->nnz, MPI_DOUBLE, MPI_SUM,
                  A->comm);
#if 0
    for (i = 0; i < A->com_count; i++) {
            printf("%lf ", out_val[i]);
        }
#endif
    // sum receive buffer into values
    m2s(out_val, x, A->com_indices, A->lcount - A->trash_pix, A->nnz);
    free(com_val);
    free(out_val);
    free(lvalues);
}
#endif

/** Perform matrix-vector multiplication, \f$y \leftarrow A x\f$.
    @param A pointer to a Mat
    @param x input vector (overlapped)
    @param y output vector (distributed)
    @ingroup matmap_group11
    @ingroup matmap_group12a */
void MatVecProd(Mat *A, const double *x, double *y) {
    // refresh output vector
    for (int j = 0; j < A->m; j++)
        y[j] = 0.0;

    // elements not present in input vector (potentially zero)
    int n_ignore = A->ignore_extra ? A->trash_pix : 0;

    int nnz = A->nnz;
    for (int j = 0; j < A->m; j++) {
        if (!(A->ignore_extra) || A->indices[j] >= A->trash_pix) {
            for (int k = 0; k < nnz; k++) {
                int map_index = nnz * (A->indices[j] - n_ignore) + k;
                y[j] += A->values[nnz * j + k] * x[map_index];
            }
        }
    }
}

/** Perform a transposed matrix-vector multiplication, \f$x \leftarrow A^t y\f$
    using a precomputed communication scheme. Before calling this routine,
    the communication structure should have been set, calling MatInit or
   MatComShape. The routine can be divided in two steps :
    - a local matrix vector multiplication
    - a collective-reduce. it consits in a sum reduce over all processes.

    The collective reduce is performed using algorithm previously defined :
   ring, butterfly ...
    @sa MatVecProd MatComShape TrMatVecProd_Naive MatInit
    @param A a pointer to a Mat
    @param y local input vector (distributed)
    @param x local output vector (overlapped)
    @ingroup matmap_group11 */
void TrMatVecProd(Mat *A, const double *y, double *x) {
    // elements not present in output vector (potentially zero)
    int n_ignore = A->ignore_extra ? A->trash_pix : 0;

    int nnz = A->nnz;
    // refresh output vector
    for (int i = 0; i < A->lcount - n_ignore; i++) {
        for (int k = 0; k < nnz; k++) {
            x[nnz * i + k] = 0.0;
        }
    }

    for (int j = 0; j < A->m; j++) {
        if (!(A->ignore_extra) || A->indices[j] >= A->trash_pix) {
            for (int k = 0; k < A->nnz; k++) {
                int map_index = nnz * (A->indices[j] - n_ignore) + k;
                x[map_index] += A->values[nnz * j + k] * y[j];
            }
        }
    }

#ifdef W_MPI
    // perform global reduce on valid pixels
    greedyreduce(A, x + (A->ignore_extra ? 0 : nnz * A->trash_pix));
#endif
}
