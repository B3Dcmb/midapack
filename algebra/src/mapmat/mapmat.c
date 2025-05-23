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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mapmat/alm.h>
#include <mapmat/bitop.h>
#include <mapmat/cindex.h>
#include <mapmat/csort.h>
#include <mapmat/mapmat.h>
#include <memutils.h>

#ifdef W_MPI
#include <mapmat/butterfly.h>
// #include <mapmat/butterfly_wrappers.h>
#include <mapmat/ring.h>
#include <mpi.h>
#endif

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
int MatInit(Mat *A, int m, int nnz, int *indices, double *values, int flag
#ifdef W_MPI
            ,
            MPI_Comm comm
#endif
) {
    int err;
    MatSetIndices(A, m, nnz, indices);

    MatSetValues(A, m, nnz, values);

    // compute lindices (local indexation)
    err = MatLocalShape(A, 3 /* counting sort */);

#ifdef W_MPI
    // build communication scheme
    err = MatComShape(A, flag, comm);
#endif
    return err;
}

/** Set column indices of the nonzero elements.
    indices tab must be allocated and contains at least m*nnz elements.
    @param A pointer to a Mat struct
    @param m number of local rows
    @param nnz number of non-zero per row
    @param indices input tab
    @return void
    @ingroup matmap_group11*/
void MatSetIndices(Mat *A, int m, int nnz, int *indices) {
    A->m = m;             // set number of local rows
    A->nnz = nnz;         // set number of non-zero values per row
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
    int err;
    A->m = m;           // set number of local rows
    A->nnz = nnz;       // set number of non-zero values per row
    A->values = values; // point to values
}

//===================Part added by Sebastien Cayrols to get amount of memory
// needed by communication algoritms
#if W_MPI
void CommInfo(Mat *A) {
    int i = 0, size, rank;
    double maxSizeR = 0.0;
    double maxSizeS = 0.0;
    double amountSizeR = 0.0;
    double amountSizeS = 0.0;
    double stepSum = 0.0, stepAvg = 0.0;
    // this value is based on data sent
    double minStep = 0.0, maxStep = 0.0;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    double *s = SAFEMALLOC(sizeof(double) * 4);
    double *r = SAFEMALLOC(sizeof(double) * 4 * 3);
    double *amountSizeByStep = SAFEMALLOC(sizeof(double) * A->steps);
    switch (A->flag) {
    case NONE:
        break;
    case BUTTERFLY:
        for (i = 0; i < A->steps; i++) {
            amountSizeR += A->nR[i];
            amountSizeS += A->nS[i];
            if (A->nR[i] > maxSizeR)
                maxSizeR = A->nR[i];
            if (A->nS[i] > maxSizeS)
                maxSizeS = A->nS[i];
        }
        break;
    //==========================Modification added by Sebastien Cayrols :
    // 01/09/2015 , Berkeley
    case BUTTERFLY_BLOCKING_1:
        for (i = 0; i < A->steps; i++) {
            amountSizeR += A->nR[i];
            amountSizeS += A->nS[i];
            if (A->nR[i] > maxSizeR)
                maxSizeR = A->nR[i];
            if (A->nS[i] > maxSizeS)
                maxSizeS = A->nS[i];
        }
        break;
    case BUTTERFLY_BLOCKING_2:
        for (i = 0; i < A->steps; i++) {
            amountSizeR += A->nR[i];
            amountSizeS += A->nS[i];
            if (A->nR[i] > maxSizeR)
                maxSizeR = A->nR[i];
            if (A->nS[i] > maxSizeS)
                maxSizeS = A->nS[i];
        }
        break;
    case NOEMPTYSTEPRING:
        for (i = 0; i < A->steps; i++) {
            amountSizeR += A->nR[i];
            amountSizeS += A->nS[i];
            if (A->nR[i] > maxSizeR)
                maxSizeR = A->nR[i];
            if (A->nS[i] > maxSizeS)
                maxSizeS = A->nS[i];
        }
        break;
    //==========================End modification
    case RING:
        for (i = 0; i < A->steps; i++) {
            amountSizeR += A->nR[i];
            amountSizeS += A->nS[i];
            if (A->nR[i] > maxSizeR)
                maxSizeR = A->nR[i];
            if (A->nS[i] > maxSizeS)
                maxSizeS = A->nS[i];
        }
        break;
    case NONBLOCKING:
        for (i = 0; i < A->steps; i++) {
            amountSizeR += A->nR[i];
            amountSizeS += A->nS[i];
            if (A->nR[i] > maxSizeR)
                maxSizeR = A->nR[i];
            if (A->nS[i] > maxSizeS)
                maxSizeS = A->nS[i];
        }
        break;
    case NOEMPTY:
        for (i = 0; i < A->steps; i++) {
            amountSizeR += A->nR[i];
            amountSizeS += A->nS[i];
            if (A->nR[i] > maxSizeR)
                maxSizeR = A->nR[i];
            if (A->nS[i] > maxSizeS)
                maxSizeS = A->nS[i];
        }
        break;
    case ALLTOALLV: // added -- rs 2015/02/04
        for (i = 0; i < A->steps; i++) {
            amountSizeR += A->nR[i];
            amountSizeS += A->nS[i];
        }
        break;
    case ALLREDUCE:
        amountSizeR = A->com_count;
        amountSizeS = A->com_count;
        maxSizeR = A->com_count;
        maxSizeS = A->com_count;
        break;
    }

    if (A->flag != ALLREDUCE && A->flag != ALLTOALLV) {
        double *t = SAFEMALLOC(sizeof(double) * A->steps);
        // Copy int array into double array
        for (i = 0; i < A->steps; i++)
            t[i] = A->nS[i];

        MPI_Reduce(t, amountSizeByStep, A->steps, MPI_DOUBLE, MPI_SUM, 0, comm);

        FREE(t);

        if (rank == 0) {
            stepSum = minStep = maxStep = amountSizeByStep[0];
            printf("\n[MEMORY]Step n°%4d, message size : %e", 0,
                   amountSizeByStep[0]);
            for (i = 1; i < A->steps; i++) {
                printf("\n[MEMORY]Step n°%4d, message size : %e", i,
                       amountSizeByStep[i]);
                if (minStep > amountSizeByStep[i])
                    minStep = amountSizeByStep[i];
                else if (maxStep < amountSizeByStep[i])
                    maxStep = amountSizeByStep[i];
                stepSum += amountSizeByStep[i];
            }
            stepAvg = stepSum / A->steps;
        }
    }
    s[0] = amountSizeR;
    s[1] = amountSizeS;
    s[2] = maxSizeR;
    s[3] = maxSizeS;
    MPI_Reduce(s, r, 4, MPI_DOUBLE, MPI_SUM, 0, comm);
    if (rank == 0)
        for (i = 0; i < 4; i++)
            r[i] /= size;
    MPI_Reduce(s, &r[4], 4, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(s, &r[8], 4, MPI_DOUBLE, MPI_MAX, 0, comm);
    if (rank == 0) {
        printf("\n[MEMORY]Step average            : %e\t[%e,%e]", stepAvg,
               minStep, maxStep);
        printf("\n[MEMORY]Amount of data received : %e\t[%e,%e]", r[0], r[4],
               r[8]);
        printf("\n[MEMORY]Amount of data sent     : %e\t[%e,%e]", r[1], r[5],
               r[9]);
        printf("\n[MEMORY]Message size received   : %e\t[%e,%e]", r[2], r[6],
               r[10]);
        printf("\n[MEMORY]Message size sent       : %e\t[%e,%e]\n", r[3], r[7],
               r[11]);
    }
    FREE(s);
    FREE(r);
    FREE(amountSizeByStep);
}
#endif

//===================End

/** Free allocated tabs of a matrix structure including local indices tab and
   communication tabs. Do not free indices and values which user is responsible
   for.
    @param A pointer to a Mat struct
    @return void
    @sa MatInit MatLocalShape
    @ingroup matmap_group11 */
void MatFree(Mat *A) {

    // get information about communication size
    // CommInfo(A);

    FREE(A->lindices);
#if W_MPI
    int rank;
    switch (A->flag) {
    case NONE:
        break;
    case BUTTERFLY:
#if 0
        MPI_Comm_rank(A->comm, &rank);
        free_butterfly_superstruct(A->bstruct, rank);
        break;
#endif
    //==========================Modification added by Sebastien Cayrols :
    // 01/09/2015 , Berkeley
    case BUTTERFLY_BLOCKING_1:
    case BUTTERFLY_BLOCKING_2:
        FREE(A->com_indices);
        FREE(A->R);
        FREE(A->nR);
        FREE(A->S);
        FREE(A->nS);
        break;
    case NOEMPTYSTEPRING:
    case RING:
    case NONBLOCKING:
    case NOEMPTY:
    case ALLTOALLV: // Added: rs 2015/02/04
        FREE(A->R);
        FREE(A->nR);
        FREE(A->S);
        FREE(A->nS);
        break;
    case ALLREDUCE:
        FREE(A->com_indices);
        break;
    }
#endif
}

/** Load matrix from a file.
    This is MatSave dual routine which loads data into matrix reading a
   specified file (or several specified files). Number of files should equal
   number of processor. File format should be ascii files (for examples look at
   files generated by MapMatSave routine).
    @warning Does not include gap samples flags functionality
    @todo Implement to read several file formats as basic ASCII, XML, HDF5...
    @param mat pointer to the Mat
    @param filename basename of a file, actually data are loaded from  several
   files denotes by "basename + processor number"
    @return error code
    @ingroup matmap_group11 */
int MatLoad(Mat *mat, char *filename) {
    int err;
    int rank;
#if W_MPI
    MPI_Comm_rank(mat->comm, &rank);
#else
    rank = 0;
#endif
    FILE *in;
    char fn[100];
    int i = 0;
    sprintf(fn, "%s_%d.dat", filename, rank);
    printf("%s", fn);
    in = fopen(fn, "r");
    if (in == NULL) {
        printf("cannot open file %s", fn);
        return 1;
    }
    while (feof(in) == 0 && i < (mat->m * mat->nnz)) {
        if (mat->nnz == 1) {
            fscanf(in, "%d %lf", &(mat->indices[i]), &(mat->values[i]));
        } else if (mat->nnz == 2) {
            fscanf(in, "%d %lf %d %lf", &(mat->indices[i]), &(mat->values[i]),
                   &(mat->indices[i + 1]), &(mat->values[i + 1]));
        } else {
            return 1; //(nnz > 2) not implement
        }
        i += mat->nnz;
    }
    if (i != mat->m * mat->nnz) {
        printf("WARNNING data size doesn't fit\n");
    }
    fclose(in);
    return 0;
}

/** Write matrix into files
    This is the dual routine of MatLoad.
    It saves matrix data into files, and can be usefull to check that data are
   well stored. Obviously it performs IO, moreover with one file per processor.
    Therfore, just call this function in a development phase, with few data and
   few processors.
    @warning Does not include gap samples flags functionality.
    @param A pointer to the Mat
    @param filename file basename, for instance passing "toto" should produce
   the output files named "toto_$(rank)"
    @return error code
    @ingroup matmap_group11 */
int MatSave(Mat *mat, char *filename) {
    FILE *out;
    char fn[100];
    int i, j;
    int rank;
#if W_MPI
    MPI_Comm_rank(mat->comm, &rank);
#else
    sprintf(fn, "%s_%d.dat", filename, rank);
#endif
    out = fopen(fn, "w");
    if (out == NULL) {
        printf("cannot open file %s", fn);
        return 1;
    }
    for (i = 0; i < (mat->nnz * mat->m); i += mat->nnz) {
        for (j = 0; j < mat->nnz; j++) {
            fprintf(out, "%d ", mat->indices[i + j]);
            fprintf(out, "%f ", mat->values[i + j]);
        }
        fprintf(out, "\n");
    }
    fclose(out);
    return 0;
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
int MatLocalShape(Mat *A, int sflag) {
    // allocate a tmp copy of indices tab to sort
    size_t count = A->m * A->nnz;
    int *tmp_indices = SAFEMALLOC(sizeof(int) * count);
    memcpy(tmp_indices, A->indices, sizeof(int) * count);

    // sort tmp_indices
    // A->lcount = omp_psort(tmp_indices, A->m * A->nnz, sflag);
    A->lcount = ssort(tmp_indices, count, sflag);
    A->lindices = SAFEMALLOC(sizeof(int) * A->lcount);
    memcpy(A->lindices, tmp_indices, A->lcount * sizeof(int));
    FREE(tmp_indices);

    sindex(A->lindices, A->lcount, A->indices, A->nnz * A->m);

    // count extra pixels
    int c = 0;
    while (A->lindices[c] < 0) {
        c += 1;
    }

    // store information in the pointing matrix
    A->trash_pix = c / A->nnz;

    return 0;
}

#if W_MPI

void _allocate_Mat_buffers(Mat *A) {
    // sending maps tab
    A->S = SAFEMALLOC(sizeof(int *) * A->steps);
    // receiving maps tab
    A->R = SAFEMALLOC(sizeof(int *) * A->steps);
    // sending map sizes tab
    A->nS = SAFEMALLOC(sizeof(int) * A->steps);
    // receiving map size tab
    A->nR = SAFEMALLOC(sizeof(int) * A->steps);
}

/** Transform the matrix data structure, identifying columns shared by several
   processors
    @warning [MPI ONLY!] this function does not exist in Midapack sequential
   version
    @sa MatLocalShape MatInit TrMatVecProd
    @ingroup matmap_group11 */
int MatComShape(Mat *A, int flag, MPI_Comm comm) {
    int rank, size;
    int i, min, max, j;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if ((flag == BUTTERFLY || flag == BUTTERFLY_BLOCKING_1 ||
         flag == BUTTERFLY_BLOCKING_2) &&
        is_pow_2(size) != 0) {
        if (rank == 0) {
            printf("BUTTERFLY_BLOCKING_1 or BUTTERFLY_BLOCKING_2 strategies "
                   "can only be used with 2^n processes -- switching to "
                   "BUTTERFLY.\n");
        }
        flag = ALLREDUCE; // FIXME: when generalized butterfly is implemented,
                          // switch back to BUTTERFLY
    }

    if (rank == 0) {
        printf("[MatComShape] communication strategy = %d\n", flag);
    }

    // set communicator and flag
    A->comm = comm;
    A->flag = flag;

    // prepare communication
    int ntrash = A->nnz * A->trash_pix;
    switch (A->flag) {
    case BUTTERFLY:
#if 0
        // Use Magdy's generalization to any number of processes
        A->bstruct = SAFEMALLOC(sizeof(struct Butterfly_superstruct));
        prepare_butterfly_communication(A->lindices + ntrash, // indices_in
                                        A->lcount - ntrash,   // count_in
                                        A->lindices + ntrash, // indices_out
                                        A->lcount - ntrash,   // count_out
                                        0,                    // classic
                                        A->bstruct, // initialized struct
                                        A->comm     // communicator
        );
        break;
#endif
    case BUTTERFLY_BLOCKING_1:
    case BUTTERFLY_BLOCKING_2:
        // butterfly-like schemes
        A->steps = log_2(size);

        // init butterfly-like communication
        _allocate_Mat_buffers(A);
        butterfly_init(A->lindices + ntrash, A->lcount - ntrash, A->R, A->nR,
                       A->S, A->nS, &(A->com_indices), &(A->com_count),
                       A->steps, A->comm);
        break;
    case NOEMPTYSTEPRING:
    case RING:
    case NONBLOCKING:
    case NOEMPTY:
    case ALLTOALLV:
        // ring-like schemes
        A->steps = size;

        // init ring-like communication
        _allocate_Mat_buffers(A);
        ring_init(A->lindices + ntrash, A->lcount - ntrash, A->R, A->nR, A->S,
                  A->nS, A->steps, A->comm);
        A->com_count = A->lcount - (A->nnz) * (A->trash_pix);
        A->com_indices = A->lindices + (A->nnz) * (A->trash_pix);
        break;
    case ALLREDUCE:
        // maximum index
        MPI_Allreduce(&(A->lindices[A->lcount - 1]), &max, 1, MPI_INT, MPI_MAX,
                      A->comm);
        MPI_Allreduce(&(A->lindices[ntrash]), &min, 1, MPI_INT, MPI_MIN,
                      A->comm);
        A->com_count = (max - min + 1);
        A->com_indices = SAFEMALLOC(sizeof(int) * (A->lcount - ntrash));
        i = ntrash;
        j = 0;
        while (j < A->com_count && i < A->lcount) {
            // same as subsetmap for a contiguous set
            if (min + j < A->lindices[i]) {
                j++;
            } else {
                A->com_indices[i - ntrash] = j;
                i++;
                j++;
            }
        }
        break;
    }
    return 0;
}
#endif

/** Perform matrix-vector multiplication, \f$y \leftarrow A x\f$.
    @param A pointer to a Mat
    @param x input vector (overlapped)
    @param y output vector (distributed)
    @ingroup matmap_group11
    @ingroup matmap_group12a */
void MatVecProd(const Mat *A, const double *x, double *y) {
    // refresh output vector
    for (int j = 0; j < A->m; j++)
        y[j] = 0.0;

    // elements not present in input vector (potentially zero)
    int extra = A->trash_pix * A->nnz;
    int off_extra = A->flag_ignore_extra ? extra : 0;

    for (int j = 0; j < A->m; j++) {
        int jnnz = j * A->nnz;
        if (!(A->flag_ignore_extra) || A->indices[jnnz] >= extra) {
            for (int k = 0; k < A->nnz; k++) {
                int map_index = A->indices[jnnz + k] - off_extra;
                y[j] += A->values[jnnz + k] * x[map_index];
            }
        }
    }
}

#ifdef W_MPI
/** Perform transposed matrix-vector multiplication, \f$x \leftarrow A^t y\f$.
    This naive version does not require a precomputed communication structure.
    But communication volumes may be significant.
    Consequently in most of the cases is not optimized.
    @warning [MPI ONLY!] this function does not exist in Midapack sequential
   version
    @sa TrMatVecProd MatLocalShape
    @param mat pointer
    @param y local input vector (distributed)
    @param x local output vector (overlapped)
    @ingroup matmap_group11
    @ingroup matmap_group12b */
int TrMatVecProd_Naive(Mat *A, double *y, double *x, int pflag) {
    int i, j, e, rank, size;
    int *rbuf, rbufcount;
    double *rbufvalues, *lvalues;
    int p, rp, sp, tag;
    MPI_Request s_request, r_request;
    MPI_Status status;

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);

    // allocate and set local values to 0.0
    lvalues = SAFECALLOC(A->lcount, sizeof(double));

    e = 0;
    for (i = 0; i < A->m; i++) {
        // local transform reduces
        for (j = 0; j < A->nnz; j++) {
            lvalues[A->indices[i * A->nnz + j]] +=
                (A->values[i * A->nnz + j]) * y[i];
        }
    }

    // copy local values into the result
    memcpy(x, lvalues, (A->lcount) * sizeof(double));
    // find the max communication buffer sizes
    MPI_Allreduce(&(A->lcount), &(rbufcount), 1, MPI_INT, MPI_MAX, A->comm);

    rbuf = SAFEMALLOC(sizeof(int) * rbufcount);
    rbufvalues = SAFEMALLOC(sizeof(double) * rbufcount);

    tag = 0;
    for (p = 1; p < size; p++) {
        // loop : collective global reduce in ring-like fashion
        rp = (size + rank - p) % size;
        sp = (rank + p) % size;
        // exchange sizes
        MPI_Send(&(A->lcount), 1, MPI_INT, sp, 0, A->comm);
        MPI_Recv(&rbufcount, 1, MPI_INT, rp, 0, A->comm, &status);
        tag++;
        // exchange local indices
        MPI_Irecv(rbuf, rbufcount, MPI_INT, rp, tag, A->comm, &r_request);
        MPI_Isend(A->lindices, A->lcount, MPI_INT, sp, tag, A->comm,
                  &s_request);
        MPI_Wait(&r_request, &status);
        MPI_Wait(&s_request, &status);
        tag++;
        // exchange local values
        MPI_Irecv(rbufvalues, rbufcount, MPI_DOUBLE, rp, tag, A->comm,
                  &r_request);
        MPI_Isend(lvalues, A->lcount, MPI_DOUBLE, sp, tag, A->comm, &s_request);
        tag++;
        MPI_Wait(&r_request, &status);
        // sum in the result
        m2m_sum(rbufvalues, rbuf, rbufcount, x, A->lindices, A->lcount);
        MPI_Wait(&s_request, &status);
    }
    FREE(lvalues);
    return 0;
}
#endif

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
void TrMatVecProd(const Mat *A, const double *y, double *x) {
    // elements not present in output vector (potentially zero)
    int extra = A->trash_pix * A->nnz;
    int off_extra = A->flag_ignore_extra ? extra : 0;

    // refresh output vector
    for (int i = 0; i < A->lcount - off_extra; i++) {
        x[i] = 0.0;
    }

    for (int j = 0; j < A->m; j++) {
        int jnnz = j * A->nnz;
        if (!(A->flag_ignore_extra) || A->indices[jnnz] >= extra) {
            for (int k = 0; k < A->nnz; k++) {
                int map_index = A->indices[jnnz + k] - off_extra;
                x[map_index] += A->values[jnnz + k] * y[j];
            }
        }
    }
#ifdef W_MPI
    // perform global reduce on valid pixels
    greedyreduce(A, x + (A->flag_ignore_extra ? 0 : extra));
#endif
}

#ifdef W_MPI
/** @brief Print information about a matrix.
    @n Usefull function to check, debug or bench. It prints matrix array sizes.
    @sa MatSave
    @param A pointer to the Mat
    @ingroup matmap_group11*/
int MatInfo(Mat *mat, int verbose, char *filename) {
    FILE *out;
    int *n;
    int *sr;
    int *s;
    int nnzline, sparsity, maxstep, maxsize, sumline, total;
    int i, j, k;
    char fn[100];
    int rank, size;
    int master = 0;
    MPI_Comm_rank(mat->comm, &rank);
    MPI_Comm_size(mat->comm, &size);

    if (rank == master) {
        // master process saves data into filename_info.txt
        sprintf(fn, "%s_%s", filename, "info.txt");
        out = fopen(fn, "w");
        if (out == NULL) {
            printf("cannot open file %s\n", fn);
            return 1;
        }
        printf("open file %s ...", fn);
        fprintf(out, "flag %d\n",
                mat->flag); // print matirx main description : flag
        // (communication scheme),
        fprintf(out, "rows %d\n ", mat->m); // rows per process,
        fprintf(out, "nnz %d\n", mat->nnz); // nnz (number of non zero per row).
        fprintf(out, "\n");                 // separator
    }

    /*n = (int* ) calloc(mat->lcount,sizeof(int));		//allocate
    //printf("before gather %d\n", rank);
    MPI_Gather(&(mat->lcount), 1, MPI_INT, n, 1, MPI_INT, master, mat->comm);
    //gather nbnonempty cols
    //printf("after gather %d\n", rank);

    if(rank==master){			//master process saves data into
    filename_info.txt fprintf(out, "cols :\n");	//nnz (number of non zero per
    row). for(i=0; i<size; i++)		// fprintf(out, "%d ", n[i]);
    //non-empty columns per process. fprintf(out, "\n"); //
    }
    free(n); */
    // free allocated tabs

    nnzline = 0; // compute communication sparsity and maximum message size
    sumline = 0;
    for (i = 0; i < mat->steps; i++) {
        //
        sumline += mat->nS[i];
        if (mat->nS[i] == 0) {
            //
            nnzline += 1; //
        } //
    } //
    MPI_Reduce(&nnzline, &sparsity, 1, MPI_INT, MPI_SUM, 0,
               mat->comm);                                           // sparsity
    MPI_Reduce(&sumline, &total, 1, MPI_INT, MPI_SUM, 0, mat->comm); // sparsity
    if (rank == master) {
        // master process saves data into filename_info.txt
        fprintf(out, "sparsity %d\n", sparsity); //
        fprintf(out, "total %d\n", total);       //
    }

    maxsize = 0;
    for (i = 0; i < mat->steps; i++) {
        //
        MPI_Reduce(&(mat->nS[i]), &maxstep, 1, MPI_INT, MPI_MAX, 0,
                   mat->comm); // maximum message size
        maxsize += maxstep;    //
    } //
    if (rank == master) {
        // master process saves data into filename_info.txt
        fprintf(out, "maxsize %d\n ", maxsize); //
        fprintf(out, "\n");                     // separator
    } //

    /* s = (int* ) calloc((mat->steps),sizeof(int));	//allocate steps
     MPI_Reduce(mat->nS, s, mat->steps, MPI_INT, MPI_SUM, 0, mat->comm);
     //imaximum message size

     if(rank==master){			//master process saves data into
     filename_info.txt
         fprintf(out, "sumsteps :\n");	//nnz (number of non zero per row).
         for(i=0; i<mat->steps; i++)		//
             fprintf(out, "%d ", s[i]);	//non-empty columns per process.
         fprintf(out, "\n");			//
     }
     free(s);

     if(verbose==1){
         sr = (int* ) calloc((mat->steps)*size,sizeof(int));	//allocate
     send/receive matrix
         //printf("before gather %d\n", rank);
         MPI_Gather(mat->nS, mat->steps, MPI_INT, sr, mat->steps, MPI_INT,
     master, mat->comm);
     //gather nbnonempty cols
         //printf("after gather %d\n", rank);

         if(rank==master){			//master process saves data into
     filename_info.txt fprintf(out, "send/receive matrix\n");	//separator
             for(i=0; i<size; i++){ 		//print collective description :
                 if(mat->flag==BUTTERFLY){		//send-receive matrix
                     for(j=0; j<size; j++){ 		//print send/receive
     matrix if(j>i){ if(is_pow_2(j-i)==0) fprintf(out,"%d ",
     sr[i*(mat->steps)+log_2(j-i)]); else fprintf(out,"%d ", 0);
                         }
                         else if(i>j){
                             if(is_pow_2(size+j-i)==0)
                                 fprintf(out,"%d ",
     sr[i*(mat->steps)+log_2(size+j-i)]); else fprintf(out,"%d ", 0);
                         }
                         else{
                             fprintf(out,"%d ", 0);
                         }
                     }
                     fprintf(out, "\n");
                 }
                 else{
                     for(j=0; j<size; j++){ 		//print send/receive
     matrix if(j>i){ fprintf(out,"%d ", sr[i*(mat->steps)+j-i]);
                         }
                         else if(i>j){
                             fprintf(out,"%d ", sr[(i+1)*(mat->steps)-i+j]);
                         }
                         else{
                             fprintf(out,"%d ", 0);
                         }
                     }
                     fprintf(out, "\n");
                 }
             }
         }
         free(sr);
     }*/

    if (rank == master) {
        // master process saves data into filename_info.txt
        fclose(out);
        printf("close %s\n", fn);
    }
    return 0;
}
#endif

int _array_max_or_zero(const int *arr, int size) {
    int max = 0;
    for (int i = 0; i < size; i++)
        if (arr[i] > max)
            max = arr[i];
    return max;
}

int _array_sum(const int *arr, int size) {
    int sum = 0;
    for (int i = 0; i < size; i++)
        sum += arr[i];
    return sum;
}

#if W_MPI
int greedyreduce(const Mat *A, double *x) {
    // don't communicate trash pixels
    const int n_trash = A->nnz * A->trash_pix;
    const int n_good = A->lcount - n_trash;
    int *lindices_good = A->lindices + n_trash;

    // copy local values into result values
    // /!\ assumes x contains only good pixels
    double *lvalues = SAFEMALLOC(sizeof *lvalues * n_good);
    memcpy(lvalues, x, sizeof *x * n_good);

    int nSmax, nRmax, nStot, nRtot, ne;

    double *com_val;
    double *out_val;

    switch (A->flag) {
    case BUTTERFLY:
        // max communication buffer size
        nRmax = _array_max_or_zero(A->nR, A->steps);
        nSmax = _array_max_or_zero(A->nS, A->steps);
        // allocate communication buffer
        com_val = SAFECALLOC(A->com_count, sizeof *com_val);
        m2m(lvalues, lindices_good, n_good, com_val, A->com_indices,
            A->com_count);
        butterfly_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val,
                         A->steps, A->comm);
        m2m(com_val, A->com_indices, A->com_count, x, lindices_good, n_good);
        // free communication buffer
        FREE(com_val);
        break;
    case BUTTERFLY_BLOCKING_1:
    case BUTTERFLY_BLOCKING_2:
        nRmax = _array_max_or_zero(A->nR, A->steps);
        nSmax = _array_max_or_zero(A->nS, A->steps);
        com_val = SAFECALLOC(A->com_count, sizeof *com_val);
        m2m(lvalues, lindices_good, n_good, com_val, A->com_indices,
            A->com_count);
        butterfly_blocking_1instr_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax,
                                         com_val, A->steps, A->comm);
        m2m(com_val, A->com_indices, A->com_count, x, lindices_good, n_good);
        FREE(com_val);
        break;
    case NOEMPTYSTEPRING:
    case RING:
        nRmax = _array_max_or_zero(A->nR, A->steps);
        nSmax = nRmax;
        ring_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, x,
                    A->steps, A->comm);
        break;
    case NONBLOCKING:
        ring_nonblocking_reduce(A->R, A->nR, A->S, A->nS, lvalues, x, A->steps,
                                A->comm);
        break;
    case NOEMPTY:
        ne = 0;
        for (int k = 1; k < A->steps; k++)
            if (A->nR[k] != 0)
                ne++;

        ring_noempty_reduce(A->R, A->nR, ne, A->S, A->nS, ne, lvalues, x,
                            A->steps, A->comm);
        break;
    case ALLREDUCE:
        com_val = SAFECALLOC(A->com_count, sizeof *com_val);
        out_val = SAFECALLOC(A->com_count, sizeof *out_val);
        s2m(com_val, lvalues, A->com_indices, n_good);
        MPI_Allreduce(com_val, out_val, A->com_count, MPI_DOUBLE, MPI_SUM,
                      A->comm);
        m2s(out_val, x, A->com_indices, n_good);
        FREE(com_val);
        FREE(out_val);
        break;
    case ALLTOALLV:
        // compute buffer sizes
        nRtot = _array_sum(A->nR, A->steps);
        nStot = _array_sum(A->nS, A->steps);
        alltoallv_reduce(A->R, A->nR, nRtot, A->S, A->nS, nStot, lvalues, x,
                         A->steps, A->comm);
        break;
    default:
        // do nothing
        break;
    }
    FREE(lvalues);
    return 0;
}
#endif
