/** @file   butterfly.c
    @brief  Implementation of routines for butterfly-like communication scheme.
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

#ifdef W_MPI
#include "mapmat.h"
#include <mpi.h>
#include <stdlib.h>
#include <string.h>


/** @brief Initialize tables for butterfly-like communication scheme
    This routine set up needed tables for the butterfly communication scheme.
    Sending and receiving tabs should be well allocated(at least size of number
   of steps in butterfly scheme). Double pointer are partially allocated, the
   last allocation is performed inside the routine. com_indices and com_count
   are also allocated inside the routine, thus they are passing by reference.
   They represent indices which have to be communicated an their number.
    Algotithm is based 2 parts.
    The first one identify intersection between processors indices, using 3
   successives butterfly communication schemes : bottom up, top down, and top
   down again. The second part works locally to build sets of indices to
   communicate
    @param indices set of indices(monotony) handle by a process.
    @param count number of elements
    @param R pointer to receiving maps
    @param nR array of number of elements in each receiving map
    @param S pointer to sending maps
    @param nS array of number of elements in each sending map
    @param com_indices set of indices(monotony) communicated by a process
    @param com_count number of elements
    @param steps number of communication exchange in the butterfly scheme
    @param comm MPI communicator
    @return 0 if no error
    @ingroup matmap_group22*/
int butterfly_init(int *indices, int count, int **R, int *nR, int **S, int *nS,
                   int **com_indices, int *com_count, int steps,
                   MPI_Comm comm) {

    int         i, k, p2k;
    int         rank, size, rk, sk;
    int         tag;
    MPI_Request s_request, r_request;
    int         nbuf, *buf;
    int       **I, *nI;
    int       **J, *nJ;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    I   = (int **) malloc(steps * sizeof(int *));
    nI  = (int *) malloc(steps * sizeof(int));
    tag = 0;
    p2k = size / 2;

    for (k = 0; k < steps;
         k++) { // butterfly first pass : bottom up (fill tabs nI and I)
        sk = (rank + size - p2k) % size;
        rk = (rank + p2k) % size;

        if (k == 0) { // S^0 := A
            nS[k] = count;
            S[k]  = (int *) malloc(nS[k] * sizeof(int));
            memcpy(S[k], indices, nS[k] * sizeof(int));
        } else { // S^k := S^{k-1} \cup R^{k-1}
            nS[k] = card_or(S[k - 1], nS[k - 1], I[steps - k], nI[steps - k]);
            S[k]  = (int *) malloc(nS[k] * sizeof(int));
            set_or(S[k - 1], nS[k - 1], I[steps - k], nI[steps - k], S[k]);
        }

        MPI_Irecv(&nI[steps - k - 1], 1, MPI_INT, rk, tag, comm,
                  &r_request); // receive number of indices
        MPI_Isend(&nS[k], 1, MPI_INT, sk, tag, comm,
                  &s_request); // send number of indices
        MPI_Wait(&r_request, MPI_STATUS_IGNORE);
        MPI_Wait(&s_request, MPI_STATUS_IGNORE);

        I[steps - k - 1] = (int *) malloc(nI[steps - k - 1] * sizeof(int));

        tag++;
        MPI_Irecv(I[steps - k - 1], nI[steps - k - 1], MPI_INT, rk, tag, comm,
                  &r_request); // receive indices
        MPI_Isend(S[k], nS[k], MPI_INT, sk, tag, comm,
                  &s_request); // send indices
        MPI_Wait(&r_request, MPI_STATUS_IGNORE);
        MPI_Wait(&s_request, MPI_STATUS_IGNORE);

        p2k /= 2;
        tag++;
    }

    J  = (int **) malloc(steps * sizeof(int *));
    nJ = (int *) malloc(steps * sizeof(int));

    tag = 0;
    p2k = 1;
    for (k = 0; k < steps;
         k++) { // buuterfly second pass : top down (fill tabs nJ and J)
        free(S[k]);
        sk = (rank + p2k) % size;
        rk = (rank + size - p2k) % size;
        if (k == 0) {
            nJ[k] = count;
            J[k]  = (int *) malloc(nJ[k] * sizeof(int));
            memcpy(J[k], indices, nJ[k] * sizeof(int));
        } else {
            nJ[k] = card_or(J[k - 1], nJ[k - 1], R[k - 1], nR[k - 1]);
            J[k]  = (int *) malloc(nJ[k] * sizeof(int));
            set_or(J[k - 1], nJ[k - 1], R[k - 1], nR[k - 1],
                   J[k]); // J^k=R^k-1 \cup J^k-1
            free(R[k - 1]);
        }
        if (k != steps - 1) {
            MPI_Irecv(&nR[k], 1, MPI_INT, rk, tag, comm, &r_request);
            MPI_Isend(&nJ[k], 1, MPI_INT, sk, tag, comm, &s_request);
            MPI_Wait(&r_request, MPI_STATUS_IGNORE);
            MPI_Wait(&s_request, MPI_STATUS_IGNORE);

            R[k] = (int *) malloc(nR[k] * sizeof(int));
            tag++;

            MPI_Irecv(R[k], nR[k], MPI_INT, rk, tag, comm, &r_request);
            MPI_Isend(J[k], nJ[k], MPI_INT, sk, tag, comm, &s_request);
            MPI_Wait(&r_request, MPI_STATUS_IGNORE);
            MPI_Wait(&s_request, MPI_STATUS_IGNORE);
        }
        p2k *= 2;
        tag++;
    }


    tag = 0;
    p2k = 1;
    for (k = 0; k < steps; k++) { // butterfly last pass : know that Sending tab
                                  // is S = I \cap J, so send S and we'll get R
        sk = (rank + p2k) % size;
        rk = (rank + size - p2k) % size;

        nS[k] = card_and(I[k], nI[k], J[k], nJ[k]);
        S[k]  = (int *) malloc(nJ[k] * sizeof(int));
        set_and(I[k], nI[k], J[k], nJ[k], S[k]); // S^k=I^k \cap J^k

        free(I[k]);
        free(J[k]);

        MPI_Irecv(&nR[k], 1, MPI_INT, rk, tag, comm,
                  &r_request); // receive size
        MPI_Isend(&nS[k], 1, MPI_INT, sk, tag, comm, &s_request); // send size
        MPI_Wait(&r_request, MPI_STATUS_IGNORE);
        MPI_Wait(&s_request, MPI_STATUS_IGNORE);

        R[k] = (int *) malloc(nR[k] * sizeof(int));
        tag++;

        MPI_Irecv(R[k], nR[k], MPI_INT, rk, tag, comm,
                  &r_request); // receive indices
        MPI_Isend(S[k], nS[k], MPI_INT, sk, tag, comm,
                  &s_request); // send indices
        MPI_Wait(&r_request, MPI_STATUS_IGNORE);
        MPI_Wait(&s_request, MPI_STATUS_IGNORE);

        p2k *= 2;
        tag++;
    }

    // Now we work locally
    int **USR, *nUSR, **U, *nU;

    USR  = (int **) malloc(steps * sizeof(int *));
    nUSR = (int *) malloc(steps * sizeof(int));
    U    = (int **) malloc(steps * sizeof(int *));
    nU   = (int *) malloc(steps * sizeof(int));

    for (k = 0; k < steps; k++) {
        nUSR[k] = card_or(S[k], nS[k], R[k], nR[k]);
        USR[k]  = (int *) malloc(nUSR[k] * sizeof(int));
        set_or(S[k], nS[k], R[k], nR[k], USR[k]);
    }
    for (k = 0; k < steps; k++) {
        if (k == 0) {
            nU[k] = nUSR[k];
            U[k]  = (int *) malloc(nU[k] * sizeof(int));
            memcpy(U[k], USR[k], nU[k] * sizeof(int));
        } else {
            nU[k] = card_or(U[k - 1], nU[k - 1], USR[k], nUSR[k]);
            U[k]  = (int *) malloc(nU[k] * sizeof(int *));
            set_or(U[k - 1], nU[k - 1], USR[k], nUSR[k], U[k]);
        }
    }
    *com_count   = nU[steps - 1];
    *com_indices = (int *) malloc(*com_count * sizeof(int));
    memcpy(*com_indices, U[steps - 1], *com_count * sizeof(int));
    //====================================================================

    for (k = 0; k < steps; k++) {
        subset2map(*com_indices, *com_count, S[k], nS[k]);
        subset2map(*com_indices, *com_count, R[k], nR[k]);
    }
    free(USR);
    free(U);

    return 0;
}


/** @brief Perform a sparse sum reduction (or mapped reduction) using a
   butterfly-like communication scheme
    @param R pointer to receiving maps
    @param nR array of number of elements in each receiving map
    @param nRmax maximum size of received message
    @param S pointer to sending maps
    @param nS array of number of elements in each sending map
    @param nSmax maximum size of sent message
    @param val set of values (typically values associated to communicated
   indices)
    @param steps number of communication exchange in the butterfly scheme
    @param comm MPI communicator
    @return 0 if no error
    @ingroup matmap_group22*/
int butterfly_reduce(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax,
                     double *val, int steps, MPI_Comm comm) {
    // double st, t;
    // t=0.0;
    int         k, p2k, tag;
    int         rank, size, rk, sk;
    MPI_Request s_request, r_request;
    double     *sbuf, *rbuf;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    sbuf = (double *) malloc(nSmax * sizeof(double));
    rbuf = (double *) malloc(nRmax * sizeof(double));
    tag  = 0;
    p2k  = 1;

    for (k = 0; k < steps; k++) {
        // st=MPI_Wtime();
        rk = (rank + size - p2k) % size;
        MPI_Irecv(rbuf, nR[k], MPI_DOUBLE, rk, tag, comm, &r_request);
        sk = (rank + p2k) % size;
        m2s(val, sbuf, S[k], nS[k]); // fill the sending buffer
        MPI_Isend(sbuf, nS[k], MPI_DOUBLE, sk, tag, comm, &s_request);
        MPI_Wait(&r_request, MPI_STATUS_IGNORE);
        s2m_sum(val, rbuf, R[k],
                nR[k]); // sum receive buffer into values //nR[k] floating sum
        p2k *= 2;
        tag++;
        MPI_Wait(&s_request, MPI_STATUS_IGNORE);
        // t=t+MPI_Wtime()-st;
    }
    free(sbuf);
    free(rbuf);
    return 0;
}

//===============================================Modification of the code by
//Sebastien Cayrols : 01/09/2015 ; Berkeley

/** @brief Perform a sparse sum reduction (or mapped reduction) using a
   butterfly-like communication scheme
    @param R pointer to receiving maps
    @param nR array of number of elements in each receiving map
    @param nRmax maximum size of received message
    @param S pointer to sending maps
    @param nS array of number of elements in each sending map
    @param nSmax maximum size of sent message
    @param val set of values (typically values associated to communicated
   indices)
    @param steps number of communication exchange in the butterfly scheme
    @param comm MPI communicator
    @return 0 if no error
    @ingroup matmap_group22*/
int butterfly_blocking_2instr_reduce(int **R, int *nR, int nRmax, int **S,
                                     int *nS, int nSmax, double *val, int steps,
                                     MPI_Comm comm) {
    // double st, t;
    // t=0.0;
    int        k, p2k, tag;
    int        rank, size, rk, sk;
    double    *sbuf, *rbuf;
    MPI_Status status;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    sbuf = (double *) malloc(nSmax * sizeof(double));
    rbuf = (double *) malloc(nRmax * sizeof(double));
    tag  = 0;
    p2k  = 1;

    for (k = 0; k < steps; k++) {
        // st=MPI_Wtime();
        sk = (rank + p2k) % size;
        m2s(val, sbuf, S[k], nS[k]); // fill the sending buffer
        MPI_Send(sbuf, nS[k], MPI_DOUBLE, sk, tag, comm);
        rk = (rank + size - p2k) % size;
        MPI_Recv(rbuf, nR[k], MPI_DOUBLE, rk, tag, comm, &status);
        s2m_sum(val, rbuf, R[k],
                nR[k]); // sum receive buffer into values //nR[k] floating sum
        p2k *= 2;
        tag++;
        // t=t+MPI_Wtime()-st;
    }
    free(sbuf);
    free(rbuf);
    return 0;
}

/** @brief Perform a sparse sum reduction (or mapped reduction) using a
   butterfly-like communication scheme
    @param R pointer to receiving maps
    @param nR array of number of elements in each receiving map
    @param nRmax maximum size of received message
    @param S pointer to sending maps
    @param nS array of number of elements in each sending map
    @param nSmax maximum size of sent message
    @param val set of values (typically values associated to communicated
   indices)
    @param steps number of communication exchange in the butterfly scheme
    @param comm MPI communicator
    @return 0 if no error
    @ingroup matmap_group22*/
int butterfly_blocking_1instr_reduce(int **R, int *nR, int nRmax, int **S,
                                     int *nS, int nSmax, double *val, int steps,
                                     MPI_Comm comm) {
    // double st, t;
    // t=0.0;
    int        k, p2k, tag;
    int        rank, size, rk, sk;
    double    *sbuf, *rbuf;
    MPI_Status status;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    sbuf = (double *) malloc(nSmax * sizeof(double));
    rbuf = (double *) malloc(nRmax * sizeof(double));
    tag  = 0;
    p2k  = 1;

    for (k = 0; k < steps; k++) {
        // st=MPI_Wtime();
        sk = (rank + p2k) % size;
        rk = (rank + size - p2k) % size;
        m2s(val, sbuf, S[k], nS[k]); // fill the sending buffer
        MPI_Sendrecv(sbuf, nS[k], MPI_DOUBLE, sk, tag, rbuf, nR[k], MPI_DOUBLE,
                     rk, tag, comm, &status);
        s2m_sum(val, rbuf, R[k],
                nR[k]); // sum receive buffer into values //nR[k] floating sum
        p2k *= 2;
        tag++;
        // t=t+MPI_Wtime()-st;
    }
    free(sbuf);
    free(rbuf);
    return 0;
}
#endif
