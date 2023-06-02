/** @file   ring.c
    @brief Implementation of routines for ring-like communication scheme.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** @brief Initialize tables for ring-like communication scheme.

    @n This routine set up needed tables for the ring communication scheme.
    Sending and receiving tabs should be well allocated(at least size of number
   of steps in ring scheme). Double pointer are partially allocated, the last
   allocation is performed inside the routine (only for R S are just pointer).
    @param indices set of indices(monotony) handle by a process.
    @param count number of elements
    @param R pointer to receiving maps
    @param nR array of number of elements in each receiving map
    @param S pointer to sending maps
    @param nS array of number of elements in each sending map
    @param com_indices set of indices(monotony) communicated by a process
    @param com_count number of elements
    @param steps number of communication exchange in the ring scheme
    @param comm MPI communicator
    @todo Ring loop and ring table are set from index 1 to size. Should be shift
   and be set from index 0 to size-1.
    @return 0 if no error
    @ingroup matmap_group22*/
int ring_init(int *indices, int count, int **R, int *nR, int **S, int *nS,
              int steps, MPI_Comm comm) {
    int         err, p, tag;
    int         size, rank, sp, rp;
    int        *buf, nbuf;
    MPI_Request s_request, r_request;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    MPI_Allreduce(&count, &nbuf, 1, MPI_INT, MPI_MAX,
                  comm); // compute the buffer size : max(count)_{comm}
    buf = (int *) malloc(nbuf * sizeof(int)); // allocate buffer
    tag = 0;
    for (p = 1; p < steps; p++) { // communication phase to get nb shared
                                  // indices between peer of preocesses
        sp = (rank + p) % size;
        rp = (rank + size - p) % size;
        MPI_Isend(&count, 1, MPI_INT, sp, 0, comm,
                  &s_request); // send my number of indices
        MPI_Irecv(&nbuf, 1, MPI_INT, rp, 0, comm,
                  &r_request); // receive a number of indices
        tag++;
        MPI_Wait(&r_request, MPI_STATUS_IGNORE);
        MPI_Irecv(buf, nbuf, MPI_INT, rp, tag, comm,
                  &r_request); // receive indices tab
        MPI_Wait(&s_request, MPI_STATUS_IGNORE);
        MPI_Isend(indices, count, MPI_INT, sp, tag, comm,
                  &s_request); // send indices tab
        tag++;

        MPI_Wait(&r_request, MPI_STATUS_IGNORE);
        nR[p]         = card_and(indices, count, buf,
                                 nbuf); // compute number of shared indices
        nS[steps - p] = nR[p];
        R[p] = (int *) malloc(nR[p] * sizeof(int));   // allocate receiving tab
        S[steps - p] = (int *) malloc(nS[steps - p]
                                      * sizeof(int)); // allocate sanding tab
        map_and(indices, count, buf, nbuf, R[p]);     // fill receiving tab
        S[steps - p] = R[p];                          //
    }
    free(buf);
    nS[0] = 0; //
    nR[0] = 0; //
    return 0;
}

/** @brief Perform a sparse sum reduction (or mapped reduction) using a
   ring-like communication scheme

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
int ring_reduce(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax,
                double *val, double *res_val, int steps, MPI_Comm comm) {
    int         tag, rank, size, p;
    MPI_Request s_request, r_request;
    int         sp, rp;
    double     *sbuf, *rbuf;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    tag = 0;

    rbuf = (double *) malloc(nRmax * sizeof(double));
    sbuf = (double *) malloc(nSmax * sizeof(double));

    for (p = 1; p < steps; p++) {
        rp = (rank + size - p) % size;
        MPI_Irecv(rbuf, nR[p], MPI_DOUBLE, rp, tag, comm, &r_request);
        sp = (rank + p) % size;
        m2s(val, sbuf, S[p], nS[p]); // fill the sending buffer
        MPI_Isend(sbuf, nS[p], MPI_DOUBLE, sp, tag, comm, &s_request);

        tag++;

        MPI_Wait(&r_request, MPI_STATUS_IGNORE);
        s2m_sum(res_val, rbuf, R[p], nR[p]); // sum receive buffer into values

        MPI_Wait(&s_request, MPI_STATUS_IGNORE);
    }
    free(sbuf);
    free(rbuf);
    return 0;
}


/** @brief Perform a sparse sum reduction (or mapped reduction) using an
 MPI-Alltoallv call

 @param R pointer to receiving maps
 @param nR array of number of elements in each receiving map
 @param nRtot size of the receive buffer
 @param S pointer to sending maps
 @param nS array of number of elements in each sending map
 @param nStot size of the send buffer
 @param val set of values (typically values associated to communicated indices)
 @param steps number of communication exchange in the butterfly scheme
 @param comm MPI communicator
 @return 0 if no error
 @ingroup matmap_group22*/
int alltoallv_reduce(int **R, int *nR, int nRtot, int **S, int *nS, int nStot,
                     double *val, double *res_val, int steps, MPI_Comm comm) {
    int         rank, size, p;
    MPI_Request s_request, r_request;
    int         sp, rp, *rindx, *sindx, *rdisp, *sdisp;
    double     *sbuf, *rbuf;


    MPI_Comm_size(comm,
                  &size); // N.B. size and steps must be equal, shall we check
                          // for this ?! -- rs
    MPI_Comm_rank(comm, &rank);

    rbuf = (double *) malloc(nRtot * sizeof(double));
    sbuf = (double *) malloc(nStot * sizeof(double));

    rindx = (int *) calloc(size, sizeof(int));
    sindx = (int *) calloc(size, sizeof(int));

    rdisp = (int *) calloc(size, sizeof(int));
    sdisp = (int *) calloc(size, sizeof(int));

    // compute shifts ...

    for (p = 0; p < steps; p++) { // starts with 0 !
        rp        = (rank + size - p) % size;
        rindx[rp] = nR[p];
        sp        = (rank + p) % size;
        sindx[sp] = nS[p];
    }

    for (p = 1; p < size; p++) {
        sdisp[p] = sdisp[p - 1] + sindx[p - 1];
        rdisp[p] = rdisp[p - 1] + rindx[p - 1];
    }

    // prepare data to send ...

    for (p = 0; p < steps; p++) {
        sp = (rank + p) % size;
        m2s(val, &sbuf[sdisp[sp]], S[p], nS[p]); // fill the sending buffer
    }

    MPI_Alltoallv(sbuf, sindx, sdisp, MPI_DOUBLE, rbuf, rindx, rdisp,
                  MPI_DOUBLE, comm);

    // accumulate contributions ...

    for (p = 0; p < steps; p++) {
        rp = (rank + size - p) % size;
        s2m_sum(res_val, &rbuf[rdisp[rp]], R[p],
                nR[p]); // sum receive buffer into values
    }

    free(sdisp);
    free(rdisp);
    free(sindx);
    free(rindx);
    free(sbuf);
    free(rbuf);

    return 0;
}

/** @brief Perform a sparse sum reduction (or mapped reduction) using a
   ring-like non-blocking communication scheme

    @param R pointer to receiving maps
    @param nR array of number of elements in each receiving map
    @param S pointer to sending maps
    @param nS array of number of elements in each sending map
    @param val set of values (typically values associated to communicated
   indices)
    @param steps number of communication exchange in the butterfly scheme
    @param comm MPI communicator
    @return 0 if no error
    @ingroup matmap_group22*/
int ring_nonblocking_reduce(int **R, int *nR, int **S, int *nS, double *val,
                            double *res_val, int steps, MPI_Comm comm) {
    int          tag, rank, size, p;
    MPI_Request *s_request, *r_request;
    int          sp, rp;
    double     **sbuf, **rbuf;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    // printf("\n non_blocking rank %d", rank);

    s_request = (MPI_Request *) malloc((steps - 1) * sizeof(MPI_Request));
    r_request = (MPI_Request *) malloc((steps - 1) * sizeof(MPI_Request));

    rbuf = (double **) malloc((steps - 1) * sizeof(double *));
    sbuf = (double **) malloc((steps - 1) * sizeof(double *));

    for (p = 1; p < steps; p++) {
        // printf("\n buf alloc %d", p);
        rbuf[p - 1] = (double *) malloc(nR[p] * sizeof(double));
        sbuf[p - 1] = (double *) malloc(nS[p] * sizeof(double));
        m2s(val, sbuf[p - 1], S[p], nS[p]); // fill the sending buffer
    }

    tag = 0;
    for (p = 1; p < steps; p++) {
        // printf("\n isend  %d", p);
        sp = (rank + p) % size;
        rp = (rank + size - p) % size;

        MPI_Irecv(rbuf[p - 1], nR[p], MPI_DOUBLE, rp, tag, comm,
                  &r_request[p - 1]);
        MPI_Isend(sbuf[p - 1], nS[p], MPI_DOUBLE, sp, tag, comm,
                  &s_request[p - 1]);
        tag++;
    }
    MPI_Waitall(size - 1, r_request, MPI_STATUSES_IGNORE);

    for (p = 1; p < steps; p++) {
        s2m_sum(res_val, rbuf[p - 1], R[p],
                nR[p]); // sum receive buffer into values
    }
    MPI_Waitall(size - 1, s_request, MPI_STATUSES_IGNORE);
    free(r_request);
    free(s_request);
    free(sbuf);
    free(rbuf);
    return 0;
}

/** @brief Perform a sparse sum reduction (or mapped reduction) using a
   ring-like non-blocking no-empty communication scheme
    @param R pointer to receiving maps
    @param nR array of number of elements in each receiving map
    @param nneR number of no-empty receiving messages
    @param S pointer to sending maps
    @param nS array of number of elements in each sending map
    @param nneS number of no-empty sending messages
    @param val set of values (typically values associated to communicated
   indices)
    @param steps number of communication exchange in the butterfly scheme
    @param comm MPI communicator
    @return 0 if no error
    @ingroup matmap_group22*/
int ring_noempty_reduce(int **R, int *nR, int nneR, int **S, int *nS, int nneS,
                        double *val, double *res_val, int steps,
                        MPI_Comm comm) {
    int          tag, rank, size, p;
    MPI_Request *s_request, *r_request;
    int          sp, rp, nesi, neri;
    double     **sbuf, **rbuf;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    // printf("\n non_blocking rank %d", rank);

    s_request = (MPI_Request *) malloc(nneS * sizeof(MPI_Request));
    r_request = (MPI_Request *) malloc(nneR * sizeof(MPI_Request));

    rbuf = (double **) malloc(nneR * sizeof(double *));
    sbuf = (double **) malloc(nneS * sizeof(double *));

    nesi = 0;
    for (p = 1; p < steps; p++) {
        if (nS[p] != 0) {
            sbuf[nesi] = (double *) malloc(nS[p] * sizeof(double));
            m2s(val, sbuf[nesi], S[p], nS[p]); // fill the sending buffer
            nesi++;
        }
    }

    tag  = 0;
    nesi = 0;
    neri = 0;
    for (p = 1; p < steps; p++) {
        sp = (rank + p) % size;
        rp = (rank + size - p) % size;
        if (nR[p] != 0) {
            rbuf[neri] = (double *) malloc(nR[p] * sizeof(double));
            MPI_Irecv(rbuf[neri], nR[p], MPI_DOUBLE, rp, tag, comm,
                      &r_request[neri]);
            neri++;
        }
        if (nS[p] != 0) {
            MPI_Isend(sbuf[nesi], nS[p], MPI_DOUBLE, sp, tag, comm,
                      &s_request[nesi]);
            nesi++;
        }
        tag++;
    }
    MPI_Waitall(nneR, r_request, MPI_STATUSES_IGNORE);

    neri = 0;
    for (p = 1; p < steps; p++) {
        if (nR[p] != 0) {
            s2m_sum(res_val, rbuf[neri], R[p],
                    nR[p]); // sum receive buffer into values
            neri++;
        }
    }
    MPI_Waitall(nneS, s_request, MPI_STATUSES_IGNORE);
    free(r_request);
    free(s_request);
    free(sbuf);
    free(rbuf);
    return 0;
}

//=======================================================Modification added by
//Sebastien Cayrols : 01/09/2015 ; Berkeley

/** @brief Perform a sparse sum reduction (or mapped reduction) using a
   ring-like communication scheme

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
int ring_noempty_step_reduce(int **R, int *nR, int nRmax, int **S, int *nS,
                             int nSmax, double *val, double *res_val, int steps,
                             MPI_Comm comm) {
    int         tag, rank, size, p;
    MPI_Request s_request, r_request;
    int         sp, rp;
    double     *sbuf, *rbuf;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    tag = 0;

    rbuf = (double *) malloc(nRmax * sizeof(double));
    sbuf = (double *) malloc(nSmax * sizeof(double));

    for (p = 1; p < steps; p++) {
        rp = (rank + size - p) % size;
        if (nR[p] != 0)
            MPI_Irecv(rbuf, nR[p], MPI_DOUBLE, rp, tag, comm, &r_request);
        sp = (rank + p) % size;
        if (nS[p] != 0) {
            m2s(val, sbuf, S[p], nS[p]); // fill the sending buffer
            MPI_Isend(sbuf, nS[p], MPI_DOUBLE, sp, tag, comm, &s_request);
        }
        tag++;

        if (nR[p] != 0) {
            MPI_Wait(&r_request, MPI_STATUS_IGNORE);
            s2m_sum(res_val, rbuf, R[p],
                    nR[p]); // sum receive buffer into values
        }
        if (nS[p] != 0) MPI_Wait(&s_request, MPI_STATUS_IGNORE);
    }
    free(sbuf);
    free(rbuf);
    return 0;
}

#endif
