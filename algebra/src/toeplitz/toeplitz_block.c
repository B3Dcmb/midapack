/**
@file toeplitz_block.c version 1.1b, July 2012
@brief Contains routines related to the Toeplitz blocks diagonal routine for
Toeplitz algebra
@author  Frederic Dauvergne, Radek Stompor
**
** Project:  Midapack library, ANR MIDAS'09 - Toeplitz Algebra module
** Purpose:  Provide Toeplitz algebra tools suitable for Cosmic Microwave
Background (CMB)
**           data analysis.
**
***************************************************************************
@note Copyright (c) 2010-2012 APC CNRS Université Paris Diderot
@note
@note This program is free software; you can redistribute it and/or modify it
under the terms
@note of the GNU Lesser General Public License as published by the Free Software
Foundation;
@note either version 3 of the License, or (at your option) any later version.
This program is
@note distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even
@note the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU
@note Lesser General Public License for more details.
@note
@note You should have received a copy of the GNU Lesser General Public License
along with this
@note program; if not, see http://www.gnu.org/licenses/lgpl.html
@note For more information about ANR MIDAS'09 project see :
@note http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
@note ACKNOWLEDGMENT: This work has been supported in part by the French
National Research
@note Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
***************************************************************************
** Log: toeplitz*.c
**
** Revision 1.0b  2012/05/07  Frederic Dauvergne (APC)
** Official release 1.0beta. The first installement of the library is the
Toeplitz algebra
** module.
**
** Revision 1.1b  2012/07/-  Frederic Dauvergne (APC)
** - mpi_stbmm allows now rowi-wise order per process datas and no-blocking
communications.
** - OMP improvment for optimal cpu time.
** - bug fixed for OMP in the stmm_basic routine.
** - distcorrmin is used to communicate only lambda-1 datas when it is needed.
** - new reshaping routines using transformation functions in stmm. Thus, only
one copy
**   at most is needed.
** - tpltz_init improvement using define_nfft and define_blocksize routines.
** - add Block struture to define each Toeplitz block.
** - add Flag structure and preprocessing parameters to define the computational
strategy.
**   All the flag parameters are then available directly from the API.
**
***************************************************************************
**
*/

#include "toeplitz.h"

#include <stdlib.h>

#define max(a, b)               \
    ({                          \
        __typeof__(a) _a = (a); \
        __typeof__(b) _b = (b); \
        _a > _b ? _a : _b;      \
    })

#define min(a, b)               \
    ({                          \
        __typeof__(a) _a = (a);                                                \
        __typeof__(b) _b = (b);                                                \
        _a < _b ? _a : _b;                                                     \
    })

// r1.1 - Frederic Dauvergne (APC)
// This is the routines related to the Toeplitz blocks diagonal routine.
// There is a sequential equivalent routine in the file toeplitz_seq.c

// todo:
//- remove the nooptimize communication

//=========================================================================
#ifdef W_MPI
/// Performs the multiplication of a symmetric, Toeplitz block-diagonal matrix,
/// T, by an arbitrary matrix, V, distributed over processes in the generalized
/// column-wise way.
/** @ingroup group12
    Each process performs the multiplication sequentially for each diagonal
   block and based on the sliding window algorithm. Prior to that MPI calls are
   used to exchange data between neighboring process. Each of the diagonal
   blocks is a symmetric, band-diagonal Toeplitz matrix, which can be different
   for each block. The parameters are : \param V \b [input] distributed data
   matrix (with the convention V(i,j)=V[i+j*n]) ; \b [out] result of the product
   TV \param nrow number of rows of the global data matrix V \param m number of
   columns for the data matrix V in the global rowwise order \param m_rowwise
   number of columns for the data matrix V in the rowwise order per processor
    \param tpltzblocks list of the toeplitz blocks struture with its own
   parameters (idv, n, T_block, lambda) :
    - idv is the global row index defining for each Toeplitz block as stored in
   the vector T ;
    - n size of each Toeplitz block
    - T_block location of each Toeplitz matrix data composed of the non-zero
   entries of the first row ;
    - lambda size of each Toeplitz block data T_block. The bandwith size is then
   equal to lambda*2-1 \param nb_blocks_all number of all Toeplitz block on the
   diagonal of the full Toeplitz matrix \param nb_blocks_local number of
   Toeplitz blocks as stored in T \param idp global index of the first element
   of the local part of V \param local_V_size a number of all elements in local
   V \param flag_stgy flag strategy for the product computation \param comm MPI
   communicator
*/
int mpi_stbmm(double **V, int64_t nrow, int m, int m_rowwise,
              Block *tpltzblocks, int nb_blocks_local, int nb_blocks_all,
              int64_t idp, int local_V_size, Flag flag_stgy, MPI_Comm comm) {
#else // for sequential use only
int mpi_stbmm(double **V, int64_t nrow, int m, int m_rowwise,
              Block *tpltzblocks, int nb_blocks_local, int nb_blocks_all,
              int64_t idp, int local_V_size, Flag flag_stgy) {
#endif


    // MPI parameters
    int rank; // process rank
    int size; // process number

#ifdef W_MPI
    MPI_Status status;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

#else
    rank = 0;
    size = 1;
#endif

    PRINT_RANK = rank;

    FILE *file;
    file = stdout;

    int i, j, k; // some indexes


    // identification of the mpi neighbours process to communicate when there is
    // a shared block
    int right = rank + 1;
    int left  = rank - 1;


    // Define the indices for each process
    int idv0, idvn; // indice of the first and the last block of V for each
                    // processes

    int *nnew;
    nnew = (int *) calloc(nb_blocks_local, sizeof(int));
    int64_t idpnew;
    int     local_V_size_new;
    int     n_rowwise = local_V_size;

    int status_params = get_overlapping_blocks_params(
            nb_blocks_local, tpltzblocks, local_V_size, nrow, idp, &idpnew,
            &local_V_size_new, nnew, &idv0, &idvn);


    if (PRINT_RANK == 0 && VERBOSE > 2)
        printf("status_params=%d\n", status_params);

    if (status_params == 0) {
        free(nnew);
        return (0); // no work to be done
    }

    if (tpltzblocks[idv0].lambda == 0 || tpltzblocks[idvn].lambda == 0)
        return print_error_message(2, __FILE__, __LINE__);


    if (PRINT_RANK == 0 && VERBOSE > 2) { // print on screen news parameters
                                          // definition if VERBOSE
        fprintf(file, "new parameters caracteristics:\n");
        fprintf(file, "[%d] idp=%ld ; idpnew=%ld\n", rank, idp, idpnew);
        fprintf(file, "[%d] local_V_size=%d ; local_V_size_new=%d\n", rank,
                local_V_size, local_V_size_new);
        for (i = 0; i < nb_blocks_local; i++)
            fprintf(file, "[%d] n[%d]=%d ; nnew[%d]=%d\n", rank, i,
                    (tpltzblocks[i].n), i, nnew[i]);
        for (i = 0; i < nb_blocks_local; i++)
            fprintf(file, "[%d] tpltzblocks[%d].idv=%ld\n", rank, i,
                    tpltzblocks[i].idv);
    }


    int vShft = idpnew - idp; // new first element of relevance in V

    // Define the column indices:
    // index of the first and the last column of V for the current process
    int idvm0 = idpnew / nrow;
    int idvmn = (idpnew + local_V_size_new - 1) / nrow;
    // number of columns of V for the current process
    int ncol_rank = idvmn - idvm0 + 1;
    // number of blocks for the current process with possibly repetitions
    int nb_blocks_rank;

    if (ncol_rank == 1) // Empty process not allowed
        nb_blocks_rank = idvn - idv0 + 1;
    else
        nb_blocks_rank =
                (ncol_rank - 2) * nb_blocks_local + (nb_blocks_local - idv0)
                + (idvn + 1); // in this case nb_blocks_local = nblocs_all

    if (PRINT_RANK == 0 && VERBOSE > 2)
        fprintf(file, "[%d] nb_blocks_rank=%d, nb_blocks_local=%d\n", rank,
                nb_blocks_rank, nb_blocks_local);

    // Define the indices for the first and the last element in each blocks
    int idvp0 = idpnew % nrow
              - tpltzblocks[idv0].idv; // index of the first element of the
                                       // process in the first block
    int idvpn; // reverse index of the last element of the process in the last
               // block
    // It's the number of remaining elements needed to fully complete the last
    // block
    idvpn = tpltzblocks[idvn].idv + nnew[idvn] - 1
          - (idpnew + local_V_size_new - 1) % nrow;


    // Define the offsets for the first and last blocks of the process for V1
    int offset0, offsetn;
    int distcorrmin_idv0 = (tpltzblocks[idv0].lambda) - 1;
    int distcorrmin_idvn = (tpltzblocks[idvn].lambda) - 1;

    // if(idvp0 != 0)
    offset0 = min(idvp0, distcorrmin_idv0);
    // if(idvpn != 0)
    offsetn = min(idvpn, distcorrmin_idvn);


    int toSendLeft  = 0;
    int toSendRight = 0;

#ifdef W_MPI
    if (offset0 != 0) {
        toSendLeft = min(tpltzblocks[idv0].idv + nnew[idv0] - idpnew % nrow,
                         distcorrmin_idv0);
    }
    if (offsetn != 0) {
        toSendRight =
                min((idpnew + local_V_size_new) % nrow - tpltzblocks[idvn].idv,
                    distcorrmin_idvn);
    }

    int flag_optimlambda = 1; // to allocate only the memory place needed

    int     lambdaOut_offset;
    int     lambdaIn_offset;
    double *LambdaOut;
    int     lambdaOut_size, lambdaIn_size;

    if (flag_optimlambda == 1) {
        LambdaOut = (double *) calloc((toSendLeft + toSendRight) * m_rowwise,
                                      sizeof(double));
        lambdaOut_offset = toSendLeft * m_rowwise;
        lambdaIn_offset  = offset0 * m_rowwise;
        lambdaOut_size   = (toSendLeft + toSendRight) * m_rowwise;
        lambdaIn_size    = (offset0 + offsetn) * m_rowwise;
    } else {
        LambdaOut = (double *) calloc(
                (tpltzblocks[idv0].lambda + tpltzblocks[idvn].lambda)
                        * m_rowwise,
                sizeof(double));
        lambdaOut_offset = tpltzblocks[idv0].lambda * m_rowwise;
        lambdaIn_offset  = tpltzblocks[idv0].lambda * m_rowwise;
        lambdaOut_size   = (tpltzblocks[idv0].lambda + tpltzblocks[idvn].lambda)
                       * m_rowwise;
        lambdaIn_size = (tpltzblocks[idv0].lambda + tpltzblocks[idvn].lambda)
                      * m_rowwise;
    }


    if (offset0 != 0) {
        for (j = 0; j < m_rowwise; j++)
            for (i = 0; i < toSendLeft; i++)
                LambdaOut[i + j * toSendLeft] =
                        (*V)[i + j * n_rowwise]; // good because toSendLeft=0 if
                                                 // it
    } // doesnt start on a the first block.
    if (offsetn != 0) {
        for (j = 0; j < m_rowwise; j++)
            for (i = 0; i < toSendRight; i++)
                LambdaOut[i + j * toSendRight + lambdaOut_offset] =
                        (*V)[i + j * n_rowwise + local_V_size - toSendRight];
    } // good too using same argument than for offset0!=0
    // if local_V_size!=local_V_size_new+vShft mean there is extra
    // terms a the end and so offsetn=0
    // idpnew+local_V_size_new = idp+local_V_size and vShft=idpnew-idp
    // so local_V_size=vShft+local_V_size_new
    if (rank == 0 || offset0 == 0) left = MPI_PROC_NULL;
    if (rank == size - 1 || offsetn == 0) right = MPI_PROC_NULL;

    double *LambdaIn = (double *) calloc(lambdaIn_size, sizeof(double));


    int         flag_blockingcomm = 0; // to use blocking comm
    MPI_Request requestLeft_r, requestLeft_s;
    MPI_Request requestRight_r, requestRight_s;

    if (flag_blockingcomm == 1) {
        // send and receive data
        //! [communication blocking example]
        // to the Left
        MPI_Sendrecv(LambdaOut, toSendLeft * m_rowwise, MPI_DOUBLE, left,
                     MPI_USER_TAG, (LambdaIn + lambdaIn_offset),
                     offsetn * m_rowwise, MPI_DOUBLE, right, MPI_USER_TAG, comm,
                     &status);

        // to the Right
        MPI_Sendrecv((LambdaOut + lambdaOut_offset), toSendRight * m_rowwise,
                     MPI_DOUBLE, right, MPI_USER_TAG, LambdaIn,
                     offset0 * m_rowwise, MPI_DOUBLE, left, MPI_USER_TAG, comm,
                     &status);
        //! [communication blocking example]
    } else {
        //! [communication non-blocking example]
        // to the Left
        MPI_Irecv((LambdaIn + lambdaIn_offset), offsetn * m_rowwise, MPI_DOUBLE,
                  right, MPI_USER_TAG, comm, &requestLeft_r);
        MPI_Isend(LambdaOut, toSendLeft * m_rowwise, MPI_DOUBLE, left,
                  MPI_USER_TAG, comm, &requestLeft_s);

        // to the Right
        MPI_Irecv(LambdaIn, offset0 * m_rowwise, MPI_DOUBLE, left, MPI_USER_TAG,
                  comm, &requestRight_r);
        MPI_Isend((LambdaOut + lambdaOut_offset), toSendRight * m_rowwise,
                  MPI_DOUBLE, right, MPI_USER_TAG, comm, &requestRight_s);
        //! [communication non-blocking example]
    }

#endif


    // size of the first and the last block for the current process
    int v0rank_size, vnrank_size;
    if (nb_blocks_rank == 1) { // only one block
        v0rank_size = ((idpnew + local_V_size_new - 1) % nrow + 1)
                    - idpnew % nrow + offset0 + offsetn;
        vnrank_size = 0; // just for convenience - no really need it
    } else {             // more than one block
        v0rank_size =
                tpltzblocks[idv0].idv + nnew[idv0] - idpnew % nrow + offset0;
        vnrank_size = ((idpnew + local_V_size_new - 1) % nrow + 1)
                    - tpltzblocks[idvn].idv + offsetn;
    }


#ifdef W_MPI

    if (flag_blockingcomm != 1) {
        // MPI_Wait for lambda comm
        //! [communication Wait example]
        MPI_Wait(&requestLeft_r, &status);
        MPI_Wait(&requestLeft_s, &status);
        MPI_Wait(&requestRight_r, &status);
        MPI_Wait(&requestRight_s, &status);
        //! [communication Wait example]
    }


    free(LambdaOut);

#endif


    //---------------------------------------
    // initialization for the blocks loop

    int idv1 = 0; // old index of *V1
    int idv2 = 0; // index


    int mid; // local number of column for the current block
    // index of the first element of the process inside the first block
    int offset_id0;
    offset_id0 = idvp0;

    // fftw variables
    fftw_complex *V_fft, *T_fft;
    double       *V_rfft;
    fftw_plan     plan_f, plan_b;
    // init local block vector
    double *V1block;
    //  int lambdaShft;


    //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // loop on the blocks inside the process
    int nfft, blocksize;
    int iblock; // index for the loop on the blocks
    //  int loopindex;
    int id; // indice of the current block

    int vblock_size;
    int id0block;

    int jj;


    for (iblock = idv0; iblock < idv0 + nb_blocks_rank; iblock++) {
        id = iblock % nb_blocks_local; // index of current block


        if (nnew[id] > 0) { // the block is ok

#ifdef W_MPI
            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            // first case : First block of the process
            if (iblock == idv0) {
                if (PRINT_RANK == 0 && VERBOSE > 2)
                    fprintf(file, "[%d] First block...\n", rank);

                vblock_size = v0rank_size;
                id0block    = (offset_id0 - offset0);

                V1block = (double *) calloc(vblock_size * m_rowwise,
                                            sizeof(double));

                for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                    for (i = 0; i < offset0; i++)
                        V1block[i + j * vblock_size] =
                                LambdaIn[i + j * offset0];
                }
                // note: check if copyblock could be used instead.


                // if (nb_blocks_rank == 1)
                // currentsize_middlepart=vblock_size-offset0-offsetn =
                // local_V_size_new else
                // currentsize_middlepart=vblock_size-offset0
                int currentsize_middlepart =
                        min(vblock_size - offset0, local_V_size_new);

                for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                    for (i = 0; i < currentsize_middlepart; i++)
                        V1block[offset0 + i + j * vblock_size] =
                                (*V)[i + vShft + j * n_rowwise];
                }

                if (nb_blocks_rank == 1) {
                    for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                        for (i = 0; i < offsetn; i++) {
                            V1block[vblock_size - offsetn + i
                                    + j * vblock_size] =
                                    LambdaIn[i + lambdaIn_offset + j * offsetn];
                        }
                    }
                }


                // init Toeplitz arrays
                tpltz_init(vblock_size, tpltzblocks[id].lambda, &nfft,
                           &blocksize, &T_fft, (tpltzblocks[id].T_block),
                           &V_fft, &V_rfft, &plan_f, &plan_b, flag_stgy);

                // Toeplitz computation
                if (PRINT_RANK == 0 && VERBOSE > 2)
                    fprintf(file,
                            "[%d] Before stmm_main call : nfft = %d, blocksize "
                            "= %d\n",
                            rank, nfft, blocksize);
                stmm_main(&V1block, vblock_size, m_rowwise, 0,
                          m_rowwise * vblock_size, (tpltzblocks[id].T_block),
                          T_fft, tpltzblocks[id].lambda, V_fft, V_rfft, plan_f,
                          plan_b, blocksize, nfft, flag_stgy);

                tpltz_cleanup(&T_fft, &V_fft, &V_rfft, &plan_f, &plan_b);


                int currentsize = min(vblock_size - offset0, local_V_size_new);
                for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                    for (i = 0; i < currentsize; i++)
                        (*V)[vShft + i + j * n_rowwise] =
                                V1block[offset0 + i + j * vblock_size];
                }

                free(V1block);

            } // end (First case)


            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            // Generic case : Generic block of the process
            else if (iblock != idv0 && iblock != idv0 + nb_blocks_rank - 1) {
#endif

                if (PRINT_RANK == 0 && VERBOSE > 2)
                    fprintf(file, "[%d] generic block...\n", rank);

                vblock_size = nnew[id];
                id0block    = 0;

                V1block = (double *) calloc(vblock_size * m_rowwise,
                                            sizeof(double));

                idv1 = (tpltzblocks[id].idv) - idp % nrow - vShft + offset0
                     + nrow * ((iblock / nb_blocks_local)); // no need
                //  idv2 = idv[id]-idp%nrow + nrow*( (iblock/nb_blocks_local) );
                idv2 = (tpltzblocks[id].idv) - (idpnew) % nrow + vShft
                     + nrow * ((iblock / nb_blocks_local));

                for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                    for (i = 0; i < vblock_size; i++)
                        V1block[i + j * vblock_size] =
                                (*V)[i + idv2 + j * n_rowwise];
                    //    V1block[i] = (*V)[i+idv1-offset0+vShft];
                }

                // init Toeplitz arrays
                tpltz_init(nnew[id], tpltzblocks[id].lambda, &nfft, &blocksize,
                           &T_fft, (tpltzblocks[id].T_block), &V_fft, &V_rfft,
                           &plan_f, &plan_b, flag_stgy);

                // Toeplitz computation
                if (PRINT_RANK == 0 && VERBOSE > 2)
                    fprintf(file,
                            "[%d] Before stmm_main call : nfft = %d, blocksize "
                            "= %d\n",
                            rank, nfft, blocksize);
                stmm_main(&V1block, vblock_size, m_rowwise, 0,
                          m_rowwise * vblock_size, (tpltzblocks[id].T_block),
                          T_fft, tpltzblocks[id].lambda, V_fft, V_rfft, plan_f,
                          plan_b, blocksize, nfft, flag_stgy);


                tpltz_cleanup(&T_fft, &V_fft, &V_rfft, &plan_f, &plan_b);


                for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                    for (i = 0; i < vblock_size; i++) {
                        (*V)[i + idv2 + j * n_rowwise] =
                                V1block[i + j * vblock_size];
                    }
                }


                free(V1block);

#ifdef W_MPI
            } // end (Generic case)

              //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
            // Last case : Last block of the process
            else if (iblock == idv0 + nb_blocks_rank - 1 && iblock != idv0) {
                if (PRINT_RANK == 0 && VERBOSE > 2)
                    fprintf(file, "[%d] last block...\n", rank);

                vblock_size = vnrank_size;
                id0block    = 0;

                V1block = (double *) calloc(vblock_size * m_rowwise,
                                            sizeof(double));

                idv1 = (tpltzblocks[id].idv) - idp % nrow - vShft + offset0
                     + nrow * ((iblock / nb_blocks_local));
                idv2 = (tpltzblocks[id].idv) - (idpnew) % nrow + vShft
                     + nrow * ((iblock / nb_blocks_local));


                for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                    for (i = 0; i < vblock_size - offsetn; i++)
                        V1block[i + j * vblock_size] =
                                (*V)[i + idv2 + j * n_rowwise];
                    //    V1block[i] = (*V)[i+idv1-offset0+vShft];
                }

                for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                    for (i = 0; i < offsetn; i++)
                        V1block[vblock_size - offsetn + i + j * vblock_size] =
                                LambdaIn[i + lambdaIn_offset + j * offsetn];
                }


                // init Toeplitz arrays
                tpltz_init(vblock_size, tpltzblocks[id].lambda, &nfft,
                           &blocksize, &T_fft, (tpltzblocks[id].T_block),
                           &V_fft, &V_rfft, &plan_f, &plan_b, flag_stgy);

                // Toeplitz computation
                if (PRINT_RANK == 0 && VERBOSE > 2)
                    fprintf(file,
                            "[%d] Before stmm_main call : nfft = %d, blocksize "
                            "= %d\n",
                            rank, nfft, blocksize);

                stmm_main(&V1block, vblock_size, m_rowwise, 0,
                          vblock_size * m_rowwise, (tpltzblocks[id].T_block),
                          T_fft, tpltzblocks[id].lambda, V_fft, V_rfft, plan_f,
                          plan_b, blocksize, nfft, flag_stgy);

                tpltz_cleanup(&T_fft, &V_fft, &V_rfft, &plan_f, &plan_b);

                for (j = 0; j < m_rowwise; j++) {
#pragma omp parallel for // num_threads(NB_OMPTHREADS_STBMM)
                    for (i = 0; i < vnrank_size - offsetn; i++) {
                        (*V)[idv2 + i + j * n_rowwise] =
                                V1block[i + j * vblock_size];
                    }
                }


                free(V1block);

            } // end of last block
            else {
                break;
            } // error  //we can put the generic case here instead of between
              // first and last cases
#endif
            //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        } // end of if(nnew[id]>0)
    }     // end of loop over the blocks


    free(LambdaIn);


    return 0;
}

// #endif

//====================================================================

/// Select the useful flotting blocks for the local data of the current
/// processor. parameters (idpnew, local_V_size_new, nnew) for the computation.
/// ide the local range are set with a size nnew equal to zero.
/** @ingroup group22
    This compute the right parameters (idpnew, local_V_size_new, nnew) for the
   computation. All the block outside the local range are set with a size nnew
   equal to zero. local_V_size_new correspond to the size without the shift
   between the global rank index and the global index of the first flotting
   block. idnew is then set to the index of this first flotting block.
*/
int get_overlapping_blocks_params(int nbloc, Block *tpltzblocks,
                                  int local_V_size, int64_t nrow, int64_t idp,
                                  int64_t *idpnew, int *local_V_size_new,
                                  int *nnew, int *ifirstBlock,
                                  int *ilastBlock) {
    int     ib, nblockOK = 0, nfullcol_data;
    int64_t firstrow, lastrow;
    int64_t idptmp;


    // check how many full columns input data have
    nfullcol_data = max(0, (local_V_size - (nrow - idp % nrow) % nrow
                            - (idp + local_V_size) % nrow)
                                   / nrow);

    if (nfullcol_data > 0) {

        for (ib = 0; ib < nbloc; ib++) {
            if (tpltzblocks[ib].idv < nrow) {
                nnew[ib] = min(tpltzblocks[ib].n,
                               nrow - tpltzblocks[ib].idv); // block used for
                                                            // the product
                nblockOK++;
            }
        }

    } else { // no full column observed

        firstrow = idp % nrow;
        lastrow  = (idp + local_V_size - 1) % nrow;

        if (firstrow < lastrow) { // just one column partially observed

            for (ib = 0; ib < nbloc; ib++) {
                if ((tpltzblocks[ib].idv + tpltzblocks[ib].n > firstrow)
                    && (tpltzblocks[ib].idv < lastrow + 1)) {
                    nnew[ib] =
                            min(tpltzblocks[ib].n,
                                nrow - tpltzblocks[ib].idv); // block used for
                                                             // the product
                    nblockOK++;
                }
            }

        } else { // two columns partially observed

            for (ib = 0; ib < nbloc; ib++) {
                if ((tpltzblocks[ib].idv + tpltzblocks[ib].n > firstrow)
                    && (tpltzblocks[ib].idv
                        < nrow)) { // intersects first partial column
                    nnew[ib] =
                            min(tpltzblocks[ib].n,
                                nrow - tpltzblocks[ib].idv); // block used for
                                                             // the product
                    nblockOK++;
                }

                if ((tpltzblocks[ib].idv < lastrow + 1)
                    && (tpltzblocks[ib].idv + tpltzblocks[ib].n
                        > 0)) { // intersects second partial column
                    nnew[ib] =
                            min(tpltzblocks[ib].n,
                                nrow - tpltzblocks[ib].idv); // block used for
                                                             // the product
                    nblockOK++; // may overcount but we do not care
                }               // could use else insteed!
            }
        }
    }
    if (PRINT_RANK == 0 && VERBOSE > 2) printf("nblockOK=%d\n", nblockOK);


    if (nblockOK == 0) return (0); // no blocks overlapping with the data

    // find the first and last relevant blocks for the begining and end of the
    // local data  V

    // first block
    idptmp = idp;

    for (*ifirstBlock = -1; *ifirstBlock == -1;) {
        for (ib = 0; ib < nbloc; ib++) {
            if (nnew[ib] != 0 && idptmp % nrow < tpltzblocks[ib].idv + nnew[ib])
                break;
        }

        if (ib < nbloc && tpltzblocks[ib].idv <= idptmp % nrow) {
            *ifirstBlock = ib;
            *idpnew      = idptmp;
        } else if (ib < nbloc && tpltzblocks[ib].idv > idptmp % nrow) {
            *ifirstBlock = ib;
            //   int64_t extrabegining = tpltzblocks[ib].idv-idp%nrow;  //note I
            //   put int64 just to be sure. Never used
            //      *idpnew = idp+extrabegining;//tpltzblocks[ib].idv;
            int idvfirstcolumn = idptmp / nrow;
            *idpnew            = tpltzblocks[ib].idv + idvfirstcolumn * nrow;
        } else {                            // ib=nb_blocs
            idptmp += nrow - idptmp % nrow; //(int) (nrow-idptmp%nrow);
            //          idtmp = (int) ceil((1.0*idpnew)/(1.0*nrow))*nrow; // go
            //          to the first element of the next column
        }
    }


    // last block
    idptmp = idp + local_V_size - 1;

    for (*ilastBlock = -1; *ilastBlock == -1;) {
        for (ib = nbloc - 1; ib >= 0; ib--) {
            if (nnew[ib] != 0 && tpltzblocks[ib].idv <= idptmp % nrow) break;
        }


        if (ib >= 0 && idptmp % nrow < tpltzblocks[ib].idv + nnew[ib]) {
            *ilastBlock       = ib;
            *local_V_size_new = local_V_size - (*idpnew) + idp;
        } else if (ib >= 0 && tpltzblocks[ib].idv + nnew[ib] <= idptmp % nrow) {
            *ilastBlock = ib;
            // int64_t extraend =
            // (local_V_size-1+idp)%nrow+1-(tpltzblocks[ib].idv+nnew[ib]);
            // //note I put int64 just to be sure *local_V_size_new =
            //(local_V_size+idp)%nrow-(idv[*ilastBlock]+nnew[*ilastBlock]);
            // idv[*ilastBlock]+nnew[*ilastBlock]-(*idpnew);
            //*local_V_size_new = local_V_size-(*idpnew)+idp-extraend; //compute
            //twice ... ? remove this one

            int idvlastcolumn = idptmp / nrow;
            *local_V_size_new = tpltzblocks[ib].idv + nnew[ib]
                              + idvlastcolumn * nrow - (*idpnew);

        } else {
            idptmp = idptmp - (idptmp % nrow)
                   - 1; //(int) idptmp - (idptmp%nrow)-1;
            //        idtmp = (int) floor( (1.0*idpnew)/(1.0*nrow))*nrow-1; //
            //        go to the last element of the previous column
        }
    }

    return (1);
}
