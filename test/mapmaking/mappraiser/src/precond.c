// MAPPRAISER preconditioner vdev
// Routines for computing the diagonal, block-diagonal Jacobi, and Two-level preconditioners for the PCG
// The routines also deal with degenerate pixels to ensure numerical stability of the system

/** @file   precond.c
        @author Hamza El Bouhargani
        @date   May 2019
        @credit  Adapted from work by Frederic Dauvergne
        @Last_update January 2021 by Aygul Jamal */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <assert.h>

// choose header based on compilation option
#ifdef W_MKL
#include <mkl.h>
#else
#include <lapacke.h>
#endif

#include "midapack.h"
#include "mappraiser.h"

#define eps 1.0e-15

// do the local Atdiag(Nm1)A with as output a block-diagonal matrix (stored as a vector) in the pixel domain
int getlocalW(const Mat *A, const Tpltz *Nm1, double *vpixBlock, int *lhits)
{
    int i, j, k, l;                                 // some indexes
    int m = Nm1->local_V_size;                       // number of local time samples
    int nnz = A->nnz;                               // number of non-zero entries
    int n_valid = A->lcount - (A->trash_pix) * nnz; // size of map-domain objects

    // Define the indices for each process
    int idv0, idvn; // indice of the first and the last block of V for each processes
    int *nnew;
    nnew = (int *)calloc(Nm1->nb_blocks_loc, sizeof(int));
    int64_t idpnew;
    int local_V_size_new;

    // get idv0 and idvn
    get_overlapping_blocks_params(Nm1->nb_blocks_loc, Nm1->tpltzblocks, Nm1->local_V_size, Nm1->nrow, Nm1->idp, &idpnew, &local_V_size_new, nnew, &idv0, &idvn);
    // double *vpixDiag;
    // vpixDiag = (double *) malloc(A->lcount *sizeof(double));

    int istart, il, istartn;
    for (i = 0; i < nnz * n_valid; i++)
        vpixBlock[i] = 0.0;

    int vShft = idpnew - Nm1->idp; //=Nm1->tpltzblocks[idv0].idv-Nm1->idp in principle
    /*
        printf("Nm1->idp=%d, idpnew=%d, vShft=%d\n", Nm1->idp, idpnew, vShft);
        printf("idv0=%d, idvn=%d\n", idv0, idvn);
        printf("Nm1->nb_blocks_loc=%d, Nm1->local_V_size=%d\n", Nm1->nb_blocks_loc, Nm1->local_V_size);

        for(i=0; i < Nm1->nb_blocks_loc; i++)
        printf("Nm1->tpltzblocks[%d].idv=%d\n", i, Nm1->tpltzblocks[i].idv);
    */

    // go until the first piecewise stationary period
    int vpix_i, vpix_j = 0;
    for (i = 0; i < vShft; i++)
    {
        // only do something for valid pixels
        if (!A->trash_pix || (A->trash_pix && (A->indices[i * nnz] != 0)))
        {
            vpix_i = A->indices[i * nnz] / nnz - A->trash_pix;
            lhits[vpix_i] += 1;
            for (j = 0; j < nnz; j++)
            {
                vpix_j = A->indices[i * nnz + j] - nnz * A->trash_pix;
                for (k = 0; k < nnz; k++)
                {
                    // value (j, k) of 3*3 block
                    vpixBlock[nnz * vpix_j + k] += A->values[i * nnz + j] * A->values[i * nnz + k];
                }
            }
        }
    }

    // temporary buffer for one diag value of Nm1
    double diagNm1;
    // loop on the blocks
    for (k = idv0; k < (idv0 + Nm1->nb_blocks_loc); k++)
    {
        if (nnew[idv0] > 0)
        { // if nnew==0, this is a wrong defined block

            if (k + 1 < idv0 + Nm1->nb_blocks_loc) // if there is a next block, compute his next first indice
                istartn = Nm1->tpltzblocks[k + 1].idv - Nm1->idp;
            else
                istartn = Nm1->local_V_size;
            // istartn = 0;

            istart = max(0, Nm1->tpltzblocks[k].idv - Nm1->idp);
            il = Nm1->tpltzblocks[k].n; // added this line to default code

            // if block cut from the left:
            if (k == idv0)
                il = min(Nm1->tpltzblocks[k].n, Nm1->tpltzblocks[k].idv + Nm1->tpltzblocks[k].n - Nm1->idp);
            // if block cut from the right:
            if (k == idv0 + Nm1->nb_blocks_loc - 1)
                il = min(il, (Nm1->idp + Nm1->local_V_size) - Nm1->tpltzblocks[k].idv);
            // if block alone in the middle, and cut from both sides
            if (Nm1->nb_blocks_loc == 1)
                il = min(il, Nm1->local_V_size);

            // get the diagonal value of the Toeplitz
            diagNm1 = Nm1->tpltzblocks[k].T_block[0];
            /*
            printf("istart=%d, il=%d, istartn=%d\n", istart, il, istartn);
            printf("Nm1->tpltzblocks[k=%d].idv=%d, Nm1->tpltzblocks[k=%d].n=%d, Nm1->idp=%d\n", k, Nm1->tpltzblocks[k].idv, k, Nm1->tpltzblocks[k].n, Nm1->idp);
            */

            // a piecewise stationary period
            vpix_i = 0;
            vpix_j = 0;
            for (i = istart; i < istart + il; i++)
            {
                // only do something for valid pixels
                if (!A->trash_pix || (A->trash_pix && (A->indices[i * nnz] != 0)))
                {
                    vpix_i = A->indices[i * nnz] / nnz - A->trash_pix;
                    lhits[vpix_i] += 1;
                    for (j = 0; j < nnz; j++)
                    {
                        vpix_j = A->indices[i * nnz + j] - nnz * A->trash_pix;
                        for (l = 0; l < nnz; l++)
                        {
                            vpixBlock[nnz * vpix_j + l] += A->values[i * nnz + j] * A->values[i * nnz + l] * diagNm1;
                        }
                    }
                }
            }

            // continue until the next period if exist or to the last line of V
            vpix_i = 0;
            vpix_j = 0;
            for (i = istart + il; i < istartn; i++)
            {
                // only do something for valid pixels
                if (!A->trash_pix || (A->trash_pix && (A->indices[i * nnz] != 0)))
                {
                    vpix_i = A->indices[i * nnz] / nnz - A->trash_pix;
                    lhits[vpix_i] += 1;
                    for (j = 0; j < nnz; j++)
                    {
                        vpix_j = A->indices[i * nnz + j] - nnz * A->trash_pix;
                        for (l = 0; l < nnz; l++)
                        {
                            vpixBlock[nnz * vpix_j + l] += A->values[i * nnz + j] * A->values[i * nnz + l];
                        }
                    }
                }
            }
        }
    } // end of the loop over the blocks

    return 0;
}

// do the local diag( At diag(Nm1) A ) with as output a vector in the pixel domain
//  This is an old deprecated routine
int getlocDiagN(Mat *A, Tpltz Nm1, double *vpixDiag)
{
    int i, j, k; // some indexes
    int m;

    m = Nm1.local_V_size; // number of local time samples
    //  int nnz=(A->nnz);

    // Define the indices for each process
    int idv0, idvn; // indice of the first and the last block of V for each processes
    int *nnew;
    nnew = (int *)calloc(Nm1.nb_blocks_loc, sizeof(int));
    int64_t idpnew;
    int local_V_size_new;
    // get idv0 and idvn
    get_overlapping_blocks_params(Nm1.nb_blocks_loc, Nm1.tpltzblocks, Nm1.local_V_size, Nm1.nrow, Nm1.idp, &idpnew, &local_V_size_new, nnew, &idv0, &idvn);
    // double *vpixDiag;
    // vpixDiag = (double *) malloc(A->lcount *sizeof(double));

    int istart, il, istartn;
    for (i = 0; i < A->lcount; i++)
        vpixDiag[i] = 0.0; // 0.0;

    int vShft = idpnew - Nm1.idp; //=Nm1.tpltzblocks[idv0].idv-Nm1.idp in principle
    /*
        printf("Nm1.idp=%d, idpnew=%d, vShft=%d\n", Nm1.idp, idpnew, vShft);
        printf("idv0=%d, idvn=%d\n", idv0, idvn);
        printf("Nm1.nb_blocks_loc=%d, Nm1.local_V_size=%d\n", Nm1.nb_blocks_loc, Nm1.local_V_size);

        for(i=0; i < Nm1.nb_blocks_loc; i++)
        printf("Nm1.tpltzblocks[%d].idv=%d\n", i, Nm1.tpltzblocks[i].idv);
    */

    // go until the first piecewise stationary period
    for (i = 0; i < vShft; i++)
    {
        for (j = 0; j < (A->nnz); j++)
            vpixDiag[A->indices[i * (A->nnz) + j]] += (A->values[i * (A->nnz) + j] * A->values[i * (A->nnz) + j]);
    }

    // temporary buffer for one diag value of Nm1
    int diagNm1;
    // loop on the blocks
    for (k = idv0; k < (idv0 + Nm1.nb_blocks_loc); k++)
    {
        if (nnew[idv0] > 0)
        { // if nnew==0, this is a wrong defined block

            if (k + 1 < idv0 + Nm1.nb_blocks_loc) // if there is a next block, compute his next first indice
                istartn = Nm1.tpltzblocks[k + 1].idv - Nm1.idp;
            else
                istartn = 0;

            istart = max(0, Nm1.tpltzblocks[k].idv - Nm1.idp);

            // if block cut from the left:
            if (k == idv0)
                il = min(Nm1.tpltzblocks[k].n, Nm1.tpltzblocks[k].idv + Nm1.tpltzblocks[k].n - Nm1.idp);
            // if block cut from the right:
            if (k == idv0 + Nm1.nb_blocks_loc - 1)
                il = min(il, (Nm1.idp + Nm1.local_V_size) - Nm1.tpltzblocks[k].idv);
            // if block alone in the middle, and cut from both sides
            if (Nm1.nb_blocks_loc == 1)
                il = min(il, Nm1.local_V_size);

            // get the diagonal value of the Toeplitz
            diagNm1 = Nm1.tpltzblocks[k].T_block[0];

            /*
            printf("istart=%d, il=%d, istartn=%d\n", istart, il, istartn);
            printf("Nm1.tpltzblocks[k].idv=%d, Nm1.tpltzblocks[k].n=%d, Nm1.idp=%d\n", Nm1.tpltzblocks[k].idv, Nm1.tpltzblocks[k].n, Nm1.idp);
            */
            // a piecewise stationary period
            for (i = istart; i < istart + il; i++)
            {
                for (j = 0; j < (A->nnz); j++)
                    vpixDiag[A->indices[i * (A->nnz) + j]] += (A->values[i * (A->nnz) + j] * A->values[i * (A->nnz) + j]) * diagNm1;
            }

            // continue until the next period if exist
            for (i = istart + il; i < istartn; i++)
            {
                for (j = 0; j < (A->nnz); j++)
                    vpixDiag[A->indices[i * (A->nnz) + j]] += (A->values[i * (A->nnz) + j] * A->values[i * (A->nnz) + j]);
            }
        }
    } // end of the loop over the blocks

    return 0;
}

// communication scheme in the pixel domain for the vector vpixDiag
// extract from a Madmap routine
int commScheme(Mat *A, double *vpixDiag, int pflag)
{
    int i, j, k;
    int nSmax, nRmax, nStot, nRtot;
    double *lvalues, *com_val, *out_val;

    int nbr_values = A->lcount - (A->nnz) * (A->trash_pix);

#if W_MPI
    lvalues = (double *)malloc(nbr_values * sizeof(double)); /*<allocate and set to 0.0 local values*/
    memcpy(lvalues, vpixDiag, nbr_values * sizeof(double));  /*<copy local values into result values*/

    nRmax = 0;
    nSmax = 0;

    if (A->flag == BUTTERFLY)
    { /*<branch butterfly*/
        // memcpy(out_values, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
        for (k = 0; k < A->steps; k++) /*compute max communication buffer size*/
            if (A->nR[k] > nRmax)
                nRmax = A->nR[k];
        for (k = 0; k < A->steps; k++)
            if (A->nS[k] > nSmax)
                nSmax = A->nS[k];

        com_val = (double *)malloc(A->com_count * sizeof(double));
        for (i = 0; i < A->com_count; i++)
        {
            com_val[i] = 0.0;
        }
        // already done    memcpy(vpixDiag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
        m2m(lvalues, A->lindices + (A->nnz) * (A->trash_pix), nbr_values, com_val, A->com_indices, A->com_count);
        butterfly_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
        m2m(com_val, A->com_indices, A->com_count, vpixDiag, A->lindices + (A->nnz) * (A->trash_pix), nbr_values);
        free(com_val);
    }
    else if (A->flag == BUTTERFLY_BLOCKING_1)
    {
        for (k = 0; k < A->steps; k++) // compute max communication buffer size
            if (A->nR[k] > nRmax)
                nRmax = A->nR[k];
        for (k = 0; k < A->steps; k++)
            if (A->nS[k] > nSmax)
                nSmax = A->nS[k];
        com_val = (double *)malloc(A->com_count * sizeof(double));
        for (i = 0; i < A->com_count; i++)
            com_val[i] = 0.0;
        m2m(lvalues, A->lindices + (A->nnz) * (A->trash_pix), nbr_values, com_val, A->com_indices, A->com_count);
        butterfly_blocking_1instr_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
        m2m(com_val, A->com_indices, A->com_count, vpixDiag, A->lindices + (A->nnz) * (A->trash_pix), nbr_values);
        free(com_val);
    }
    else if (A->flag == BUTTERFLY_BLOCKING_2)
    {
        for (k = 0; k < A->steps; k++) // compute max communication buffer size
            if (A->nR[k] > nRmax)
                nRmax = A->nR[k];
        for (k = 0; k < A->steps; k++)
            if (A->nS[k] > nSmax)
                nSmax = A->nS[k];
        com_val = (double *)malloc(A->com_count * sizeof(double));
        for (i = 0; i < A->com_count; i++)
            com_val[i] = 0.0;
        m2m(lvalues, A->lindices + (A->nnz) * (A->trash_pix), nbr_values, com_val, A->com_indices, A->com_count);
        butterfly_blocking_1instr_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
        m2m(com_val, A->com_indices, A->com_count, vpixDiag, A->lindices + (A->nnz) * (A->trash_pix), nbr_values);
        free(com_val);
    }
    else if (A->flag == NOEMPTYSTEPRING)
    {
        for (k = 1; k < A->steps; k++) // compute max communication buffer size
            if (A->nR[k] > nRmax)
                nRmax = A->nR[k];
        nSmax = nRmax;
        ring_noempty_step_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, vpixDiag, A->steps, A->comm);
    }
    else if (A->flag == RING)
    {
        // already done    memcpy(vpixDiag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
        for (k = 1; k < A->steps; k++) /*compute max communication buffer size*/
            if (A->nR[k] > nRmax)
                nRmax = A->nR[k];

        nSmax = nRmax;
        ring_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, vpixDiag, A->steps, A->comm);
    }
    else if (A->flag == NONBLOCKING)
    {
        // already done    memcpy(vpixDiag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
        ring_nonblocking_reduce(A->R, A->nR, A->S, A->nS, lvalues, vpixDiag, A->steps, A->comm);
    }
    else if (A->flag == NOEMPTY)
    {
        // already done    memcpy(vpixDiag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
        int ne = 0;
        for (k = 1; k < A->steps; k++)
            if (A->nR[k] != 0)
                ne++;

        ring_noempty_reduce(A->R, A->nR, ne, A->S, A->nS, ne, lvalues, vpixDiag, A->steps, A->comm);
    }
    else if (A->flag == ALLREDUCE)
    {
        com_val = (double *)malloc(A->com_count * sizeof(double));
        out_val = (double *)malloc(A->com_count * sizeof(double));
        for (i = 0; i < A->com_count; i++)
        {
            com_val[i] = 0.0;
            out_val[i] = 0.0;
        }
        s2m(com_val, lvalues, A->com_indices, nbr_values);
        /*for(i=0; i < A->com_count; i++){
            printf("%lf ", com_val[i]);
            } */
        MPI_Allreduce(com_val, out_val, A->com_count, MPI_DOUBLE, MPI_SUM, A->comm); // maximum index
        /*for(i=0; i < A->com_count; i++){
            printf("%lf ", out_val[i]);
            } */
        m2s(out_val, vpixDiag, A->com_indices, nbr_values); // sum receive buffer into values
        free(com_val);
        free(out_val);
    }
    else if (A->flag == ALLTOALLV)
    {
        nRtot = nStot = 0;
        for (k = 0; k < A->steps; k++)
        {                      // compute buffer sizes
            nRtot += A->nR[k]; // to receive
            nStot += A->nS[k]; // to send
        }
        alltoallv_reduce(A->R, A->nR, nRtot, A->S, A->nS, nStot, lvalues, vpixDiag, A->steps, A->comm);
    }
    else
    {
        printf("\n\n####### ERROR! : Incorrect communication scheme #######\n\n");
        exit(1);
    }
#endif
    free(lvalues);
    return 0;
}

/** @brief Compute Diag(A' diag(Nm1) A). This an old deprecated routine
        @param out_values local output array of doubles*/
int DiagAtA(Mat *A, double *diag, int pflag)
{
    int i, j, k;
    int nSmax, nRmax;
    double *lvalues;

    lvalues = (double *)malloc(A->lcount * sizeof(double)); /*<allocate and set to 0.0 local va
                                    lues*/
    for (i = 0; i < A->lcount; i++)
        lvalues[i] = 0.0;

    // Naive computation with a full defined diag(Nm1):
    for (i = 0; i < A->m; i++)
        for (j = 0; j < A->nnz; j++)                                                                              /*<dot products */
            lvalues[A->indices[i * (A->nnz) + j]] += (A->values[i * (A->nnz) + j] * A->values[i * (A->nnz) + j]); //*vdiagNm1[i];

#if W_MPI
    nRmax = 0;
    nSmax = 0;

    if (A->flag == BUTTERFLY)
    { /*<branch butterfly*/
        // memcpy(out_values, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
        for (k = 0; k < A->steps; k++) /*compute max communication buffer size*/
            if (A->nR[k] > nRmax)
                nRmax = A->nR[k];
        for (k = 0; k < A->steps; k++)
            if (A->nS[k] > nSmax)
                nSmax = A->nS[k];

        double *com_val;
        com_val = (double *)malloc(A->com_count * sizeof(double));
        for (i = 0; i < A->com_count; i++)
        {
            com_val[i] = 0.0;
        }
        memcpy(diag, lvalues, (A->lcount) * sizeof(double)); /*<copy local values into result values*/
        m2m(lvalues, A->lindices, A->lcount, com_val, A->com_indices, A->com_count);
        butterfly_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
        m2m(com_val, A->com_indices, A->com_count, diag, A->lindices, A->lcount);
        free(com_val);
    }
    else if (A->flag == RING)
    {
        memcpy(diag, lvalues, (A->lcount) * sizeof(double)); /*<copy local values into result values*/
        for (k = 1; k < A->steps; k++)                       /*compute max communication buffer size*/
            if (A->nR[k] > nRmax)
                nRmax = A->nR[k];

        nSmax = nRmax;
        ring_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, diag, A->steps, A->comm);
    }
    else if (A->flag == NONBLOCKING)
    {
        memcpy(diag, lvalues, (A->lcount) * sizeof(double)); /*<copy local values into result values*/
        ring_nonblocking_reduce(A->R, A->nR, A->S, A->nS, lvalues, diag, A->steps, A->comm);
    }
    else if (A->flag == NOEMPTY)
    {
        memcpy(diag, lvalues, (A->lcount) * sizeof(double)); /*<copy local values into result values*/
        int ne = 0;
        for (k = 1; k < A->steps; k++)
            if (A->nR[k] != 0)
                ne++;

        ring_noempty_reduce(A->R, A->nR, ne, A->S, A->nS, ne, lvalues, diag, A->steps, A->comm);
    }
    else
    {
        return 1;
    }
#endif
    free(lvalues);
    return 0;
}

int get_pixshare_pond(Mat *A, double *pixpond)
{

    int i, j, k; // some indexes
    int m, n, rank, size;
    double localreduce; // reduce buffer
    double st, t;       // timers

    MPI_Comm_rank(A->comm, &rank); //
    MPI_Comm_size(A->comm, &size); //

    m = A->m;                                  // number of local time samples
    n = A->lcount - (A->nnz) * (A->trash_pix); // number of local pixels

    // create an eyes local vector
    for (i = 0; i < n; i++)
        pixpond[i] = 1.;

    // communicate with the others processes to have the global reduce
    commScheme(A, pixpond, 2);

    // compute the inverse vector
    for (i = 0; i < n; i++)
        pixpond[i] = 1. / pixpond[i];

    return 0;
}

// Block diagonal jacobi preconditioner with degenerate pixels pre-processing
int precondblockjacobilike(Mat *A, Tpltz *Nm1, Mat *BJ_inv, Mat *BJ, double *b, double *cond, int *lhits)
{
    int i, j, k; // some indexes
    int m, m_cut, n, rank, size, nnz;
    int *indices_new, *tmp1;
    double *vpixBlock, *vpixBlock_inv, *vpixBlock_loc, *hits_proc, *tmp2, *tmp3, *tmp4;
    double det, invdet;
    int info, nb, lda;
    double anorm, rcond;

    int iw[3];
    double w[18];
    double x[9];
    nb = 3;
    lda = 3;

    m = A->m;
    nnz = A->nnz;
    m_cut = m;
    // at this stage there can be some invalid (flagged) pixels
    n = (A->lcount) - nnz * (A->trash_pix);
    MPI_Comm_rank(A->comm, &rank); //
    MPI_Comm_size(A->comm, &size);

    vpixBlock_loc = (double *)malloc(n * sizeof(double));
    vpixBlock = (double *)malloc(n * nnz * sizeof(double));
    vpixBlock_inv = (double *)malloc(n * nnz * sizeof(double));
    hits_proc = (double *)malloc(n * sizeof(double));
    if (indices_new == NULL || vpixBlock_loc == NULL || vpixBlock == NULL || vpixBlock_inv == NULL || hits_proc == NULL)
    {
        printf("Out of memory: memory allocation failed");
        exit(1);
    }

    // Compute local Atdiag(N^1)A
    getlocalW(A, Nm1, vpixBlock, lhits);

    // sum hits globally: dumb chunk of code that needs to be rewritten, but works for now ...
    for (i = 0; i < n; i += 3)
    {
        hits_proc[i] = lhits[(int)i / 3];
        hits_proc[i + 1] = lhits[(int)i / 3];
        hits_proc[i + 2] = lhits[(int)i / 3];
    }
    commScheme(A, hits_proc, 2);
    for (i = 0; i < n; i += 3)
    {
        lhits[(int)i / 3] = (int)hits_proc[i];
    }
    free(hits_proc);

    // communicate with the other processes to have the global reduce
    // TODO : This should be done in a more efficient way
    for (i = 0; i < n * nnz; i += nnz * nnz)
    {
        for (j = 0; j < nnz; j++)
        {
            vpixBlock_loc[(i / nnz) + j] = vpixBlock[i + j];
        }
    }
    commScheme(A, vpixBlock_loc, 2);
    for (i = 0; i < n * nnz; i += nnz * nnz)
    {
        for (j = 0; j < nnz; j++)
        {
            vpixBlock[i + j] = vpixBlock_loc[(i / nnz) + j];
        }
    }

    for (i = 3; i < n * nnz; i += nnz * nnz)
    {
        for (j = 0; j < nnz; j++)
            vpixBlock_loc[(i - 3) / nnz + j] = vpixBlock[i + j];
    }
    commScheme(A, vpixBlock_loc, 2);
    for (i = 3; i < n * nnz; i += nnz * nnz)
    {
        for (j = 0; j < nnz; j++)
            vpixBlock[i + j] = vpixBlock_loc[(i - 3) / nnz + j];
    }

    for (i = 6; i < n * nnz; i += nnz * nnz)
    {
        for (j = 0; j < nnz; j++)
            vpixBlock_loc[(i - 6) / nnz + j] = vpixBlock[i + j];
    }
    commScheme(A, vpixBlock_loc, 2);
    for (i = 6; i < n * nnz; i += nnz * nnz)
    {
        for (j = 0; j < nnz; j++)
            vpixBlock[i + j] = vpixBlock_loc[(i - 6) / nnz + j];
    }
    free(vpixBlock_loc);

    // Compute the inverse of the global Atdiag(N^-1)A blocks (beware this part is only valid for nnz = 3)
    int uncut_pixel_index = 0;
    for (i = 0; i < n * nnz; i += nnz * nnz)
    {
        // init 3x3 block
        double block[3][3];
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
            {
                block[j][k] = vpixBlock[i + (j * 3) + k];
                x[k + 3 * j] = block[j][k];
            }
        }

        // Compute the reciprocal of the condition number of the block

        /* Computes the norm of x */
        anorm = dlange_("1", &nb, &nb, x, &lda, w);

        /* Modifies x in place with a LU decomposition */
        dgetrf_(&nb, &nb, x, &lda, iw, &info);
        // if (info != 0) fprintf(stderr, "failure with error %d\n", info);

        /* Computes the reciprocal norm */
        dgecon_("1", &nb, x, &lda, &anorm, &rcond, w, iw, &info);
        // if (info != 0) fprintf(stderr, "failure with error %d\n", info);

        cond[(int)i / 9] = rcond;

        // Compute det
        // TODO: This should take into account the fact that the blocks are symmetric
        //  and generalize beyond the 3x3 case
        det = block[0][0] * (block[1][1] * block[2][2] - block[2][1] * block[1][2]) -
              block[0][1] * (block[1][0] * block[2][2] - block[1][2] * block[2][0]) +
              block[0][2] * (block[1][0] * block[2][1] - block[1][1] * block[2][0]);

        // Remove the degenerate pixels from the map-making

        if (rcond > 1e-1)
        {
            // The pixel is well enough observed, inverse the preconditioner block

            invdet = 1 / det;

            // Compute the inverse coeffs
            // TODO: This should take into account the fact that the blocks are symmetric
            // and generalize beyond the 3x3 case
            vpixBlock_inv[i] = (block[1][1] * block[2][2] - block[2][1] * block[1][2]) * invdet;
            vpixBlock_inv[i + 1] = (block[0][2] * block[2][1] - block[0][1] * block[2][2]) * invdet;
            vpixBlock_inv[i + 2] = (block[0][1] * block[1][2] - block[0][2] * block[1][1]) * invdet;
            vpixBlock_inv[i + 3] = (block[1][2] * block[2][0] - block[1][0] * block[2][2]) * invdet;
            vpixBlock_inv[i + 4] = (block[0][0] * block[2][2] - block[0][2] * block[2][0]) * invdet;
            vpixBlock_inv[i + 5] = (block[1][0] * block[0][2] - block[0][0] * block[1][2]) * invdet;
            vpixBlock_inv[i + 6] = (block[1][0] * block[2][1] - block[2][0] * block[1][1]) * invdet;
            vpixBlock_inv[i + 7] = (block[2][0] * block[0][1] - block[0][0] * block[2][1]) * invdet;
            vpixBlock_inv[i + 8] = (block[0][0] * block[1][1] - block[1][0] * block[0][1]) * invdet;
        }
        else
        {
            // The pixel is not well enough observed
            // Remove the poorly conditioned pixel from the map
            // and point the associated gap samples to trash pixel

            // Raise pointing matrix flag to signal presence of trash pixels
            A->trash_pix = 1;

            // Search for the corresponding gap samples
            j = A->id_last_pix[uncut_pixel_index]; // last index of time sample pointing to degenerate pixel
            
            // Point the last gap sample to trash pixel
            for (k = 0; k < nnz; k++)
            {
                A->indices[j * nnz + k] = k - nnz;
                A->values[j * nnz + k] = 0;
            }

            // Set the corresponding signal time stream sample to zero
            b[j] = 0;

            // Point all the preceding gap samples to trash pixel and set them to zero in the TOD
            while (A->ll[j] != -1)
            {
                b[A->ll[j]] = 0;
                for (k = 0; k < nnz; k++)
                {
                    A->indices[A->ll[j] * nnz + k] = k - nnz;
                    A->values[A->ll[j] * nnz + k] = 0;
                }
                j = A->ll[j];
            }

            // Remove degenerate pixel from vpixBlock, lhits, and cond
            memmove(vpixBlock + i, vpixBlock + i + nnz * nnz, (n * nnz - nnz * nnz - i) * sizeof(double));
            memmove(lhits + (int)i / (nnz * nnz), lhits + (int)i / (nnz * nnz) + 1, ((int)n / nnz - 1 - (int)i / (nnz * nnz)) * sizeof(int));
            memmove(cond + (int)i / (nnz * nnz), cond + (int)i / (nnz * nnz) + 1, ((int)n / nnz - 1 - (int)i / (nnz * nnz)) * sizeof(double));

            // Shrink effective size of vpixBlock
            n -= nnz;
            i -= nnz * nnz;
        }
        ++uncut_pixel_index;
    }
    // free memory
    free(A->id_last_pix);
    free(A->ll);

    // Reallocate memory for preconditioner blocks and redefine pointing matrix in case of the presence of degenerate pixels
    if (A->trash_pix)
    {
        // Reallocate memory of vpixBlock by shrinking its memory size to its effective size (no degenerate pixel)
        tmp2 = (double *)realloc(vpixBlock, n * nnz * sizeof(double));
        tmp4 = (double *)realloc(vpixBlock_inv, n * nnz * sizeof(double));
        tmp1 = (int *)realloc(lhits, (int)n / nnz * sizeof(int));
        tmp3 = (double *)realloc(cond, (int)n / nnz * sizeof(double));
        if (tmp1 == NULL || tmp2 == NULL || tmp3 == NULL || tmp4 == NULL)
        {
            printf("[rank %d] Out of memory: reallocation of preconditioner blocks, hits or rcond maps failed", rank);
            exit(1);
        }
        else
        {
            vpixBlock = tmp2;
            vpixBlock_inv = tmp4;
            lhits = tmp1;
            cond = tmp3;
        }
    }

    // map local indices to global indices in indices_cut
    for (i = 0; i < m * nnz; i++)
    {
        // switch to global indices
        // except degenerate pixels, which have index < 0
        if (A->indices[i] >= 0)
            A->indices[i] = A->lindices[A->indices[i]];
    }

    // free  memory of original pointing matrix and synchronize
    MatFree(A);

    // Define new pointing matrix (marginalized over degenerate pixels and free from gap samples)
    MatInit(A, m, nnz, A->indices, A->values, A->flag, MPI_COMM_WORLD);

    // Define Block-Jacobi preconditioner indices
    indices_new = (int *)malloc(n * nnz * sizeof(int));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < nnz; j++)
        {
            indices_new[i * nnz + j] = A->lindices[nnz * (A->trash_pix) + nnz * ((int)i / nnz) + j];
        }
    }

    // Init Block-Jacobi inv preconditioner
    MatSetIndices(BJ_inv, n, nnz, indices_new);
    MatSetValues(BJ_inv, n, nnz, vpixBlock_inv);
    BJ_inv->trash_pix = 0;
    MatLocalShape(BJ_inv, 3);

    // Init Block-Jacobi preconditioner
    MatSetIndices(BJ, n, nnz, indices_new);
    MatSetValues(BJ, n, nnz, vpixBlock);
    BJ->trash_pix = 0;
    MatLocalShape(BJ, 3);

    return 0;
}

int precondjacobilike_avg(Mat A, Tpltz Nm1, double *c)
{

    int i, j, k; // some indexes
    int m, n, rank, size;
    double localreduce; // reduce buffer
    double st, t;       // timers

    m = A.m;                      // number of local time samples
    n = A.lcount;                 // number of local pixels
    MPI_Comm_rank(A.comm, &rank); //
    MPI_Comm_size(A.comm, &size); //

    double diagNm1;

    // Compute diag( AtA )
    DiagAtA(&A, c, 2);

    // multiply by the diagonal Toeplitz
    diagNm1 = Nm1.tpltzblocks[0].T_block[0];

    printf("diagNm1 = %f \n", diagNm1);
    for (j = 0; j < n; j++)
        c[j] = diagNm1 * c[j];

    // compute c inverse vector
    for (j = 0; j < n; j++)
        c[j] = 1. / c[j];

    return 0;
}

int precondjacobilike(Mat A, Tpltz Nm1, int *lhits, double *cond, double *vpixDiag)
{

    int i, j, k; // some indexes
    int m, n, rank, size;
    double localreduce; // reduce buffer
    double st, t;       // timers

    MPI_Comm_rank(A.comm, &rank); //
    MPI_Comm_size(A.comm, &size); //

    m = A.m;      // number of local time samples
    n = A.lcount; // number of local pixels

    // Compute local diag( At diag(N^1) A )
    getlocDiagN(&A, Nm1, vpixDiag);

    // communicate with the other processes to have the global reduce
    commScheme(&A, vpixDiag, 2);

    // compute the inverse vector
    for (i = 0; i < n; i++)
    {
        if (i % 3 == 0)
        {
            lhits[(int)i / 3] = (int)vpixDiag[i];
            cond[(int)i / 3] = vpixDiag[i + 1] + vpixDiag[i + 2];
        }
        vpixDiag[i] = 1. / vpixDiag[i];
    }
    return 0;
}

// low-level routine: to be moved somwhere else
void transpose_nn(double *A, int n)
{
    int i, j;
    double temp;

    for (i = 0; i < n - 1; i++)
        for (j = i + 1; j < n; j++)
        {
            temp = A[i * n + j];
            A[i * n + j] = A[j * n + i];
            A[j * n + i] = temp;
        }
}

void inverse_svd(int m, int n, int lda, double *a)
{

    int info = 0;
    int i = 0;
    int rank = 0;
    // lapack_int LAPACKE_dgelss( int matrix_order, lapack_int m, lapack_int n, lapack_int nrhs, double* a, lapack_int lda, double* b, lapack_int ldb, double* s, double rcond, lapack_int* rank );
    double *b = calloc(m * n, sizeof(double));
    int nsv = m < n ? m : n;
    double *s = malloc(nsv * sizeof(double));

    for (i = 0; i < nsv; i++)
        b[i * n + i] = 1;

    info = LAPACKE_dgelss(LAPACK_ROW_MAJOR, m, n, n, a, n, b, n, s, eps, &rank);

    if (info != 0)
        printf("LAPACK_dgelss does not work.\n");
    memcpy(a, b, (m * n) * sizeof(double));
    free(b);
}

// Deflation subspace matrix constructor
void build_Z(const Mat *A, int Zn, double ***out_Z)
{
    int i, j, g, e, k, rank, size, group;
    int p, rp, sp, tag = 0;
    MPI_Request s_request, r_request;
    MPI_Status status;

    int lcount_max = 0, rlcount;
    int *count, *tcount;
    int *rcount, *rindices;
    double *rZ;
    double **Z; // Zn * A->lcount, pointers to columns (ie column-major)

    MPI_Comm_rank(A->comm, &rank); // get rank and size of the communicator
    MPI_Comm_size(A->comm, &size);

    // Compute lcount_max (all procs)
    MPI_Allreduce(&(A->lcount), &(lcount_max), 1, MPI_INT, MPI_MAX, A->comm);

    // If number of columns of Z >= number of processes, each process
    // will compute a group of columns of Z instead of a single column
    if (Zn > size)
    {
        group = Zn / size;
        assert(group * size == Zn);
    }
    else
    {
        group = 1;
    }

    // Allocate buffers
    count = (int *)calloc(group * A->lcount, sizeof(int)); // the number of appereance in local processor
    tcount = (int *)calloc(A->lcount, sizeof(int));
    rcount = (int *)malloc(group * lcount_max * sizeof(int)); // the number of appereance in neighbor processors
    rindices = (int *)malloc(lcount_max * sizeof(int));       // real indices in neighbor processors

    // Compute local count for a given pixel
    for (g = 0; g < group; g++)
        for (i = 0; i < A->m / group * A->nnz; i++)
            count[g * A->lcount + A->indices[g * (A->m / group * A->nnz) + i]]++;

    // Copy local count in total count
    for (g = 0; g < group; g++)
        for (i = 0; i < A->lcount; i++)
            tcount[i] += count[g * A->lcount + i];

    // Compute total counts
    for (p = 1; p < size; p++)
    { // loop : collective global reduce in ring-like fashion
        rp = (size + rank - p) % size;
        sp = (rank + p) % size;
        MPI_Send(&(A->lcount), 1, MPI_INT, sp, 0, A->comm); // exchange sizes
        MPI_Recv(&rlcount, 1, MPI_INT, rp, 0, A->comm, &status);
        tag++;
        MPI_Irecv(rindices, rlcount, MPI_INT, rp, tag, A->comm, &r_request); // exchange global indices
        MPI_Isend(A->lindices, A->lcount, MPI_INT, sp, tag, A->comm, &s_request);
        MPI_Wait(&r_request, &status);
        MPI_Wait(&s_request, &status);
        tag++;
        MPI_Irecv(rcount, group * rlcount, MPI_INT, rp, tag, A->comm, &r_request); // exchange local count/rcount values
        MPI_Isend(count, group * A->lcount, MPI_INT, sp, tag, A->comm, &s_request);
        tag++;
        MPI_Wait(&r_request, &status);
        for (g = 0; g < group; g++)
            m2m_sum_i(rcount + g * rlcount, rindices, rlcount, tcount, A->lindices, A->lcount); // sum in the result
        MPI_Wait(&s_request, &status);
    }

    // Free no longer used buffers
    free(rcount);

    // Allocate Z
    Z = calloc(group * size, sizeof(double *));

    // Compute the current process' Z
    for (g = 0; g < group; g++)
    {
        Z[rank * group + g] = calloc(A->lcount, sizeof(double));
        for (i = 0; i < A->lcount; i += A->nnz)
            Z[rank * group + g][i] =
                (double)((double)count[g * group + i] / (double)tcount[i]);
    }

    // Free no longer used buffers
    free(count);
    free(tcount);

    // Allocate the buffer to exchange Z
    rZ = (double *)malloc(lcount_max * sizeof(double));

    // Exchange Z
    for (p = 1; p < size; p++)
    { // loop : collective global reduce in ring-like fashion
        rp = (size + rank - p) % size;
        sp = (rank + p) % size;
        MPI_Send(&(A->lcount), 1, MPI_INT, sp, 0, A->comm); // exchange sizes
        MPI_Recv(&rlcount, 1, MPI_INT, rp, 0, A->comm, &status);
        tag++;
        MPI_Irecv(rindices, rlcount, MPI_INT, rp, tag, A->comm, &r_request); // exchange global indices
        MPI_Isend(A->lindices, A->lcount, MPI_INT, sp, tag, A->comm, &s_request);
        MPI_Wait(&r_request, &status);
        MPI_Wait(&s_request, &status);
        tag++;

        for (g = 0; g < group; g++)
        {
            MPI_Irecv(rZ, rlcount, MPI_DOUBLE, rp, tag, A->comm, &r_request); // exchange local values
            MPI_Isend(Z[rank * group + g], A->lcount, MPI_DOUBLE, sp, tag, A->comm, &s_request);
            tag++;
            MPI_Wait(&r_request, &status);
            Z[rp * group + g] = calloc(A->lcount, sizeof(double));
            m2m(rZ, rindices, rlcount, Z[rp * group + g], A->lindices, A->lcount); // copy the interesting value the corresponding Z[rp] value
        }
        MPI_Wait(&s_request, &status);
    }

    double **ZZ;

    // If number of columns of Z < number of processes, shrink Z
    if (Zn < size)
    {
        int ratio = size / Zn;
        assert(Zn * ratio == size);

        ZZ = calloc(Zn, sizeof(double *));
        for (j = 0; j < Zn; j++)
        {
            ZZ[j] = Z[j * ratio];
            for (k = 1; k < ratio; k++)
            {
                for (i = 0; i < A->lcount; i++)
                {
                    ZZ[j][i] += Z[j * ratio + k][i];
                }
                free(Z[j * ratio + k]);
            }
        }
        free(Z);

        // Otherwise, the result is just Z
    }
    else
    {
        ZZ = Z;
    }

    for (j = 0; j < Zn; j++)
    {
        for (i = 0; i < A->lcount - A->nnz * A->trash_pix; i++)
        {
            ZZ[j][i] = ZZ[j][i + (A->nnz) * (A->trash_pix)];
        }
    }

    // Free no longer used buffers
    free(rindices);
    free(rZ);

    *out_Z = ZZ;
}

// in vector x overlapped,
// in P distributed among processors
// Z should be contructed in place
// Calculation of E, then E^{-1}, we have the algorithm
// Pt N^{-1} P Z (Zt Pt N^{-1} P Z)^{-1} Zt x = A Z E^{-1} Zt x = A Q x
// v1 = P * Zi column by column
// v2 = N^{-1} * v1
// v3 = Pt * v2
// Ei = Zt * v3, Zt is entire matrix and Ei is computed column by column, Zi is a column

void build_Em1(const Mat *A, double **Z, double **AZ, const double *pixpond, int Zn, int n, double **out_E)
{
    int i, j, k, rank;
    double *E, *EO;

    MPI_Comm_rank(A->comm, &rank);

    E = calloc(Zn * Zn, sizeof(double));

    // Ei = Zt * Pt * N^{-1} * P * Zi
    // Em1 = E^{-1}

    for (i = 0; i < Zn; i++)
    {
        // E = Zt * v3
        for (j = 0; j < Zn; j++)
        {
            for (k = 0; k < n; k++)
            {
                E[i * Zn + j] += (double)(Z[j][k] * AZ[i][k] * pixpond[k]);
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, E, Zn * Zn, MPI_DOUBLE, MPI_SUM, A->comm);

    int info;
    double anorm, rcond, *w;
    int *iw;
    iw = calloc(Zn, sizeof(int));
    w = calloc(Zn * Zn * 2, sizeof(double));

    inverse_svd(Zn, Zn, Zn, E);

    // EO = calloc(Zn * Zn, sizeof(double));
    //
    // memcpy(EO, E, sizeof(double) * Zn * Zn);
    //
    // /* Computes the norm of x */
    // anorm = dlange_("1", &Zn, &Zn, EO, &Zn, w);
    //
    // /* Modifies x in place with a LU decomposition */
    // dgetrf_(&Zn, &Zn, EO, &Zn, iw, &info);
    // // if (info != 0) fprintf(stderr, "failure with error %d\n", info);
    //
    // /* Computes the reciprocal norm */
    // dgecon_("1", &Zn, EO, &Zn, &anorm, &rcond, w, iw, &info);
    // // if (info != 0) fprintf(stderr, "failure with error %d\n", info);
    //
    // //printf("condition number of Einv = %25.18e\n", rcond);
    //
    // free(EO);

    *out_E = E;
}

// AZ constructor to speed-up iterations with the 2lvl preconditioners
void build_AZ(Mat *A, const Tpltz *Nm1, double **Z, int Zn, int n, double ***out_AZ)
{
    double **AZ;
    double *v;
    int i, j, k;

    v = calloc(A->m, sizeof(double)); // P * Zi

    AZ = calloc(Zn, sizeof(double *));
    for (k = 0; k < Zn; k++)
    {
        AZ[k] = calloc(n, sizeof(double));
    }

    for (i = 0; i < Zn; i++)
    {
        // AZ[i] = Pt N^{-1} P * Z[i]
        MatVecProd(A, Z[i], v, 0);
        stbmmProd(Nm1, v); // In-place
        TrMatVecProd(A, v, AZ[i], 0);
    }

    free(v);

    *out_AZ = AZ;
}

void build_Qtx(const Mat *A, double **Z, const double *Em1, const double *x, const double *pixpond, double *w, double *Qtx, int Zn, int n)
{
    int i, j, k;
    int rank;

    MPI_Comm_rank(A->comm, &rank);

    // Pt N^{-1} P Z (Zt Pt N^{-1} P Z)^{-1} Zt x
    // Pt N^{-1} P Z (E)^{-1} Zt x
    // Pt N^{-1} P Q x

    // w = Zt x
    for (j = 0; j < Zn; j++)
    {
        w[j] = 0.0;
        for (k = 0; k < n; k++)
        {
            w[j] += Z[j][k] * x[k] * pixpond[k];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, w, Zn, MPI_DOUBLE, MPI_SUM, A->comm);

    // Qtx = Em1 * w (dense matrix E (Zn*Zn) times a vector w (Zn*1)
    for (i = 0; i < Zn; i++)
    {
        Qtx[i] = 0.0;
        for (j = 0; j < Zn; j++)
        {
            Qtx[i] += Em1[i * Zn + j] * w[j];
        }
    }
}

void mul_ZQtx(double **Z, const double *Qtx, double *vec, int Zn, int n)
{
    int k, j;

    // vec = Z * Qtx (overlapped times dense);
    // for (k = 0; k < n; k++) {
    //   vec[k] = 0.0;
    //   for (j = 0; j < Zn; j++) {
    //     vec[k] += Z[j][k] * Qtx[j];
    //   }
    // }
    for (k = 0; k < n; k++)
        vec[k] = 0.0;

    for (j = 0; j < Zn; j++)
    {
        for (k = 0; k < n; k++)
        {
            vec[k] += Z[j][k] * Qtx[j];
        }
    }
}

// Lanczos procedure to build the deflation subspace of the "a posteriori" 2lvl preconditioner
void Lanczos_eig(Mat *A, const Tpltz *Nm1, const Mat *BJ_inv, const Mat *BJ, double *x, const double *b, const double *noise, double tol, const double *pixpond, int K, double ***out_Ritz_vectors, double ***out_Ritz_vectors_AZ)
{
    int i, j, k; // some indexes
    int m, n, rank, size;
    double st, t; // timers
    double solve_time = 0.0;
    double beta, alpha, result, dot;
    int info = 0, lwork = -1;

    double *Av = NULL, *_g = NULL;
    double *Tt = NULL, *T = NULL;
    double *w = NULL, *v = NULL, *vold = NULL;
    double *V = NULL, *AmulV = NULL;

    double *Ritz_values = NULL;
    double *Ritz_vectors_out = NULL;
    double *Ritz_vectors_out_r = NULL;
    double **Ritz_vectors = NULL;
    double **Ritz_vectors_AZ = NULL;

    double *work = NULL;
    double wkopt = 0.0;

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);

    m = A->m;                                    // Number of local time samples
    n = (A->lcount) - (A->nnz) * (A->trash_pix); // Number of local pixels

    st = MPI_Wtime();

    // Map domain
    Av = (double *)malloc(n * sizeof(double));
    _g = (double *)malloc(m * sizeof(double));

    T = (double *)calloc((K + 1) * (K + 1), sizeof(double));
    Tt = (double *)calloc(K * K, sizeof(double));

    Ritz_vectors_out = (double *)calloc(n * K, sizeof(double));
    // Ritz_vectors_out_r = (double *)calloc(n*K,  sizeof(double));
    w = (double *)malloc(n * sizeof(double));
    v = (double *)calloc(n, sizeof(double));
    vold = (double *)calloc(n, sizeof(double));
    V = (double *)calloc(n * (K + 1), sizeof(double));

    AmulV = (double *)calloc(n * (K + 1), sizeof(double));

    Ritz_values = (double *)calloc(K, sizeof(double));

    Ritz_vectors = calloc(K, sizeof(double *));
    for (i = 0; i < K; i++)
        Ritz_vectors[i] = calloc(n, sizeof(double));

    Ritz_vectors_AZ = calloc(K, sizeof(double *));
    for (i = 0; i < K; i++)
        Ritz_vectors_AZ[i] = calloc(n, sizeof(double));

    for (i = 0; i < n; ++i)
        w[i] = 0.0;

    MatVecProd(A, x, _g, 0);

    for (i = 0; i < m; i++)
        _g[i] = b[i] + noise[i] - _g[i];

    stbmmProd(Nm1, _g);

    TrMatVecProd(A, _g, w, 0);

    // beta = sqrt(dot(w, w))
    dot = 0.0;
    for (i = 0; i < n; i++)
        dot += w[i] * w[i] * pixpond[i];
    MPI_Allreduce(&dot, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    beta = sqrt(result);

    if (beta > eps)
    {
        for (i = 0; i < n; i++)
        {
            v[i] = w[i] / beta;
            V[i * (K + 1)] = v[i];
        }
    }

    // Av = A * v = Pt N P * v
    MatVecProd(A, v, _g, 0);
    stbmmProd(Nm1, _g);
    TrMatVecProd(A, _g, Av, 0);

    memcpy(w, Av, n * sizeof(double));

    // alpha = dot(v, Av)
    dot = 0.0;
    for (i = 0; i < n; i++)
        dot += v[i] * Av[i] * pixpond[i];
    MPI_Allreduce(&dot, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (i = 0; i < n; i++)
    {
        AmulV[i * (K + 1)] = w[i];
        w[i] = w[i] - (alpha * v[i]);
    }

    t = MPI_Wtime();
    if (rank == 0)
    {
        printf("[rank %d] Lanczos init time=%lf \n", rank, t - st);
        fflush(stdout);
    }

    st = MPI_Wtime();

    for (i = 0; i < K; i++)
    {

        dot = 0.0;
        for (j = 0; j < n; j++)
            dot += w[j] * w[j] * pixpond[j];
        MPI_Allreduce(&dot, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta = sqrt(result);

        if (beta > eps)
        {
            for (j = 0; j < n; j++)
            {
                v[j] = w[j] / beta;
                V[j * (K + 1) + i + 1] = v[j];
            }
        }

        else if (rank == 0)
            printf("division by zero in iteration %d\n", i);

        // What should we do to construct this special triangular matrix
        // T(i,i) = alpha
        // T(i,i+1) = beta
        // T(i+1,i) = beta

        T[(i * (K + 1)) + i] = alpha;
        T[(i * (K + 1)) + i + 1] = beta;
        T[((i + 1) * (K + 1)) + i] = beta;

        // Av = A * v = Pt N P * v
        MatVecProd(A, v, _g, 0);
        stbmmProd(Nm1, _g);
        TrMatVecProd(A, _g, Av, 0);

        memcpy(w, Av, n * sizeof(double));

        // alpha = dot(v, Av)
        dot = 0.0;
        for (j = 0; j < n; j++)
            dot += v[j] * Av[j] * pixpond[j];
        MPI_Allreduce(&dot, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        for (j = 0; j < n; j++)
        {
            AmulV[j * (K + 1) + i + 1] = w[j];
            w[j] = w[j] - (alpha * v[j]) - (beta * vold[j]);
            vold[j] = v[j];
        }

        t = MPI_Wtime();
        if (rank == 0)
        {
            printf("Iteration = %d, [rank %d] Lanczos iteration time=%lf \n", i, rank, t - st);
            fflush(stdout);
        }

        st = MPI_Wtime();
    }

    // Here we reduce the dimention of T from (K+1 * K+1) to K * K;
    // Here we reduce the dimention of V from (N * K+1) to (N * K)

    for (i = 0; i < K; i++)
    {
        for (j = 0; j < K; j++)
        {
            T[i * K + j] = T[i * (K + 1) + j];
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < K; j++)
        {
            AmulV[i * K + j] = AmulV[i * (K + 1) + j];
        }
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < K; j++)
        {
            V[i * K + j] = V[i * (K + 1) + j];
        }
    }

    //[Ritz_vectors, Ritz_values] = eig(T)
    // Ritz_values contains the eigenvalues of the matrix A in ascending order.

    // lapack_int LAPACKE_dsyev(int matrix_layout, char jobz, char uplo, lapack_int n, double* a, lapack_int lda, double* w);

    st = MPI_Wtime();

    transpose_nn(T, K);

    // int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info)

    dsyev_("Vectors", "Upper", &K, T, &K, Ritz_values, &wkopt, &lwork, &info);

    lwork = (int)wkopt;
    work = (double *)malloc(lwork * sizeof(double));

    dsyev_("Vectors", "Upper", &K, T, &K, Ritz_values, work, &lwork, &info);

    // free(work);

    transpose_nn(T, K);

    t = MPI_Wtime();
    if (rank == 0)
    {
        printf("[rank %d] Lanczos dsyev time=%lf \n", rank, t - st);
        fflush(stdout);
    }

    st = MPI_Wtime();

    memset(Ritz_vectors_out, 0, n * K * sizeof(double));

    for (i = 0; i < n; i++)
    {
        for (k = 0; k < K; k++)
        {
            for (j = 0; j < K; j++)
            {
                Ritz_vectors_out[i * K + j] += V[i * K + k] * T[k * K + j];
            }
        }
    }

    for (i = 0; i < n; i++)
        for (j = 0; j < K; j++)
            Ritz_vectors[j][i] = Ritz_vectors_out[i * K + j];

    memset(Ritz_vectors_out, 0, n * K * sizeof(double));

    for (i = 0; i < n; i++)
    {
        for (k = 0; k < K; k++)
        {
            for (j = 0; j < K; j++)
            {
                Ritz_vectors_out[i * K + j] += AmulV[i * K + k] * T[k * K + j];
            }
        }
    }

    for (i = 0; i < n; i++)
        for (j = 0; j < K; j++)
            Ritz_vectors_AZ[j][i] = Ritz_vectors_out[i * K + j];

    t = MPI_Wtime();
    if (rank == 0)
    {
        printf("[rank %d] Lanczos V*T multiplication time=%lf \n", rank, t - st);
        fflush(stdout);
    }

    free(Av);
    free(_g);
    free(Tt);
    free(T);
    free(w);
    free(v);
    free(vold);
    free(V);
    free(Ritz_values);
    free(Ritz_vectors_out);

    *out_Ritz_vectors = Ritz_vectors;
    *out_Ritz_vectors_AZ = Ritz_vectors_AZ;
}

// General routine for constructing a preconditioner
void build_precond(struct Precond **out_p, double **out_pixpond, int *out_n, Mat *A, Tpltz *Nm1, double **in_out_x, double *b, const double *noise, double *cond, int *lhits, double tol, int Zn, int precond)
{
    int rank, size, i;
    double st, t;
    double *x;
    // double *t1, *t2;

    struct Precond *p = calloc(1, sizeof(struct Precond));

    MPI_Comm_rank(A->comm, &rank);
    MPI_Comm_size(A->comm, &size);

    p->precond = precond;
    p->Zn = Zn;

    precondblockjacobilike(A, Nm1, &(p->BJ_inv), &(p->BJ), b, cond, lhits);

    p->n = (A->lcount) - (A->nnz) * (A->trash_pix);

    // Reallocate memory for well-conditioned map
    x = realloc(*in_out_x, p->n * sizeof(double));
    if (x == NULL)
    {
        printf("Out of memory: map reallocation failed");
        exit(1);
    }

    p->pixpond = (double *)malloc(p->n * sizeof(double));

    // Compute pixel share ponderation
    get_pixshare_pond(A, p->pixpond);

    if (precond != 0)
    {
        p->Qg = calloc(p->n, sizeof(double));
        p->AQg = calloc(p->n, sizeof(double));
        p->w = calloc(Zn, sizeof(double));   // Zt * x
        p->Qtx = calloc(Zn, sizeof(double)); // Em1 * Zt * x

        st = MPI_Wtime();

        // 2lvl a priori
        if (precond == 1)
        {
            build_Z(A, Zn, &(p->Z));
            build_AZ(A, Nm1, p->Z, Zn, p->n, &(p->AZ));
        }

        // 2lvl a posteriori
        else if (precond == 2)
            Lanczos_eig(A, Nm1, &(p->BJ_inv), &(p->BJ), x, b, noise, tol, p->pixpond, Zn, &(p->Z), &(p->AZ)); // 2lvl a posteriori preconditioner

        // Invalid precond
        else
        {
            printf("Whoops! Incorrect preconditioner parameter, please check documentation and try again: %d\n", precond);
            printf("Quick reminder: precond = 0 -> BD, precond = 1 -> 2lvl a priori, precond = 2 -> 2lvl a posteriori\n");
            exit(1);
        }

        t = MPI_Wtime();
        if (rank == 0)
        {
            printf("[rank %d] Deflation (Z) & projected deflation (AZ) subspaces computation time=%lf \n", rank, t - st);
            fflush(stdout);
        }
        st = MPI_Wtime();

        build_Em1(A, p->Z, p->AZ, p->pixpond, Zn, p->n, &(p->Em1));

        t = MPI_Wtime();
        if (rank == 0)
        {
            printf("[rank %d] 2lvl Em1 computation time=%lf \n", rank, t - st);
            fflush(stdout);
        }

        // For x initial guess is 0
        // memcpy(x, p->Z[Zn-1], p->n * sizeof(double));

        // For x initial guess is not 0
        // for (i = 0; i < p->n; i++)
        //  x[i] = x[i] +  p->Z[Zn-1][i];
    }

    *out_p = p;
    *out_pixpond = p->pixpond;
    *out_n = p->n;
    *in_out_x = x;
}

// General routine for applying the preconditioner to a map vector
void apply_precond(struct Precond *p, const Mat *A, const Tpltz *Nm1, double *g, double *Cg)
{
    int i;

    // BJ
    if (p->precond == 0)
    {
        MatVecProd(&(p->BJ_inv), g, Cg, 0);
    }

    // 2lvl
    else
    {
        build_Qtx(A, p->Z, p->Em1, g, p->pixpond, p->w, p->Qtx, p->Zn, p->n);
        mul_ZQtx(p->Z, p->Qtx, p->Qg, p->Zn, p->n);
        mul_ZQtx(p->AZ, p->Qtx, p->AQg, p->Zn, p->n);

        for (i = 0; i < p->n; ++i)
        {
            p->AQg[i] = g[i] - p->AQg[i];
        }

        MatVecProd(&(p->BJ_inv), p->AQg, Cg, 0);

        for (i = 0; i < p->n; ++i)
        {
            Cg[i] += p->Qg[i];
        }
    }
}

// Free memory from the preconditioner
void free_precond(struct Precond **in_out_p)
{
    int i;
    struct Precond *p = *in_out_p;

    free(p->pixpond);
    free(p->BJ_inv.indices);
    free(p->BJ_inv.values);
    MatFree(&(p->BJ_inv));

    // free(p->BJ.indices); // Shared with p->BJ_inv.indices
    free(p->BJ.values);
    MatFree(&(p->BJ));

    if (p->precond != 0)
    {
        for (i = 0; i < p->Zn; i++)
            free(p->Z[i]);
        free(p->Z);

        for (i = 0; i < p->Zn; i++)
            free(p->AZ[i]);
        free(p->AZ);

        free(p->Em1);
        free(p->Qg);
        free(p->AQg);
        free(p->w);
        free(p->Qtx);
    }

    free(p);
    *in_out_p = NULL;
}
