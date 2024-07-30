#include "sys/param.h"
#include "sys/stat.h"
#include "sys/types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef W_MPI
#include <mpi.h>

#include "mapmat/alm.h"
#include "mapmat/bitop.h"
#include "mapmat/butterfly.h"
#include "mapmat/butterfly_new.h"
#include "mapmat/butterfly_wrappers.h"
#include "mapmat/csort.h"
#include "mapmat/mpi_utils.h"

int mirror_butterfly(double *values_local, int *indices_local, int size_local,
                     double *values_received, int *indices_received,
                     int *size_received,
                     int flag_mirror_unmirror_size_indices_data,
                     MPI_Comm worldcomm) {
    /* Sends data of the excess processes, which are over the 2^k processes
       which will be used for the Butterfly scheme, back to the 2^k processes.
        The "mirror" step corresponds to sending the information of the excess
       processes (over 2^k) to the Butterfly processes, while the "unmirror"
       step corresponds to sending back the information of the Butterfly
       processes to the excess processes.

        The flag flag_mirror_unmirror_size_indices_data should be either
       MIRROR_SIZE, MIRROR_INDICES, MIRROR_DATA, UNMIRROR_SIZE, UNMIRROR_INDICES
       or UNMIRROR_DATA to respectively mirror/unmirror the size of the data,
       indices of the data, or the data itself

        The way this proceed is by taking the data of the N-2^k processes and
       send them respectively to the last N-2^k processes which will be used for
       the Butterfly scheme. This means the data of the process ranked 2^k will
       be send to the process ranked 2^k-1, and the data of the process ranked
       2^k+10 will be send to the process ranked 2^k-11 This function should be
       followed by reorder_indices to rearrange the indices and associated
       values, to make sure the received indices are order monotonously and the
       redundant indices and values are treated correctly */

    int rank, size;
    MPI_Comm_size(worldcomm, &size);
    MPI_Comm_rank(worldcomm, &rank);
    int i;
    int number_steps_Butterfly = log_2(size);
    int nb_rank_Butterfly = pow_2(number_steps_Butterfly);
    int *list_rank_excess_butterfly, *list_last_ranks_within_butterfly;

    int number_ranks_to_send = size - nb_rank_Butterfly;
    int tag = 100;

    if (number_ranks_to_send == 0)
        return 0;

    if (number_ranks_to_send) {
        list_rank_excess_butterfly =
            (int *)malloc(number_ranks_to_send * sizeof(int));
        list_last_ranks_within_butterfly =
            (int *)malloc(number_ranks_to_send * sizeof(int));
    }
    for (i = 0; i < number_ranks_to_send; i++) {
        list_rank_excess_butterfly[i] = nb_rank_Butterfly + i;
        list_last_ranks_within_butterfly[i] = nb_rank_Butterfly - i - 1;
    }

    int index_excess = elem_in_list_elem(rank, list_rank_excess_butterfly,
                                         number_ranks_to_send);
    int index_within = elem_in_list_elem(rank, list_last_ranks_within_butterfly,
                                         number_ranks_to_send);

    if ((index_excess == -1 && index_within == -1)) {
        if (number_ranks_to_send) {
            free(list_rank_excess_butterfly);
            free(list_last_ranks_within_butterfly);
        }
        return 0;
    }

    // First, communicate the sizes, determine which processes will receive/send
    // data
    int size_communicated = 0;
    switch (flag_mirror_unmirror_size_indices_data) {
    case MIRROR_SIZE:
        // Case 0 : MIRROR - send only size and indices to be exchanged, to have
        // the excess processes over 2^k to send their data to the last 2^k
        // processes before the Butterfly scheme
        mpi_send_data_from_list_rank(
            list_rank_excess_butterfly, list_last_ranks_within_butterfly,
            number_ranks_to_send, &size_local, 1 * sizeof(int),
            &size_communicated, tag, worldcomm);
        *size_received = size_communicated;
        break;

    case MIRROR_INDICES:
        size_communicated = *size_received;
        if (index_excess != -1)
            size_communicated = size_local;

        mpi_send_data_from_list_rank(list_rank_excess_butterfly,
                                     list_last_ranks_within_butterfly,
                                     number_ranks_to_send, indices_local,
                                     size_communicated * sizeof(int),
                                     indices_received, tag + 1, worldcomm);
        break;

    case MIRROR_DATA:
        size_communicated = *size_received;
        if (index_excess != -1)
            size_communicated = size_local;
        mpi_send_data_from_list_rank(list_rank_excess_butterfly,
                                     list_last_ranks_within_butterfly,
                                     number_ranks_to_send, values_local,
                                     size_communicated * sizeof(double),
                                     values_received, tag + 2, worldcomm);
        break;

    case UNMIRROR_SIZE:
        mpi_send_data_from_list_rank(
            list_last_ranks_within_butterfly, list_rank_excess_butterfly,
            number_ranks_to_send, &size_local, 1 * sizeof(int),
            &size_communicated, tag + 3, worldcomm);
        *size_received = size_communicated;
        break;

    case UNMIRROR_INDICES:
        size_communicated = *size_received;
        if (index_within != -1)
            size_communicated = size_local;

        mpi_send_data_from_list_rank(list_last_ranks_within_butterfly,
                                     list_rank_excess_butterfly,
                                     number_ranks_to_send, indices_local,
                                     size_communicated * sizeof(int),
                                     indices_received, tag + 4, worldcomm);
        break;

    case UNMIRROR_DATA:
        // Case unmirror : the last processes within the 2^k butterfly processes
        // send their data to the excess processes over the 2^k processes after
        // the Butterfly scheme
        size_communicated = *size_received;
        if (index_within != -1)
            size_communicated = size_local;

        mpi_send_data_from_list_rank(list_last_ranks_within_butterfly,
                                     list_rank_excess_butterfly,
                                     number_ranks_to_send, values_local,
                                     size_communicated * sizeof(double),
                                     values_received, tag + 5, worldcomm);
        break;
    }
    free(list_rank_excess_butterfly);
    free(list_last_ranks_within_butterfly);
    return 0;
}

int construct_butterfly_struct(Butterfly_struct *Butterfly_obj, int *indices_in,
                               int count_in, int *indices_out, int count_out,
                               int flag_classic_or_reshuffle_butterfly,
                               MPI_Comm worldcomm) {
    /* Initialize the butterfly communication, assuming the mirroring of the
     * data have been already done */
    int rank, size;
    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &size);

    Butterfly_obj->classic_or_reshuffle_butterfly =
        flag_classic_or_reshuffle_butterfly;
    // 0 if classic butterfly, e.g. if pixel distributions are the same
    // before/after ; 1 if reshuffle butterfly

    // Construct the butterfly communicator
    MPI_Comm comm_butterfly;
    int number_steps = log_2(size);

    mpi_create_subset(pow_2(number_steps), worldcomm, &comm_butterfly);

    Butterfly_obj->steps = number_steps;

    if (rank < pow(2, number_steps)) {
        Butterfly_obj->comm_butterfly = comm_butterfly;
        Butterfly_obj->S = (int **)malloc(Butterfly_obj->steps * sizeof(int *));
        // allocate sending maps tab
        Butterfly_obj->R = (int **)malloc(Butterfly_obj->steps * sizeof(int *));
        // allocate receiving maps tab
        Butterfly_obj->nS = (int *)malloc(Butterfly_obj->steps * sizeof(int));
        // allocate sending map sizes tab
        Butterfly_obj->nR = (int *)malloc(Butterfly_obj->steps * sizeof(int));
        // allocate receiving map size tab

        int nb_butterfly_ranks = pow_2(number_steps);
        int i;
        switch (flag_classic_or_reshuffle_butterfly) {
        case 0: // Classic butterfly
            butterfly_init(indices_in, count_in, Butterfly_obj->R,
                           Butterfly_obj->nR, Butterfly_obj->S,
                           Butterfly_obj->nS, &(Butterfly_obj->com_indices),
                           &(Butterfly_obj->com_count), Butterfly_obj->steps,
                           comm_butterfly);
            break;
        case 1: // Reshuffle butterfly
            butterfly_reshuffle_init(
                indices_in, count_in, indices_out, count_out, Butterfly_obj->R,
                Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS,
                &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count),
                Butterfly_obj->steps, comm_butterfly);
            break;
        }
    }
    return 0;
}

int prepare_butterfly_communication(
    int *indices_in, int count_in, int *indices_out, int count_out,
    int flag_classic_or_reshuffle_butterfly,
    Butterfly_superstruct *Butterfly_superstruct_obj, MPI_Comm worldcomm) {
    /* int *indices_in, int count_in must be the indices and size of what will
       be given to the Butterfly communication, while int *indices_out, int
       count_out will be the output indices and size

        Both indices_in and indices_out are expected to be ordered
    */
    int rank, nprocs;

    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &nprocs);

    int check_monotony_in = monotony(indices_in, count_in);
    int check_monotony_out = monotony(indices_out, count_out);

    if (check_monotony_in || check_monotony_out) {
        printf("ERROR : In prepare_butterfly_communication indices_in or "
               "indices_out is not monotonous : \n 1 for non-monotonous, 0 for "
               "monotonous : indices_in %d indices_out %d\n",
               check_monotony_in, check_monotony_out);
        return 1;
    }

    Butterfly_struct *Butterfly_forward_obj =
        &(Butterfly_superstruct_obj->Butterfly_obj);
    Butterfly_struct_supplement *Butterfly_mirror_supp =
        &(Butterfly_superstruct_obj->Butterfly_mirror_supp);
    Butterfly_struct_supplement *Butterfly_unmirror_supp =
        &(Butterfly_superstruct_obj->Butterfly_unmirror_supp);

    // First, exchanges sizes and indices for mirroring and unmirroring steps
    int new_count_in, size_received_mirror, size_to_send_unmirror,
        new_size_mirror, new_size_unmirror;
    int *indices_received_mirror, *indices_to_send_unmirror;
    int *ordered_indices_mirror, *ordered_indices_to_unmirror;

    size_received_mirror = 0;
    // Mirror step
    mirror_butterfly(NULL, NULL, count_in, NULL, NULL, &size_received_mirror,
                     MIRROR_SIZE, worldcomm);
    if (size_received_mirror)
        indices_received_mirror =
            (int *)malloc(size_received_mirror * sizeof(int));
    mirror_butterfly(NULL, indices_in, count_in, NULL, indices_received_mirror,
                     &size_received_mirror, MIRROR_INDICES, worldcomm);

    size_to_send_unmirror = 0;
    // Unmirror step
    mirror_butterfly(NULL, NULL, count_out, NULL, NULL, &size_to_send_unmirror,
                     MIRROR_SIZE, worldcomm);
    if (size_to_send_unmirror)
        indices_to_send_unmirror =
            (int *)malloc(size_to_send_unmirror * sizeof(int));
    mirror_butterfly(NULL, indices_out, count_out, NULL,
                     indices_to_send_unmirror, &size_to_send_unmirror,
                     MIRROR_INDICES, worldcomm);

    // Construct the new indices tab for mirorring
    new_size_mirror = modified_card_or(
        indices_in, count_in, indices_received_mirror, size_received_mirror);
    if (new_size_mirror)
        ordered_indices_mirror = (int *)malloc(new_size_mirror * sizeof(int));
    modified_set_or(indices_in, count_in, indices_received_mirror,
                    size_received_mirror, ordered_indices_mirror);

    // Construct the new indices tab for mirorring
    new_size_unmirror =
        modified_card_or(indices_out, count_out, indices_to_send_unmirror,
                         size_to_send_unmirror);
    if (new_size_unmirror)
        ordered_indices_to_unmirror =
            (int *)malloc(new_size_unmirror * sizeof(int));
    modified_set_or(indices_out, count_out, indices_to_send_unmirror,
                    size_to_send_unmirror, ordered_indices_to_unmirror);

    // Preparing mirroring
    Butterfly_mirror_supp->indices_mirror = indices_received_mirror;
    // List of indices received in the mirroring step by the local MPI process
    // -> Will be used afterwards to resend back the relevant indices when
    // unmirroring
    Butterfly_mirror_supp->size_from_mirror = size_received_mirror;
    // Size of the indices obtained from mirroring

    Butterfly_mirror_supp->ordered_indices = ordered_indices_mirror;
    // Ordered indices used for butterfly scheme, obtained from mirroring
    Butterfly_mirror_supp->new_size_local = new_size_mirror;
    // New size after mirroring when the redundant indices have been deleted

    // Preparing unmirroring
    Butterfly_unmirror_supp->indices_mirror = indices_to_send_unmirror;
    // List of indices received in the mirroring step by the local MPI process
    // -> Will be used afterwards to resend back the relevant indices when
    // unmirroring
    Butterfly_unmirror_supp->size_from_mirror = size_to_send_unmirror;
    // Size of the indices obtained from mirroring

    Butterfly_unmirror_supp->ordered_indices = ordered_indices_to_unmirror;
    // List of indices received in the mirroring step by the local MPI process
    // -> Will be used afterwards to resend back the relevant indices when
    // unmirroring
    Butterfly_unmirror_supp->new_size_local = new_size_unmirror;
    // Size of the indices obtained from mirroring

    if (new_size_mirror == 0) {
        new_size_mirror = count_in;
    }
    if (new_size_unmirror == 0) {
        new_size_unmirror = count_out;
    }

    // Construct the Butterfly_struct
    construct_butterfly_struct(Butterfly_forward_obj, ordered_indices_mirror,
                               new_size_mirror, ordered_indices_to_unmirror,
                               new_size_unmirror,
                               flag_classic_or_reshuffle_butterfly, worldcomm);

    return 0;
}

int perform_butterfly_communication(
    double *values_to_communicate, int *indices_in, int count_in,
    double *values_out, int *indices_out, int count_out,
    Butterfly_superstruct *Butterfly_superstruct_obj, MPI_Comm worldcomm) {
    /* #### Perform the butterfly communication ###
        ATTENTION : WILL ASSUME ALL MAP PIXEL VALUES ARE EQUAL INDEPENDANT OF
       THE PROCESS WHICH WILL HOLD IT */

    int check_monotony_in = monotony(indices_in, count_in);
    int check_monotony_out = monotony(indices_out, count_out);

    if (check_monotony_in || check_monotony_out) {
        printf("ERROR : In prepare_butterfly_communication indices_in or "
               "indices_out is not monotonous : \n 1 for non-monotonous, 0 for "
               "monotonous : indices_in %d indices_out %d \n",
               check_monotony_in, check_monotony_out);
        fflush(stdout);
        return 1;
    }

    Butterfly_struct *Butterfly_obj =
        &(Butterfly_superstruct_obj->Butterfly_obj);
    Butterfly_struct_supplement *Butterfly_mirror_supp =
        &(Butterfly_superstruct_obj->Butterfly_mirror_supp);
    Butterfly_struct_supplement *Butterfly_unmirror_supp =
        &(Butterfly_superstruct_obj->Butterfly_unmirror_supp);
    int k, i;
    int nSmax, nRmax;
    double *com_val, *values_received_butterfly;
    double *values_to_unmirror;
    double *values_butterfly;
    double *values_received;

    int rank, size;
    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &size);
    int number_steps = log_2(size);

    MPI_Comm comm_butterfly = Butterfly_obj->comm_butterfly;

    // Apply mirroring -> send data from excess processes to the processes which
    // will be used by the Butterfly scheme
    if (Butterfly_mirror_supp->size_from_mirror)
        values_received = (double *)malloc(
            Butterfly_mirror_supp->size_from_mirror * sizeof(double));

    mirror_butterfly(values_to_communicate, NULL, count_in, values_received,
                     NULL, &(Butterfly_mirror_supp->size_from_mirror),
                     MIRROR_DATA, worldcomm);
    if (rank < pow(2, number_steps)) {
        if (Butterfly_mirror_supp->new_size_local)
            values_butterfly = (double *)malloc(
                Butterfly_mirror_supp->new_size_local * sizeof(double));

        if (Butterfly_mirror_supp->size_from_mirror) {
            m2m(values_received, Butterfly_mirror_supp->indices_mirror,
                Butterfly_mirror_supp->size_from_mirror, values_butterfly,
                Butterfly_mirror_supp->ordered_indices,
                Butterfly_mirror_supp->new_size_local);
            free(values_received);
        }

        m2m(values_to_communicate, indices_in, count_in, values_butterfly,
            Butterfly_mirror_supp->ordered_indices,
            Butterfly_mirror_supp->new_size_local);
        // Compute max communication buffer size
        nRmax = Butterfly_obj->nR[0];
        nSmax = Butterfly_obj->nS[0];
        for (k = 1; k < Butterfly_obj->steps; k++) {
            if (Butterfly_obj->nR[k] > nRmax)
                nRmax = Butterfly_obj->nR[k];
            if (Butterfly_obj->nS[k] > nSmax)
                nSmax = Butterfly_obj->nS[k];
        }

        /* Copy value */
        com_val = (double *)calloc(Butterfly_obj->com_count, sizeof(double));

        m2m(values_butterfly, Butterfly_mirror_supp->ordered_indices,
            Butterfly_mirror_supp->new_size_local, com_val,
            Butterfly_obj->com_indices, Butterfly_obj->com_count);
        switch (Butterfly_obj->classic_or_reshuffle_butterfly) {
        case 0: // Classic butterfly
            butterfly_reduce(Butterfly_obj->R, Butterfly_obj->nR, nRmax,
                             Butterfly_obj->S, Butterfly_obj->nS, nSmax,
                             com_val, Butterfly_obj->steps,
                             Butterfly_obj->comm_butterfly);
            break;

        case 1: // Reshuffle butterfly
            butterfly_reshuffle(Butterfly_obj->R, Butterfly_obj->nR, nRmax,
                                Butterfly_obj->S, Butterfly_obj->nS, nSmax,
                                com_val, Butterfly_obj->steps,
                                Butterfly_obj->comm_butterfly);
            break;
        }

        if (Butterfly_unmirror_supp->new_size_local)
            values_received_butterfly = (double *)calloc(
                Butterfly_unmirror_supp->new_size_local, sizeof(double));

        m2m(values_butterfly, Butterfly_mirror_supp->ordered_indices,
            Butterfly_mirror_supp->new_size_local, values_received_butterfly,
            Butterfly_unmirror_supp->ordered_indices,
            Butterfly_unmirror_supp->new_size_local);
        m2m(com_val, Butterfly_obj->com_indices, Butterfly_obj->com_count,
            values_received_butterfly, Butterfly_unmirror_supp->ordered_indices,
            Butterfly_unmirror_supp->new_size_local);
        m2m(values_received_butterfly, Butterfly_unmirror_supp->ordered_indices,
            Butterfly_unmirror_supp->new_size_local, values_out, indices_out,
            count_out);

        if (Butterfly_unmirror_supp->size_from_mirror)
            values_to_unmirror = (double *)calloc(
                Butterfly_unmirror_supp->size_from_mirror, sizeof(double));
        m2m(values_received_butterfly, Butterfly_unmirror_supp->ordered_indices,
            Butterfly_unmirror_supp->new_size_local, values_to_unmirror,
            Butterfly_unmirror_supp->indices_mirror,
            Butterfly_unmirror_supp->size_from_mirror);
        if (Butterfly_mirror_supp->new_size_local)
            free(values_butterfly);
        if (Butterfly_unmirror_supp->new_size_local)
            free(values_received_butterfly);
        free(com_val);
    }
    mirror_butterfly(values_to_unmirror, NULL,
                     Butterfly_unmirror_supp->size_from_mirror, values_out,
                     NULL, &(count_out), UNMIRROR_DATA, worldcomm);
    if (Butterfly_unmirror_supp->size_from_mirror) {
        free(values_to_unmirror);
    }
    return 0;
}

int free_butterfly_struct(Butterfly_struct *Butterfly_obj, int rank) {

    if (rank < pow(2, Butterfly_obj->steps)) {
        free(Butterfly_obj->R);
        free(Butterfly_obj->nR);
        free(Butterfly_obj->S);
        free(Butterfly_obj->nS);
    }
    return 0;
}

int free_butterfly_supplement(Butterfly_struct_supplement *Butterfly_supp,
                              int rank) {
    if (Butterfly_supp->size_from_mirror)
        free(Butterfly_supp->indices_mirror);

    if (Butterfly_supp->new_size_local)
        free(Butterfly_supp->ordered_indices);

    return 0;
}

int free_butterfly_superstruct(Butterfly_superstruct *Butterfly_superstruct_obj,
                               int rank) {
    free_butterfly_supplement(
        &(Butterfly_superstruct_obj->Butterfly_mirror_supp), rank);
    free_butterfly_struct(&(Butterfly_superstruct_obj->Butterfly_obj), rank);
    free_butterfly_supplement(
        &(Butterfly_superstruct_obj->Butterfly_unmirror_supp), rank);

    return 0;
}

#endif
