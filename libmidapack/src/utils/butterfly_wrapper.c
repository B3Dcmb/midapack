#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "sys/types.h"
#include "sys/stat.h"
#include "sys/param.h"
#ifdef W_MPI
#include <mpi.h>
#include "mpi_tools.h"

#include "./mapmat/alm.h"
#include "../mapmat/csort.h"
#include "../mapmat/butterfly.h"
#include "../mapmat/bitop.h"
#include "butterfly_new.h"
// #include "midapack.h"

int mirror_butterfly(double *values_local, int *indices_local, int size_local, double *values_received, int *indices_received, int *size_received, int flag_mirror_unmirror_size_indices_data, MPI_Comm worldcomm)
{
    /* Sends data of the excess processes, which are over the 2^k processes which will be used for the Butterfly scheme, back to the 2^k processes. 
        The "mirror" step corresponds to sending the information of the excess processes (over 2^k) to the Butterfly processes, 
        while the "unmirror" step corresponds to sending back the information of the Butterfly processes to the excess processes.
    
        The flag flag_mirror_unmirror_size_indices_data should be either MIRROR_SIZE, MIRROR_INDICES, MIRROR_DATA, UNMIRROR_SIZE, UNMIRROR_INDICES or UNMIRROR_DATA
        to respectively mirror/unmirror the size of the data, indices of the data, or the data itself 
    
        The way this proceed is by taking the data of the N-2^k processes and send them respectively to the last N-2^k processes which will be used for the Butterfly scheme. 
        This means the data of the process ranked 2^k will be send to the process ranked 2^k-1, and the data of the process ranked 2^k+10 will be send to the process ranked 2^k-11 
        This function should be followed by reorder_indices to rearrange the indices and associated values, to make sure
        the received indices are order monotonously and the redundant indices and values are treated correctly */
    
    int rank, size;
    MPI_Comm_size(worldcomm, &size);
    MPI_Comm_rank(worldcomm, &rank);
    int i;
    int number_steps_Butterfly = log_2(size);
    int nb_rank_Butterfly = pow_2(number_steps_Butterfly);
    int *list_rank_excess_butterfly, *list_last_ranks_within_butterfly;

    int number_ranks_to_send = size - nb_rank_Butterfly;
    int tag=100;
    // printf("r %d ~~~~ nb_steps %d ; nb_ranks %d ; number_ranks_to_send %d \n", rank, number_steps_Butterfly, nb_rank_Butterfly, number_ranks_to_send);
    // fflush(stdout);
    if (number_ranks_to_send == 0)
        return 0;

    if (number_ranks_to_send){
        list_rank_excess_butterfly = (int *)malloc(number_ranks_to_send*sizeof(int));
        list_last_ranks_within_butterfly = (int *)malloc(number_ranks_to_send*sizeof(int));
    }
    printf("###COM-in %d ~~~~ TEST0 !!!! \n", rank); fflush(stdout);
    for (i=0; i<number_ranks_to_send; i++){
        list_rank_excess_butterfly[i] = nb_rank_Butterfly + i;
        list_last_ranks_within_butterfly[i] = nb_rank_Butterfly - i - 1;
    }
    printf("###COM-in %d ~~~~ TEST1 !!!! \n", rank); fflush(stdout);

    int index_excess = elem_in_list_elem(rank, list_rank_excess_butterfly, number_ranks_to_send);
    int index_within = elem_in_list_elem(rank, list_last_ranks_within_butterfly, number_ranks_to_send);

    printf("###COM-in %d ~~~~ index_within %d ; index_excess %d ; size_local %d ; size_received %d \n", rank, index_within, index_excess, size_local, *size_received); fflush(stdout);

    if ((index_excess==-1 && index_within==-1) || number_ranks_to_send==0){
        printf("###COM-out %d ~~~~ index_within %d ; index_excess %d ; size_local %d ; size_received %d \n", rank, index_within, index_excess, size_local, *size_received); fflush(stdout);
        if (number_ranks_to_send){
            free(list_rank_excess_butterfly);
            free(list_last_ranks_within_butterfly);
        }
        return 0;
    }
    // if (number_ranks_to_send>0)
    //     printf("r %d ~~~~ nb_steps %d ; nb_ranks %d ; nb_rank_to_send %d ; rank excess %d ; rank_within %d ; flag %d \n", rank, number_steps_Butterfly, nb_rank_Butterfly, number_ranks_to_send, list_rank_excess_butterfly[0], list_last_ranks_within_butterfly[0], flag_mirror_unmirror_size_indices_data);
    // fflush(stdout);

    // First, communicate the sizes, determine which processes will receive/send data
    int size_communicated = 0;
    switch(flag_mirror_unmirror_size_indices_data)
    {
        case MIRROR_SIZE: // Case 0 : MIRROR - send only size and indices to be exchanged, to have the excess processes over 2^k to send their data to the last 2^k processes before the Butterfly scheme
            mpi_send_data_from_list_rank(list_rank_excess_butterfly, list_last_ranks_within_butterfly, number_ranks_to_send, &size_local, 1*sizeof(int), &size_communicated, tag, worldcomm);
            *size_received = size_communicated;
            // printf("r %d ~~~~ size received : %d, %d, %d \n", rank, size_local, *size_received, size_communicated);
            break;

        case MIRROR_INDICES:
            size_communicated = *size_received;
            if (index_excess != -1)
                size_communicated = size_local;
            // printf("###COM %d ~~~~ size %d ; indice send : %d \n", rank, size_communicated, indices_local[0]);
            mpi_send_data_from_list_rank(list_rank_excess_butterfly, list_last_ranks_within_butterfly, number_ranks_to_send, indices_local, size_communicated*sizeof(int), indices_received, tag+1, worldcomm);
            if (size_communicated>0)
                // printf("###COM %d ~~~~ size %d ; indice received : %d \n", rank, size_communicated, indices_received[size_communicated-1]);
            break;

        case MIRROR_DATA:
            size_communicated = *size_received;
            if (index_excess != -1)
                size_communicated = size_local;
            mpi_send_data_from_list_rank(list_rank_excess_butterfly, list_last_ranks_within_butterfly, number_ranks_to_send, values_local, size_communicated*sizeof(double), values_received, tag+2, worldcomm);
            break;

        case UNMIRROR_SIZE:
            mpi_send_data_from_list_rank(list_last_ranks_within_butterfly, list_rank_excess_butterfly, number_ranks_to_send, &size_local, 1*sizeof(int), &size_communicated, tag+3, worldcomm);
            *size_received = size_communicated;
            // printf("UMr %d ~~~~ size received : %d, %d, %d \n", rank, size_local, *size_received, size_communicated);
            break;

        case UNMIRROR_INDICES:
            size_communicated = *size_received;
            if (index_within != -1)
                size_communicated = size_local;
            // if (size_local>0)
                // printf("###COMU %d ~~~~ size %d ; indice send : %d \n", rank, size_communicated, indices_local[0]);
            mpi_send_data_from_list_rank(list_last_ranks_within_butterfly, list_rank_excess_butterfly, number_ranks_to_send, indices_local, size_communicated*sizeof(int), indices_received, tag+4, worldcomm);
            if (size_communicated>0)
                // printf("###COMU %d ~~~~ size %d ; indice received : %d \n", rank, size_communicated, indices_received[size_communicated-1]);
            break;
        
        case UNMIRROR_DATA: // Case unmirror : the last processes within the 2^k butterfly processes send their data to the excess processes over the 2^k processes after the Butterfly scheme
            size_communicated = *size_received;
            if (index_within != -1)
                size_communicated = size_local;
            if (size_communicated>0){
                printf("###COMU %d ~~~~ index_within %d ; index_excess %d ; size %d \n", rank, index_within, index_excess, size_communicated); fflush(stdout);
                if (index_within != -1){
                    printf("###COMUn %d ~~~~ size %d ; values_local %f -", rank, size_communicated, values_local[0]); fflush(stdout);
                    for (i=1; i<size_communicated; i++)
                        printf("- %f -", values_local[i]);
                    printf("\n"); fflush(stdout);
                }
            }
            mpi_send_data_from_list_rank(list_last_ranks_within_butterfly, list_rank_excess_butterfly, number_ranks_to_send, values_local, size_communicated*sizeof(double), values_received, tag+5, worldcomm);
            if (index_excess != -1){
                printf("###COMUn %d ~~~~ size %d ; values_local %f -", rank, size_communicated, values_received[0]);
                for (i=1; i<size_communicated; i++)
                    printf("- %f -", values_received[i]);
                printf("\n"); fflush(stdout);
            }
            // if (size_communicated>0)
            //     printf("###COMU %d ~~~~ size %d ; value received : %f \n", rank, size_communicated, values_received[size_communicated-1]);
            break;
    }

    free(list_rank_excess_butterfly);
    free(list_last_ranks_within_butterfly);
    return 0;
}


int construct_butterfly_struct(Butterfly_struct *Butterfly_obj, int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, MPI_Comm worldcomm)
{
    /* Initialize the butterfly communication, assuming the mirroring of the data have been already done */
    int rank, size;
    MPI_Comm_rank(worldcomm, &rank);
    MPI_Comm_size(worldcomm, &size);

    Butterfly_obj->classic_or_reshuffle_butterfly = flag_classic_or_reshuffle_butterfly;
    // 0 if classic butterfly, e.g. if pixel distributions are the same before/after ; 1 if reshuffle butterfly
    
    // Construct the butterfly communicator
    MPI_Comm comm_butterfly;
    int number_steps = log_2(size);
    // printf("%d ~~~~ TEST-1 \n", rank); fflush(stdout);
    mpi_create_subset(pow_2(number_steps), worldcomm, &comm_butterfly);

    Butterfly_obj->steps = number_steps;
    // printf("%d ~~~~ TEST0 \n", rank); fflush(stdout);

    if (rank < pow(2,number_steps)){
        Butterfly_obj->comm_butterfly = comm_butterfly;
        Butterfly_obj->S = (int **)malloc(Butterfly_obj->steps * sizeof(int *)); // allocate sending maps tab
        Butterfly_obj->R = (int **)malloc(Butterfly_obj->steps * sizeof(int *)); // allocate receiving maps tab
        Butterfly_obj->nS = (int *)malloc(Butterfly_obj->steps * sizeof(int));   // allocate sending map sizes tab
        Butterfly_obj->nR = (int *)malloc(Butterfly_obj->steps * sizeof(int));   // allocate receiving map size tab
        // butterfly_init(A->lindices + (A->nnz) * (A->trash_pix), A->lcount - (A->nnz) * (A->trash_pix), A->R, A->nR, A->S, A->nS, &(A->com_indices), &(A->com_count), A->steps, comm);
        int nb_butterfly_ranks = pow_2(number_steps);
        int i;
        printf("%d [[]] Butterfly_init !! MODE %d ; nb_butterfly_ranks %d ; nb_steps %d ; count_in : %d ; count_out : %d \n", rank, flag_classic_or_reshuffle_butterfly, nb_butterfly_ranks, number_steps, count_in, count_out);
        printf("%d [[]] indices_in : %d -", rank, indices_in[0]);
        for (i=1; i<count_in; i++){
            printf("- %d -", indices_in[i]);
        }
        printf("\n"); fflush(stdout);
        printf("%d [[]] indices_out : %d -", rank, indices_out[0]);
        for (i=1; i<count_out; i++){
            printf("- %d -", indices_out[i]);
        }
        printf("\n"); fflush(stdout);
        // printf("%d ~~~~ TEST1 \n", rank); fflush(stdout);
        // int *indices_in_copy = (int *)malloc(count_in*sizeof(int));
        // memcpy(indices_in_copy, indices_in, count_in*sizeof(int));
        // int *indices_out_copy = (int *)malloc(count_out*sizeof(int));
        // memcpy(indices_out_copy, indices_out, count_out*sizeof(int));
        switch(flag_classic_or_reshuffle_butterfly)
        {
            case 0: // Classic butterfly
                // butterfly_reduce_init(indices_in, count_in, Butterfly_obj->R, Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS, &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count), Butterfly_obj->steps, comm_butterfly);
                butterfly_init(indices_in, count_in, Butterfly_obj->R, Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS, &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count), Butterfly_obj->steps, comm_butterfly);
                // butterfly_reduce_init(indices_out, count_out, Butterfly_obj->R, Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS, &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count), Butterfly_obj->steps, comm_butterfly);
                break;
            case 1: // Reshuffle butterfly

                butterfly_reshuffle_init(indices_in, count_in, indices_out, count_out, Butterfly_obj->R, Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS, &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count), Butterfly_obj->steps, comm_butterfly);
                break;
        }
        // free(indices_in_copy);
        // free(indices_out_copy);
    }
    printf("%d ~~~~ TEST2 \n", rank); fflush(stdout);
    return 0;
}

int prepare_butterfly_communication(int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, Butterfly_superstruct *Butterfly_superstruct_obj, MPI_Comm worldcomm)
{
    /* int *indices_in, int count_in must be the indices and size of what will be given to the Butterfly communication,
        while int *indices_out, int count_out will be the output indices and size

        Both indices_in and indices_out are expected to be ordered
    */
    int rank, nprocs;

    MPI_Comm_rank( worldcomm, &rank);
    MPI_Comm_size( worldcomm, &nprocs);

    printf("%d ~~~ Prep : Starting preparing butterfly \n", rank); fflush(stdout);
    // Butterfly_struct *Butterfly_forward_obj = (Butterfly_struct *)malloc(1*sizeof(Butterfly_struct));
    // Butterfly_struct_supplement *Butterfly_mirror_supp = (Butterfly_struct_supplement *)malloc(1*sizeof(Butterfly_struct_supplement));
    // Butterfly_struct_supplement *Butterfly_unmirror_supp = (Butterfly_struct_supplement *)malloc(1*sizeof(Butterfly_struct_supplement));
    
    // Butterfly_superstruct_obj->Butterfly_obj = Butterfly_forward_obj;
    // Butterfly_superstruct_obj->Butterfly_mirror_supp = Butterfly_mirror_supp;
    // Butterfly_superstruct_obj->Butterfly_unmirror_supp = Butterfly_unmirror_supp;

    Butterfly_struct *Butterfly_forward_obj = &(Butterfly_superstruct_obj->Butterfly_obj);
    Butterfly_struct_supplement *Butterfly_mirror_supp = &(Butterfly_superstruct_obj->Butterfly_mirror_supp);
    Butterfly_struct_supplement *Butterfly_unmirror_supp = &(Butterfly_superstruct_obj->Butterfly_unmirror_supp);


    // First, exchanges sizes and indices for mirroring and unmirroring steps
    int new_count_in, size_received_mirror, size_to_send_unmirror, new_size_mirror, new_size_unmirror;
    int *indices_received_mirror, *indices_to_send_unmirror;
    int *ordered_indices_mirror, *ordered_indices_to_unmirror;

    printf("%d ~~~ Prep : Starting mirroring \n", rank); fflush(stdout);
    size_received_mirror = 0;
    // Mirror step
    mirror_butterfly(NULL, NULL, count_in, NULL, NULL, &size_received_mirror, MIRROR_SIZE, worldcomm);
    if (size_received_mirror)
        indices_received_mirror = (int *)malloc(size_received_mirror*sizeof(int));
    mirror_butterfly(NULL, indices_in, count_in, NULL, indices_received_mirror, &size_received_mirror, MIRROR_INDICES, worldcomm);
    
    printf("%d ~~~ Prep : Starting unmirroring \n", rank); fflush(stdout);
    size_to_send_unmirror = 0;
    // Unmirror step
    mirror_butterfly(NULL, NULL, count_out, NULL, NULL, &size_to_send_unmirror, MIRROR_SIZE, worldcomm);
    if (size_to_send_unmirror)
        indices_to_send_unmirror = (int *)malloc(size_to_send_unmirror*sizeof(int));
    mirror_butterfly(NULL, indices_out, count_out, NULL, indices_to_send_unmirror, &size_to_send_unmirror, MIRROR_INDICES, worldcomm);

    printf("%d ~~~ Prep : Constructing indices mirroring \n", rank); fflush(stdout);
    // Construct the new indices tab for mirorring
    // int *ordered_indices_mirror = (int *)malloc((count_in+size_received_mirror)*sizeof(int));
    // memcpy(ordered_indices_mirror, indices_in, count_in*sizeof(int));
    // memcpy(ordered_indices_mirror+count_in, indices_received_mirror, size_received_mirror*sizeof(int));
    new_size_mirror = modified_card_or(indices_in, count_in, indices_received_mirror, size_received_mirror);
    // int *ordered_indices_mirror_temp = (int *)malloc(new_size_mirror*sizeof(int));
    if (new_size_mirror)
        ordered_indices_mirror = (int *)malloc(new_size_mirror*sizeof(int));
    // modified_set_or(indices_in, count_in, indices_received_mirror, size_received_mirror, ordered_indices_mirror_temp);
    modified_set_or(indices_in, count_in, indices_received_mirror, size_received_mirror, ordered_indices_mirror);
    // new_size_mirror = modified_set_or(indices_in, count_in, indices_received_mirror, size_received_mirror, ordered_indices_mirror);
    // new_size_mirror = ssort(ordered_indices_mirror, count_in+size_received_mirror, 0); // Argument flag=0 to use quicksort to sort the indices    
    // ordered_indices_mirror = realloc(ordered_indices_mirror, new_size_mirror*sizeof(int));

    printf("%d ~~~ Prep : Constructing indices unmirroring \n", rank); fflush(stdout);
    // Construct the new indices tab for mirorring
    // int *ordered_indices_to_unmirror = (int *)malloc((count_out+size_to_send_unmirror)*sizeof(int));
    // memcpy(ordered_indices_to_unmirror, indices_out, count_out*sizeof(int));
    // memcpy(ordered_indices_to_unmirror+count_out, indices_to_send_unmirror, size_to_send_unmirror*sizeof(int));
    // new_size_unmirror = modified_card_or(indices_out, count_out, indices_to_send_unmirror, size_to_send_unmirror);
    // int *ordered_indices_to_unmirror = (int *)malloc(new_size_unmirror*sizeof(int));
    new_size_unmirror = modified_card_or(indices_out, count_out, indices_to_send_unmirror, size_to_send_unmirror);
    if (new_size_unmirror)
        ordered_indices_to_unmirror = (int *)malloc(new_size_unmirror*sizeof(int));
    modified_set_or(indices_out, count_out, indices_to_send_unmirror, size_to_send_unmirror, ordered_indices_to_unmirror);
    // new_size_unmirror = modified_set_or(indices_out, count_out, indices_to_send_unmirror, size_to_send_unmirror, ordered_indices_to_unmirror);
    

    // new_size_unmirror = ssort(ordered_indices_to_unmirror, count_out+size_to_send_unmirror, 0); // Argument flag=0 to use quicksort to sort the indices    
    // ordered_indices_to_unmirror = realloc(ordered_indices_to_unmirror, new_size_unmirror*sizeof(int));

    printf("%d ~~~ Prep : Attributing values computed -- new_size_mirror %d new_size_unmirror %d\n", rank, new_size_mirror, new_size_unmirror); fflush(stdout);
    // Preparing mirroring
    Butterfly_mirror_supp->indices_mirror = indices_received_mirror; // List of indices received in the mirroring step by the local MPI process  -> Will be used afterwards to resend back the relevant indices when unmirroring
    Butterfly_mirror_supp->size_from_mirror = size_received_mirror; // Size of the indices obtained from mirroring
    
    // memcpy(ordered_indices_mirror, ordered_indices_mirror_temp, new_size_mirror*sizeof(int));
    Butterfly_mirror_supp->ordered_indices = ordered_indices_mirror; // Ordered indices used for butterfly scheme, obtained from mirroring
    // free(ordered_indices_mirror_temp);
    Butterfly_mirror_supp->new_size_local = new_size_mirror; // New size after mirroring when the redundant indices have been deleted

    // Preparing unmirroring
    Butterfly_unmirror_supp->indices_mirror = indices_to_send_unmirror; // List of indices received in the mirroring step by the local MPI process  -> Will be used afterwards to resend back the relevant indices when unmirroring
    Butterfly_unmirror_supp->size_from_mirror = size_to_send_unmirror; // Size of the indices obtained from mirroring

    Butterfly_unmirror_supp->ordered_indices = ordered_indices_to_unmirror; // List of indices received in the mirroring step by the local MPI process  -> Will be used afterwards to resend back the relevant indices when unmirroring
    Butterfly_unmirror_supp->new_size_local = new_size_unmirror; // Size of the indices obtained from mirroring

    if (new_size_mirror == 0){
        new_size_mirror = count_in;
    }
    if (new_size_unmirror == 0){
        new_size_unmirror = count_out;
    }
    int k;
    if (rank < pow(2,log_2(nprocs))){
        
        printf("%d ~~~ Prep : Constructing butterfly_struct \n", rank);
        printf("%d ~~~ Post-mirroring check 0 !! %d \n", rank, Butterfly_mirror_supp->new_size_local);
        printf("%d ~~~ ordered_indices_mirror : %d", rank, Butterfly_mirror_supp->ordered_indices[0]);
        for (k=1; k<Butterfly_mirror_supp->new_size_local; k++){
            printf("- %d -",Butterfly_mirror_supp->ordered_indices[k]);
        }
        printf("\n");
        printf("%d ~~~ Post-unmirroring check 0 !! %d \n", rank, Butterfly_unmirror_supp->new_size_local);
        printf("%d ~~~ ordered_indices_to_unmirror : %d", rank, Butterfly_unmirror_supp->ordered_indices[0]);
        for (k=1; k<Butterfly_unmirror_supp->new_size_local; k++){
            printf("- %d -",Butterfly_unmirror_supp->ordered_indices[k]);
        }
        printf("\n"); fflush(stdout);
    }
    // Construct the Butterfly_struct
    construct_butterfly_struct(Butterfly_forward_obj, ordered_indices_mirror, new_size_mirror, ordered_indices_to_unmirror, new_size_unmirror, flag_classic_or_reshuffle_butterfly, worldcomm);
    
    printf("%d ~~~ Out of butterfly_struct construct 0 !! \n", rank); fflush(stdout);

    if (rank < pow(2,log_2(nprocs))){
        int *com_indices = Butterfly_forward_obj->com_indices;
        printf("%d ~~~ Post-constructing butterfly_struct check 0 !! %d \n", rank, Butterfly_forward_obj->com_count); fflush(stdout);
        if(Butterfly_forward_obj->com_count){
            printf("%d ~~~ com_indices : %d", rank, com_indices[0]);
            for (k=1; k<Butterfly_forward_obj->com_count; k++){
                printf("- %d -",com_indices[k]);
            }
            printf("\n"); fflush(stdout);
        }
    }
    printf("%d ~~~~ TEST3 \n", rank); fflush(stdout);
    return 0;
}



int perform_butterfly_communication(double *values_to_communicate, int *indices_in, int count_in, double *values_out, int *indices_out, int count_out, Butterfly_superstruct *Butterfly_superstruct_obj, MPI_Comm worldcomm)
{
    /* #### Perform the butterfly communication ### 
        ATTENTION : WILL ASSUME ALL MAP PIXEL VALUES ARE EQUAL INDEPENDANT OF THE PROCESS WHICH WILL HOLD IT */ 

    Butterfly_struct *Butterfly_obj = &(Butterfly_superstruct_obj->Butterfly_obj);
    Butterfly_struct_supplement *Butterfly_mirror_supp = &(Butterfly_superstruct_obj->Butterfly_mirror_supp);
    Butterfly_struct_supplement *Butterfly_unmirror_supp = &(Butterfly_superstruct_obj->Butterfly_unmirror_supp);
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

    // Apply mirroring -> send data from excess processes to the processes which will be used by the Butterfly scheme
    if (Butterfly_mirror_supp->size_from_mirror)
        values_received = (double *)malloc(Butterfly_mirror_supp->size_from_mirror*sizeof(double));
    // int copy_size_from_mirror = Butterfly_mirror_supp->size_from_mirror;

    printf("%d (((()))) Start performing butterfly, mirroring first ! count_in %d ; count_out %d ; size_from_mirror %d \n", rank, count_in, count_out, Butterfly_mirror_supp->size_from_mirror); fflush(stdout);
    mirror_butterfly(values_to_communicate, NULL, count_in, values_received, NULL, &(Butterfly_mirror_supp->size_from_mirror), MIRROR_DATA, worldcomm);

    if (rank < pow(2,number_steps)){
        if (Butterfly_mirror_supp->new_size_local)
            values_butterfly = (double *)malloc(Butterfly_mirror_supp->new_size_local*sizeof(double));
        if (Butterfly_mirror_supp->size_from_mirror){
            m2m(values_received, Butterfly_mirror_supp->indices_mirror, Butterfly_mirror_supp->size_from_mirror, values_butterfly, Butterfly_mirror_supp->ordered_indices, Butterfly_mirror_supp->new_size_local);
            free(values_received);
        }

        m2m(values_to_communicate, indices_in, count_in, values_butterfly, Butterfly_mirror_supp->ordered_indices, Butterfly_mirror_supp->new_size_local);
        printf("%d »»»»»»»» Post-mirroring check !! %d \n", rank, Butterfly_mirror_supp->new_size_local);
        printf("%d »»»»»»»» ordered_indices values_butterfly : %d %f", rank, Butterfly_mirror_supp->ordered_indices[0], values_butterfly[0]);
        for (k=1; k<Butterfly_mirror_supp->new_size_local; k++){
            printf("- %d %f -",Butterfly_mirror_supp->ordered_indices[k], values_butterfly[k]);
        }
        printf("\n"); fflush(stdout);

        // Compute max communication buffer size
        nRmax = Butterfly_obj->nR[0];
        nSmax = Butterfly_obj->nS[0];
        for (k = 1; k < Butterfly_obj->steps; k++){
            if (Butterfly_obj->nR[k] > nRmax)
                nRmax = Butterfly_obj->nR[k];
            if (Butterfly_obj->nS[k] > nSmax)
                nSmax = Butterfly_obj->nS[k];
        }

        /* Copy value */
        com_val = (double *)calloc(Butterfly_obj->com_count,sizeof(double));
        // for (k = 0; k < Butterfly_obj->com_count; k++)
        //     com_val[k] = 0.0;

        m2m(values_butterfly, Butterfly_mirror_supp->ordered_indices, Butterfly_mirror_supp->new_size_local, com_val, Butterfly_obj->com_indices, Butterfly_obj->com_count);
        
        printf("%d <<<<<< Set-up before Butterfly !! %d \n", rank, Butterfly_obj->com_count); fflush(stdout);
        if (Butterfly_obj->com_count){
            printf("%d <<<<<< com_indices com_val : %d %f", rank, Butterfly_obj->com_indices[0], com_val[0]);
            for (k=1; k<Butterfly_obj->com_count; k++){
                printf("- %d %f -",Butterfly_obj->com_indices[k], com_val[k]);
            }
            printf("\n"); fflush(stdout);
        }

        // if (rank == 1){
        //     if (Butterfly_mirror_supp->new_size_local){
        //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        //         printf("%d ««««««« Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local);
        //         printf("%d ««««««« order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        //         }
        //         printf("\n"); fflush(stdout);
        //     }
        // }
        switch(Butterfly_obj->classic_or_reshuffle_butterfly)
        {
            case 0: // Classic butterfly
                // modified_butterfly_reduce(Butterfly_obj->R, Butterfly_obj->nR, nRmax, Butterfly_obj->S, Butterfly_obj->nS, nSmax, com_val, Butterfly_obj->steps, Butterfly_obj->comm_butterfly);
                butterfly_reduce(Butterfly_obj->R, Butterfly_obj->nR, nRmax, Butterfly_obj->S, Butterfly_obj->nS, nSmax, com_val, Butterfly_obj->steps, Butterfly_obj->comm_butterfly);
                break;

            case 1: // Reshuffle butterfly
                butterfly_reshuffle(Butterfly_obj->R, Butterfly_obj->nR, nRmax, Butterfly_obj->S, Butterfly_obj->nS, nSmax, com_val, Butterfly_obj->steps, Butterfly_obj->comm_butterfly);
                break;
        }
        printf("%d <<<<<< Results from Butterfly !! %d \n", rank, Butterfly_obj->com_count);
        if (Butterfly_obj->com_count){
            printf("%d <<<<<< com_indices com_val : %d %f", rank, Butterfly_obj->com_indices[0], com_val[0]);
            for (k=1; k<Butterfly_obj->com_count; k++){
                printf("- %d %f -",Butterfly_obj->com_indices[k], com_val[k]);
            }
            printf("\n");
        }
        // if (rank == 1){
        //     if (Butterfly_mirror_supp->new_size_local){
        //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        //         printf("%d «««««««1 Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local);
        //         printf("%d «««««««1 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        //         }
        //         printf("\n"); fflush(stdout);
        //     }
        // }

        if (Butterfly_unmirror_supp->new_size_local)
            values_received_butterfly = (double *)malloc(Butterfly_unmirror_supp->new_size_local*sizeof(double));

        m2m(values_butterfly, Butterfly_mirror_supp->ordered_indices, Butterfly_mirror_supp->new_size_local, values_received_butterfly, Butterfly_unmirror_supp->ordered_indices, Butterfly_unmirror_supp->new_size_local);        
        m2m(com_val, Butterfly_obj->com_indices, Butterfly_obj->com_count, values_received_butterfly, Butterfly_unmirror_supp->ordered_indices, Butterfly_unmirror_supp->new_size_local);
        

        // if (rank == 1){
        //     if (Butterfly_mirror_supp->new_size_local){
        //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        //         printf("%d «««««««2 Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local);
        //         printf("%d «««««««2 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        //         }
        //         printf("\n"); fflush(stdout);
        //     }
        // }

        m2m(values_received_butterfly, Butterfly_unmirror_supp->ordered_indices, Butterfly_unmirror_supp->new_size_local, values_out, indices_out, count_out);

        // if (rank == 1){
        //     if (Butterfly_mirror_supp->new_size_local){
        //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        //         printf("%d «««««««2b Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local);
        //         printf("%d «««««««2b order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        //         }
        //         printf("\n"); fflush(stdout);
        //     }
        // }

        printf("%d >>>>>>> Results sent out !! %d \n", rank, count_out);
        printf("%d >>>>>>> indices_out values_out : %d %f", rank, indices_out[0], values_out[0]);
        for (k=1; k<count_out; k++){
            printf("- %d %f -",indices_out[k], values_out[k]);
        }
        printf("\n"); fflush(stdout);
        if (Butterfly_unmirror_supp->size_from_mirror)
            values_to_unmirror = (double *)calloc(Butterfly_unmirror_supp->size_from_mirror,sizeof(int));
        // double *values_to_unmirror_2 = (double *)malloc(Butterfly_unmirror_supp->size_from_mirror*sizeof(int));
        // if (rank == 1){
        //     if (Butterfly_mirror_supp->new_size_local){
        //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        //         // free(Butterfly_mirror_supp->ordered_indices);
        //         printf("%d «««««««2c Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d -- test size %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local, Butterfly_unmirror_supp->size_from_mirror);
        //         printf("%d «««««««2c order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        //         }
        //         printf("\n"); fflush(stdout);
        //         // free(Butterfly_mirror_supp->ordered_indices);
        //         if (Butterfly_unmirror_supp->size_from_mirror){
        //             printf("%d «««««««2c indices_mirror : %d -", rank, Butterfly_unmirror_supp->indices_mirror[0]);
        //             for (i=1; i<Butterfly_unmirror_supp->size_from_mirror; i++){
        //                 printf("- %d -", Butterfly_unmirror_supp->indices_mirror[i]);
        //             }
        //         }
        //         printf("\n"); fflush(stdout);
        //     }
        // }
        // m2m(values_received_butterfly, Butterfly_unmirror_supp->ordered_indices, Butterfly_unmirror_supp->new_size_local, values_out, indices_out, count_out);
        if (Butterfly_unmirror_supp->new_size_local){
                // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
                // free(Butterfly_mirror_supp->ordered_indices);
                printf("%d «««««««2ccc Butterfly_unmirror_supp !! nb_steps %d ; new_size_local : %d -- test size_from_mirror %d \n", rank, number_steps, Butterfly_unmirror_supp->new_size_local, Butterfly_unmirror_supp->size_from_mirror);
                printf("%d «««««««2cccc order_indices mirror : %d %f -", rank, Butterfly_unmirror_supp->ordered_indices[0], values_received_butterfly[0]);
                for (i=1; i<Butterfly_unmirror_supp->new_size_local; i++){
                    printf("- %d %f -", Butterfly_unmirror_supp->ordered_indices[i], values_received_butterfly[i]);
                }
                printf("\n"); fflush(stdout);
                // free(Butterfly_mirror_supp->ordered_indices);
                if (Butterfly_unmirror_supp->size_from_mirror){
                    printf("%d «««««««2cccc indices_mirror : %d %f -", rank, Butterfly_unmirror_supp->indices_mirror[0], values_to_unmirror[0]);
                    for (i=1; i<Butterfly_unmirror_supp->size_from_mirror; i++){
                        printf("- %d %f -", Butterfly_unmirror_supp->indices_mirror[i], values_to_unmirror[i]);
                    }
                }
                printf("\n"); fflush(stdout);
            }
        m2m(values_received_butterfly, Butterfly_unmirror_supp->ordered_indices, Butterfly_unmirror_supp->new_size_local, values_to_unmirror, Butterfly_unmirror_supp->indices_mirror, Butterfly_unmirror_supp->size_from_mirror);
        if (Butterfly_unmirror_supp->new_size_local){
            // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
            // free(Butterfly_mirror_supp->ordered_indices);
            printf("%d «««««««2d00000 Butterfly_unmirror_supp !! nb_steps %d ; new_size_local : %d -- test size_from_mirror %d \n", rank, number_steps, Butterfly_unmirror_supp->new_size_local, Butterfly_unmirror_supp->size_from_mirror);
            printf("%d «««««««2d00000 order_indices mirror : %d %f -", rank, Butterfly_unmirror_supp->ordered_indices[0], values_received_butterfly[0]);
            for (i=1; i<Butterfly_unmirror_supp->new_size_local; i++){
                printf("- %d %f -", Butterfly_unmirror_supp->ordered_indices[i], values_received_butterfly[i]);
            }
            printf("\n"); fflush(stdout);
            // free(Butterfly_mirror_supp->ordered_indices);
            if (Butterfly_unmirror_supp->size_from_mirror){
                printf("%d «««««««2d00000 indices_mirror : %d %f -", rank, Butterfly_unmirror_supp->indices_mirror[0], values_to_unmirror[0]);
                for (i=1; i<Butterfly_unmirror_supp->size_from_mirror; i++){
                    printf("- %d %f -", Butterfly_unmirror_supp->indices_mirror[i], values_to_unmirror[i]);
                }
            }
            printf("\n"); fflush(stdout);
        }
        // if (rank == 1){
        //     if (Butterfly_mirror_supp->new_size_local){
        //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        //         printf("%d «««««««2d Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d -- test size %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local, Butterfly_unmirror_supp->size_from_mirror);
        //         printf("%d «««««««2d order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        //         }
        //         printf("\n"); fflush(stdout);
        //         // printf("%d «««««««2d indices_mirror : %d -", rank, Butterfly_mirror_supp->indices_mirror[0]);
        //         // for (i=1; i<Butterfly_mirror_supp->size_from_mirror; i++){
        //         //     printf("- %d -", Butterfly_mirror_supp->indices_mirror[i]);
        //         // }
        //         // printf("\n"); fflush(stdout);
        //     }
        // }
        // memcpy(values_to_unmirror, values_to_unmirror_2, Butterfly_unmirror_supp->size_from_mirror*sizeof(double));
        // mirror_butterfly(values_to_unmirror, NULL, Butterfly_unmirror_supp->size_from_mirror, NULL, NULL, &(Butterfly_unmirror_supp->size_from_mirror), UNMIRROR_DATA, worldcomm);
        // printf("%d Test perform #### Unmirror_data : %f \n", rank, values_to_unmirror[0]); fflush(stdout);

        // if (rank == 1){
        //     if (Butterfly_mirror_supp->new_size_local){
        //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        //         printf("%d «««««««3 Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local);
        //         printf("%d «««««««3 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        //         }
        //         printf("\n"); fflush(stdout);
        //     }
        // }


        // if (rank == 1){
        //     if (Butterfly_mirror_supp->new_size_local){
        //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        //         printf("%d «««««««4 Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local);
        //         printf("%d «««««««4 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        //         }
        //         printf("\n"); fflush(stdout);
        //     }
        // }

        if (Butterfly_mirror_supp->new_size_local)
            free(values_butterfly);
        if (Butterfly_unmirror_supp->new_size_local)
            free(values_received_butterfly);
        free(com_val);
    }

    // if (rank >= pow(2,number_steps)){
    //     // mirror_butterfly(NULL, NULL, 0, values_out, NULL, &(count_out), UNMIRROR_DATA, worldcomm);
    //     printf("%d Test perform #### Unmirror_data : %f \n", rank, values_out[0]);
    // }

    if (Butterfly_unmirror_supp->new_size_local){
        // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        // free(Butterfly_mirror_supp->ordered_indices);
        printf("%d «««««««2efffff Butterfly_unmirror_supp !! nb_steps %d ; new_size_local : %d -- test size_from_mirror %d \n", rank, number_steps, Butterfly_unmirror_supp->new_size_local, Butterfly_unmirror_supp->size_from_mirror); fflush(stdout);
        printf("%d «««««««2efffff order_indices mirror : %d %f -", rank, Butterfly_unmirror_supp->ordered_indices[0], values_out[0]); fflush(stdout);
        for (i=1; i<Butterfly_unmirror_supp->new_size_local; i++){
            printf("- %d %f -", Butterfly_unmirror_supp->ordered_indices[i], values_out[i]);
        }
        printf("\n"); fflush(stdout);
        // free(Butterfly_mirror_supp->ordered_indices);
    }
    if (Butterfly_unmirror_supp->size_from_mirror){
        printf("%d «««««««2efffff size %d indices_mirror values_to_unmirror : %d %f -", rank, Butterfly_unmirror_supp->size_from_mirror, Butterfly_unmirror_supp->indices_mirror[0], values_to_unmirror[0]); fflush(stdout);
        for (i=1; i<Butterfly_unmirror_supp->size_from_mirror; i++){
            printf("- %d %f -", Butterfly_unmirror_supp->indices_mirror[i], values_to_unmirror[i]);
        }
        printf("\n"); fflush(stdout);
    }
    
    // MPI_WAIT(MPI_Request *request, MPI_STATUS_IGNORE);

    printf("%d «««««««2g0 prepare unmirroring : %d \n", rank, count_out); fflush(stdout);
    mirror_butterfly(values_to_unmirror, NULL, Butterfly_unmirror_supp->size_from_mirror, values_out, NULL, &(count_out), UNMIRROR_DATA, worldcomm);
    printf("%d «««««««2g01 finish unmirroring : %d \n", rank, count_out); fflush(stdout);
    if (count_out){
        printf("%d «««««««2ggggg indices_mirror : %d %f -", rank, indices_out[0], values_out[0]); fflush(stdout);
        for (i=1; i<count_out; i++){
            printf("- %d %f -", indices_out[i], values_out[i]);
        }
        printf("\n"); fflush(stdout);
    }

    if (Butterfly_unmirror_supp->size_from_mirror){
        free(values_to_unmirror);
        
    }
    
    // if (rank == 1){
    //     if (Butterfly_mirror_supp->new_size_local){
    //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
    //         printf("%d «««««««5 Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local);
    //         printf("%d «««««««5 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
    //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
    //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
    //         }
    //         printf("\n"); fflush(stdout);
    //     }
    // }
    // // free(values_received);
    // if (rank == 1){
    //     if (Butterfly_mirror_supp->new_size_local){
    //         // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
    //         printf("%d «««««««6 Butterfly_supplement tests mirror !! nb_steps %d ; new_size_local : %d \n", rank, number_steps, Butterfly_mirror_supp->new_size_local); fflush(stdout);
    //         printf("%d «««««««6 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
    //         for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
    //             printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
    //         }
    //         printf("\n"); fflush(stdout);
    //     }
    // }
    return 0;
}


int free_butterfly_struct(Butterfly_struct *Butterfly_obj, int rank)
{

    printf("%d fffff Free 0a \n", rank); fflush(stdout);
    if (rank < pow(2, Butterfly_obj->steps)){
        free(Butterfly_obj->R);
        free(Butterfly_obj->nR);
        free(Butterfly_obj->S);
        free(Butterfly_obj->nS);
        // free(Butterfly_obj->com_indices);
    }
    printf("%d fffff Free 0b \n", rank); fflush(stdout);
    // free(Butterfly_obj);
    return 0;
}

int free_butterfly_supplement(Butterfly_struct_supplement *Butterfly_supp, int rank)
{
    printf("%d fffff Free 1a - %d \n", rank, Butterfly_supp->size_from_mirror); fflush(stdout);
    if (Butterfly_supp->size_from_mirror)
        free(Butterfly_supp->indices_mirror);
    printf("%d fffff Free 1b - %d \n", rank, Butterfly_supp->new_size_local); fflush(stdout);
    if (Butterfly_supp->new_size_local)
        free(Butterfly_supp->ordered_indices);
    printf("%d fffff Free 1c \n", rank); fflush(stdout);
    // free(Butterfly_supp);
    return 0;
}


int free_butterfly_superstruct(Butterfly_superstruct *Butterfly_superstruct_obj, int rank)
{
    printf("%d fffff Free 1 \n", rank); fflush(stdout);
    free_butterfly_supplement(&(Butterfly_superstruct_obj->Butterfly_mirror_supp), rank);
    printf("%d fffff Free 0 \n", rank); fflush(stdout);
    free_butterfly_struct(&(Butterfly_superstruct_obj->Butterfly_obj), rank);
    printf("%d fffff Free 2 \n", rank); fflush(stdout);
    free_butterfly_supplement(&(Butterfly_superstruct_obj->Butterfly_unmirror_supp), rank);
    printf("%d fffff Free 3 \n", rank); fflush(stdout);
    // free(Butterfly_superstruct_obj);
    return 0;
}

#endif
