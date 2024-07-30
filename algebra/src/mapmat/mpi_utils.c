#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mapmat/mpi_utils.h"

int elem_in_list_elem(int elem, int *list_elem, int size_list_elem) {
    /* Return the index of elem if it is in list_elem, -1 otherwise */
    int i;
    for (i = 0; i < size_list_elem; i++) {
        if (elem == list_elem[i])
            return i;
    }
    return -1;
}

int mpi_send_data_from_list_rank(int *list_rank_sender, int *list_rank_receiver,
                                 int size_list_rank, void *data_init,
                                 size_t full_size_data, void *data_out, int tag,
                                 MPI_Comm world_comm) {
    /* Redistribute data stored in list_rank_sender so that every processes with
       rank in list_rank_sender will send their data to processes in
       list_rank_receiver For a given i, this means that the process with rank
       list_rank_sender[i] will send its data to the process
       list_rank_receiver[i]

        This routine will be mostly used for butterfly scheme purposes

        BEWARE : full_size_data is expected in BYTE (typically as a result of
       the application sizeof)
    */
    int i;
    int *ranks_not_const;

    int rank, size, difference_treshold; //, tag = 0;
    void *buffer;

    MPI_Comm_rank(world_comm, &rank);
    MPI_Comm_size(world_comm, &size);
    if (full_size_data == 0)
        return 0;

    MPI_Request s_request, r_request;
    MPI_Status mpi_status;

    int index_rank_sender =
        elem_in_list_elem(rank, list_rank_sender, size_list_rank);
    int index_rank_receiver =
        elem_in_list_elem(rank, list_rank_receiver, size_list_rank);

    if ((index_rank_sender != -1) || (index_rank_receiver != -1)) {
        if (index_rank_sender != -1) {
            // MPI_Isend(buffer, full_size_data, MPI_BYTE,
            // list_rank_receiver[index_rank_sender], tag, world_comm,
            // &s_request);
            MPI_Send(data_init, full_size_data, MPI_BYTE,
                     list_rank_receiver[index_rank_sender],
                     tag + index_rank_sender, world_comm);

            // MPI_Wait(&s_request, &mpi_status);
        }
        if (index_rank_receiver != -1) {
            buffer = malloc(full_size_data);
            // MPI_Irecv(buffer, full_size_data, MPI_BYTE,
            // list_rank_sender[index_rank_receiver], tag, world_comm,
            // &r_request);
            MPI_Recv(buffer, full_size_data, MPI_BYTE,
                     list_rank_sender[index_rank_receiver],
                     tag + index_rank_receiver, world_comm, &mpi_status);
            memcpy(data_out, buffer, full_size_data);
            // fflush(stdout);
            // MPI_Wait(&r_request, &mpi_status);
            free(buffer);
        }
    }
    return 0;
}

int mpi_create_subset(int number_ranks_to_divive, MPI_Comm initcomm,
                      MPI_Comm *subset_comm) {
    /* Create a mpi communicator subset of the initial global communicator, by
     * taking the number_ranks_to_divide first ranks within it*/
    int i;
    int *ranks_not_const;
    MPI_Group global_mpi_group, mpi_subset_group;

    int tag = 0;

    ranks_not_const = (int *)malloc(number_ranks_to_divive * sizeof(int));
    for (i = 0; i < number_ranks_to_divive; i++) {
        ranks_not_const[i] = i;
    }
    const int *ranks_const = ranks_not_const;

    // MPI_Comm_split(initcomm, color, initrank, &s2hat_comm);

    // Get the group of the whole communicator
    MPI_Comm_group(initcomm, &global_mpi_group);
    // MPI_Comm_dup(initcomm, &global_mpi_group);

    // Construct group containing all ranks under number_ranks_to_divive
    // MPI_Comm_group(initcomm, &mpi_subset_group);
    MPI_Group_incl(global_mpi_group, number_ranks_to_divive, ranks_const,
                   &mpi_subset_group);

    MPI_Comm_create_group(initcomm, mpi_subset_group, tag, subset_comm);

    MPI_Group_free(&global_mpi_group);
    MPI_Group_free(&mpi_subset_group);

    return 0;
}
