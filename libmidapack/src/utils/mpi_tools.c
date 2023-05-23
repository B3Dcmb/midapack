#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
#include <unistd.h>

#include "mpi_tools.h"
// #include "midapack.h"

int elem_in_list_elem(int elem, int *list_elem, int size_list_elem)
{
  /* Return the index of elem if it is in list_elem, -1 otherwise */
  int i;
  for(i=0; i<size_list_elem; i++){
    if (elem == list_elem[i])
      return i;
  }
  return -1;
}


int mpi_send_data_from_list_rank(int *list_rank_sender, int *list_rank_receiver, int size_list_rank, void *data_init, size_t full_size_data, void *data_out, int tag, MPI_Comm world_comm)
{
    /* Redistribute data stored in list_rank_sender so that every processes with rank in list_rank_sender will send their data to processes
        in list_rank_receiver
        For a given i, this means that the process with rank list_rank_sender[i] will send its data to the process list_rank_receiver[i]

        This routine will be mostly used for butterfly scheme purposes

        BEWARE : full_size_data is expected in BYTE (typically as a result of the application sizeof)
    */
    int i;
    int *ranks_not_const;

    int rank, size, difference_treshold; //, tag = 0;
    void *buffer;
    

    MPI_Comm_rank(world_comm, &rank);
    MPI_Comm_size(world_comm, &size);
    printf("func-mpi r %d ---- Starting ! size_t %zu \n", rank, full_size_data); fflush(stdout);
    if (full_size_data == 0)
        return 0;

    MPI_Request s_request, r_request;
    MPI_Status mpi_status;

    int index_rank_sender = elem_in_list_elem(rank, list_rank_sender, size_list_rank);
    int index_rank_receiver = elem_in_list_elem(rank, list_rank_receiver, size_list_rank);
    
    printf("func-mpi r %d ---- indexes found sender %d and receiver %d ; size_list_rank %d \n", rank, index_rank_sender, index_rank_receiver, size_list_rank); fflush(stdout);
    printf("func-mpi r %d ---- ranks found sender %d and receiver %d ; size_list_rank %d \n", rank, list_rank_sender[index_rank_receiver], list_rank_receiver[index_rank_sender], size_list_rank); fflush(stdout);
    

    if ((index_rank_sender != -1) || (index_rank_receiver != -1))
        {
            // buffer = malloc(full_size_data);
            if (index_rank_sender != -1){
                // buffer = malloc(size_data*sizeof(double));
                // memcpy(buffer, data_init, full_size_data);
                // MPI_Isend(buffer, full_size_data, MPI_BYTE, list_rank_receiver[index_rank_sender], tag, world_comm, &s_request);
                MPI_Send(data_init, full_size_data, MPI_BYTE, list_rank_receiver[index_rank_sender], tag+index_rank_sender, world_comm);
                printf("func-mpi r %d ---- send to %d with index %d ; size_t %zu \n", rank, list_rank_receiver[index_rank_sender], index_rank_sender, full_size_data);
                fflush(stdout);
                // MPI_Wait(&s_request, &mpi_status);
            }
            if (index_rank_receiver != -1){
                // buffer = malloc(size_data*sizeof(double));
                // MPI_Irecv(buffer, full_size_data, MPI_BYTE, list_rank_sender[index_rank_receiver], tag, world_comm, &r_request);
                MPI_Recv(data_out, full_size_data, MPI_BYTE, list_rank_sender[index_rank_receiver], tag+index_rank_receiver, world_comm, &mpi_status);
                int number_amount;
                MPI_Get_count(&mpi_status, MPI_DOUBLE, &number_amount);
                printf("func-mpi r %d ---- receive from %d with index %d ; size_t %zu \n", rank, list_rank_sender[index_rank_receiver], index_rank_receiver, full_size_data);
                printf("func-mpi r %d ---- received %d numbers from %d, tag = %d\n", rank, number_amount, mpi_status.MPI_SOURCE, mpi_status.MPI_TAG);
                // memcpy(data_out, buffer, full_size_data);
                // printf("2func-mpi r %d ---- receive from %d with index %d \n", rank, list_rank_sender[index_rank_receiver], index_rank_receiver);
                // fflush(stdout);
                // MPI_Wait(&r_request, &mpi_status);
            }
            printf("func-mpi TEST9 r %d \n", rank); fflush(stdout);
            // free(buffer);
            // printf("TEST10 r %d \n", rank); fflush(stdout);
        }
    return 0;
}



int mpi_create_subset(int number_ranks_to_divive, MPI_Comm initcomm, MPI_Comm *subset_comm)
{
    /* Create a mpi communicator subset of the initial global communicator, by taking the number_ranks_to_divide first ranks within it*/
    int i;
    int *ranks_not_const;
    MPI_Group global_mpi_group, mpi_subset_group;

    int tag = 0;
    
    ranks_not_const = (int *)malloc(number_ranks_to_divive*sizeof(int));
    for (i=0; i<number_ranks_to_divive; i++){
        ranks_not_const[i] = i;
    }
    const int* ranks_const = ranks_not_const;

    // MPI_Comm_split(initcomm, color, initrank, &s2hat_comm);
    
    // Get the group of the whole communicator
    MPI_Comm_group(initcomm, &global_mpi_group);
    // MPI_Comm_dup(initcomm, &global_mpi_group);

    // Construct group containing all ranks under number_ranks_to_divive
    // MPI_Comm_group(initcomm, &mpi_subset_group);
    MPI_Group_incl(global_mpi_group, number_ranks_to_divive, ranks_const, &mpi_subset_group);
    
    MPI_Comm_create_group(initcomm, mpi_subset_group, tag, subset_comm);

    MPI_Group_free(&global_mpi_group);
    MPI_Group_free(&mpi_subset_group);

    return 0;
}
