#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
#include <unistd.h>

#include "s2hat.h"
// #include "midapack.h"
#include "s2hat_tools.h"


int elem_in_list_elem(int elem, int *list_elem, int size_list_elem){
  /* Return the index of elem if it is in list_elem, -1 otherwise */
  int i;
  for(i=0; i<size_list_elem; i++){
    if (elem == list_elem[i])
      return i;
  }
  return -1;
}


int mpi_send_data_from_list_rank(int *list_rank_sender, int *list_rank_receiver, int size_list_rank, void *data_init, size_t full_size_data, void *data_out, MPI_Comm world_comm)
{
    /* Redistribute data stored in list_rank_sender so that every processes with rank in list_rank_sender will send their data to processes
        in list_rank_receiver
        For a given i, this means that the process with rank list_rank_sender[i] will send its data to the process list_rank_receiver[i]

        This routine will be mostly used for butterfly scheme purposes

        BEWARE : full_size_data is expected in BYTE (typically as a result of the application sizeof)
    */
    int i;
    int *ranks_not_const;

    int rank, size, difference_treshold, tag = 0;
    void *buffer;

    MPI_Comm_rank(world_comm, &rank);
    MPI_Comm_size(world_comm, &size);
    MPI_Request s_request, r_request;

    int index_rank_sender = elem_in_list_elem(rank, list_rank_sender, size_list_rank);
    int index_rank_receiver = elem_in_list_elem(rank, list_rank_receiver, size_list_rank);


    if ((index_rank_sender != -1) && (index_rank_receiver != -1))
        {
            buffer = malloc(full_size_data);
            if (index_rank_sender != -1){
                // buffer = malloc(size_data*sizeof(double));
                memcpy(buffer, data_init, full_size_data);
                MPI_Isend(buffer, full_size_data, MPI_BYTE, list_rank_receiver[index_rank_sender], tag, world_comm, &s_request);
            }
            if (index_rank_receiver != -1){
                // buffer = malloc(size_data*sizeof(double));
                MPI_Irecv(buffer, full_size_data, MPI_BYTE, list_rank_sender[index_rank_receiver], tag, world_comm, &r_request);
                memcpy(data_out, buffer, full_size_data);
            }
            free(buffer);
        }
    return 0;
}



// int mpi_send_data_from_above_treshold(int treshold_rank, double *data_init, int size_data, double *data_out, MPI_Comm world_comm)
// {
//     /* Redistribute data stored in local processes so that every processes with rank above the treshold_rank will send their data to processes
//         symmetrically to the treshold_rank
//         This means that the process with rank treshold_rank+1 will send its data to the process treshold_rank-1

//         It is always assumed that the total number of processes in world_comm is strictly more than treshold_rank

//         This routine will be mostly used for butterfly scheme purposes
//     */
//     int i;
//     int *ranks_not_const;

//     int rank, size, difference_treshold, tag = 0;
//     double *buffer;


//     MPI_Comm_rank(world_comm, &rank);
//     MPI_Comm_size(world_comm, &size);
//     MPI_Request s_request, r_request;

//     if (rank >= treshold_rank - (size - treshold_rank))
//         {
//             buffer = malloc(size_data*sizeof(double));
//             if (rank > treshold_rank){
//                 difference_treshold = rank - treshold_rank;

//                 // buffer = malloc(size_data*sizeof(double));
//                 memcpy(buffer, data_init, size_data*sizeof(double));
//                 MPI_Isend(buffer, size_data, MPI_DOUBLE, treshold_rank-difference_treshold, tag, world_comm, &s_request);
//                 // free(buffer);
//             }
//             if (rank < treshold_rank){
//                 difference_treshold = treshold_rank - rank;

//                 // buffer = malloc(size_data*sizeof(double));
//                 MPI_Irecv(buffer, size_data, MPI_DOUBLE, treshold_rank + difference_treshold, tag, world_comm, &r_request);
//                 memcpy(data_out, buffer, size_data*sizeof(double));
//                 // free(buffer);
//             }
//             free(buffer);
//         }
    
//     return 0;
// }

// int mpi_send_data_from_below_treshold(int treshold_rank, double *data_init, int size_data, double *data_out, MPI_Comm world_comm)
// {
//     /* Redistribute data stored in local processes so that every processes with rank below the treshold_rank, within the difference size-treshold_rank, 
//         will send their data to processes symmetrically to the treshold_rank
//         This means that the process with rank treshold_rank-1 will send its data to the process treshold_rank+1

//         It is always assumed that the total number of processes in world_comm is strictly more than treshold_rank

//         This routine will be mostly used for butterfly scheme purposes
//     */
//     int i;
//     int *ranks_not_const;

//     int rank, size, difference_treshold, tag = 0;
//     double *buffer;


//     MPI_Comm_rank(world_comm, &rank);
//     MPI_Comm_size(world_comm, &size);
//     MPI_Request s_request, r_request;

//     if (rank >= treshold_rank - (size - treshold_rank))
//         {
//             buffer = malloc(size_data*sizeof(double));
//             if (rank > treshold_rank){
//                 difference_treshold = rank - treshold_rank;

//                 memcpy(buffer, data_init, size_data*sizeof(double));
//                 MPI_Irecv(buffer, size_data, MPI_DOUBLE, treshold_rank-difference_treshold, tag, world_comm, &s_request);
//             }
//             if (rank < treshold_rank){
//                 difference_treshold = treshold_rank - rank;

//                 MPI_Isend(buffer, size_data, MPI_DOUBLE, treshold_rank + difference_treshold, tag, world_comm, &r_request);
//                 memcpy(data_out, buffer, size_data*sizeof(double));
//             }
//             free(buffer);
//         }
    
//     return 0;
// }

// int butterfly_mirroring(double *values_to_send, int *indices_to_send, int number_elements_to_send, double *values_to_receive, int *indices_to_receive, MPI_Comm worldcomm){
//   /* Sends data of the excess processes, which are over the 2^k processes which will be used for the Butterfly scheme,
//      back to the 2^k processes. The way to proceed is by taking the data of the N-2^k processes and send them respectively to the last
//      N-2^k processes which will be used for the Butterfly scheme. This means the data of the process ranked 2^k will be send to the process ranked 2^k-1,
//      and the data of the process ranked 2^k+10 will be send to the process ranked 2^k-11 
//      This function should be followed by butterfly_init_reorder to rearrange the indices and associated values, to make sure
//      the received indices are order monotonously */

//   int i, rank, size;
//   MPI_Request s_request, r_request, s_request_1, r_request_1, s_request_2, r_request_2;
//   int tag = 0;

//   double *value_communicated;
//   int *indices_communicated;

//   int number_steps_Butterfly = log_2(size);
//   int number_rank_used_Butterfly = pow(2,number_steps_Butterfly);
//   int number_elements_to_receive;
//   int number_ranks_to_send = size - number_rank_used_Butterfly;
  
//   int rank_excess = fabs(rank - number_rank_used_Butterfly);

//   MPI_Comm_size(worldcomm, &size);
//   MPI_Comm_rank(worldcomm, &rank);
  
  
//   if (rank_excess <= number_ranks_to_send){
    
//     // First, send the number of elements to send/receive
//     if (rank < number_rank_used_Butterfly){
//       // MPI_Irecv(&number_elements_to_receive, 1, MPI_INT, rank + rank_excess, tag, worldcomm, &r_request);
//       MPI_Recv(&number_elements_to_receive, 1, MPI_INT, rank + rank_excess, tag, worldcomm, &r_request);
//     }
//     if (rank >= number_rank_used_Butterfly){
//       MPI_Send(&number_elements_to_send, 1, MPI_INT, rank - rank_excess - 1, tag, worldcomm);
//     }

//     // MPI_Wait(&r_request, MPI_STATUS_IGNORE);
//     // MPI_Wait(&s_request, MPI_STATUS_IGNORE);

//     // Then, send the data
//     if (rank < number_rank_used_Butterfly){
//       value_communicated = (double *)malloc(number_elements_to_receive*sizeof(double));
//       indices_communicated = (int *)malloc(number_elements_to_receive*sizeof(int));

//       MPI_Irecv(&value_communicated, number_elements_to_receive, MPI_DOUBLE, rank + rank_excess, tag+1, worldcomm, &r_request_1);
//       MPI_Irecv(&value_communicated, number_elements_to_receive, MPI_INT, rank + rank_excess, tag+2, worldcomm, &r_request_2);
//       MPI_Wait(&r_request_1, &r_request);
//       MPI_Wait(&r_request_2, &r_request);

//       memcpy(values_to_receive, value_communicated, number_elements_to_receive*sizeof(double));
//       memcpy(indices_to_receive, indices_communicated, number_elements_to_receive*sizeof(int));
//     }

//     if (rank >= number_rank_used_Butterfly){
//       value_communicated = (double *)malloc(number_elements_to_receive*sizeof(double));
//       indices_communicated = (int *)malloc(number_elements_to_receive*sizeof(int));
      
//       memcpy(value_communicated, values_to_send, number_elements_to_receive*sizeof(double));
//       memcpy(indices_communicated, indices_to_send, number_elements_to_receive*sizeof(int));
      
//       MPI_Isend(&values_to_send, number_elements_to_send, MPI_DOUBLE, rank - rank_excess - 1, tag+1, worldcomm, &s_request_1);
//       MPI_Isend(&indices_to_send, number_elements_to_send, MPI_INT, rank - rank_excess - 1, tag+2, worldcomm, &s_request_2);
//       MPI_Wait(&s_request_1, &r_request);
//       MPI_Wait(&s_request_2, &r_request);
//     }

//     // MPI_Wait(&r_request, MPI_STATUS_IGNORE);
//     // MPI_Wait(&s_request, MPI_STATUS_IGNORE);
//     free(value_communicated);
//     free(indices_communicated);
//   }
// }

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
