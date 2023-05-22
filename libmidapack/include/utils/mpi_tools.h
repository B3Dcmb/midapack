 /** @file mpi.h
    @brief <b> Declaration of the backbone routines for general MPI tools.</b>
    @author Magdy Morshed
    @date May 2023 */
 
    
#ifdef W_MPI
#include <mpi.h>
#endif

// #include "midapack.h"


#ifndef DBL_MAX
#define DBL_MAX            1.79769313486231470e+308
#endif

/* Content of mpi_tools.c */

/* Get elem in list_elem and return its index */
int elem_in_list_elem(int elem, int *list_elem, int size_list_elem);

/* Create a mpi communicator subset of the initial global communicator, by taking the number_ranks_to_divide first ranks within it*/
int mpi_create_subset(int number_ranks_to_divive, MPI_Comm initcomm, MPI_Comm *subset_comm);
/**/
/* Send data between processors of list_rank_sender and list_rank_receiver, at equivalent indexes */
int mpi_send_data_from_list_rank(int *list_rank_sender, int *list_rank_receiver, int size_list_rank, void *data_init, size_t full_size_data, void *data_out, int tag, MPI_Comm world_comm);
