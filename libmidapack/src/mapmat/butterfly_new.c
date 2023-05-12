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

#include "als.h"
#include "alm.h"
#include "butterfly.h"


int butterfly_reduce_init(int *indices, int count, int **R, int *nR, int **S, int *nS, int **com_indices, int *com_count, int steps, MPI_Comm comm){

  /* initializes butterfly communication where the pixels distribution does not change at the outset of the process *
   * This is the original MIDAPACK routine by Pierre Cargemel ...                                   - rs 2022/06/10 */

  int i, k, p2k;
  int rank, size, rk, sk;
  int tag;
  MPI_Request s_request, r_request;
  int nbuf, *buf;
  int **I, *nI;
  int **J, *nJ;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  I = (int **) malloc(steps * sizeof(int*));
  nI = (int *) malloc(steps * sizeof(int));
  tag=0;
  p2k=size/2;

  for(k=0; k<steps; k++){		/* butterfly first pass : bottom up (fill tabs nI and I) */
    sk=(rank+size-p2k)%size;
    rk=(rank+p2k)%size;

    if(k==0){     					/* S^0 := A */
      nS[k] = count;                                    /* NEEDS TO BE MODIFIED with "final" number of pixels */
      S[k] = (int *) malloc(nS[k] * sizeof(int));
      memcpy( S[k], indices, nS[k]*sizeof(int));        /* copy *my* pixel indices to the send (really receive ?!) buffer NEEDS TO BE MODIFIED with final pix numbers ?! */
    }
    else{      						/* S^k := S^{k-1} \cup R^{k-1} */
      nS[k] = modified_card_or(S[k-1], nS[k-1], I[steps-k], nI[steps-k]);
      S[k] = (int *) malloc(nS[k] * sizeof(int));
      modified_set_or(S[k-1], nS[k-1], I[steps-k], nI[steps-k], S[k]);
    }

    MPI_Irecv(&nI[steps-k-1], 1, MPI_INT, rk, tag, comm, &r_request);	/* receive number of indices */
    MPI_Isend(&nS[k], 1, MPI_INT, sk, tag, comm, &s_request);		/* send number of indices */
    MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    I[steps-k-1]= (int *) malloc(nI[steps-k-1] * sizeof(int));

    tag++;
    MPI_Irecv(I[steps-k-1], nI[steps-k-1], MPI_INT, rk, tag, comm, &r_request);	/* receive indices */
    MPI_Isend(S[k], nS[k], MPI_INT, sk, tag, comm, &s_request);			/* send indices */
    MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    p2k/=2;
    tag++;
  }

  J = (int **) malloc(steps * sizeof(int*));
  nJ = (int *) malloc(steps * sizeof(int));

  tag=0;
  p2k=1;
  for(k=0; k<steps; k++){		/* buuterfly second pass : top down (fill tabs nJ and J) */
    free(S[k]);
    sk=(rank+p2k)%size;
    rk=(rank+size-p2k)%size;
    if(k==0){
      nJ[k] = count;
      J[k] = (int *) malloc(nJ[k] * sizeof(int));
      memcpy( J[k], indices, nJ[k]*sizeof(int));
    }
    else{
      nJ[k] = modified_card_or(J[k-1], nJ[k-1], R[k-1], nR[k-1]);
      J[k] = (int *) malloc(nJ[k] * sizeof(int));
      modified_set_or(J[k-1], nJ[k-1], R[k-1], nR[k-1], J[k]);  /* J^k=R^k-1 \cup J^k-1 */
      free(R[k-1]);
    }
    if(k!=steps-1){
    MPI_Irecv(&nR[k], 1, MPI_INT, rk, tag, comm, &r_request);
    MPI_Isend(&nJ[k], 1, MPI_INT, sk, tag, comm, &s_request);
    MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    R[k]= (int *) malloc( nR[k] * sizeof(int));
    tag++;

    MPI_Irecv(R[k], nR[k], MPI_INT, rk, tag, comm, &r_request);
    MPI_Isend(J[k], nJ[k], MPI_INT, sk, tag, comm, &s_request);
    MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    MPI_Wait(&s_request, MPI_STATUS_IGNORE);
    }
    p2k*=2;
    tag++;
  }


  tag=0;
  p2k=1;
  for(k=0; k<steps; k++){		/* butterfly last pass : know that Sending tab is S = I \cap J, so send S and we'll get R */
    sk=(rank+p2k)%size;
    rk=(rank+size-p2k)%size;

    nS[k] = card_and(I[k], nI[k], J[k], nJ[k]);
    S[k] = (int *) malloc(nJ[k] *sizeof(int));
    set_and( I[k], nI[k], J[k], nJ[k], S[k]);	/* S^k=I^k \cap J^k */

    free(I[k]);
    free(J[k]);

    MPI_Irecv(&nR[k], 1, MPI_INT, rk, tag, comm, &r_request);	/* receive size */
    MPI_Isend(&nS[k], 1, MPI_INT, sk, tag, comm, &s_request);	/* send size */
    MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    R[k]= (int *) malloc( nR[k] * sizeof(int));
    tag++;

    MPI_Irecv(R[k], nR[k], MPI_INT, rk, tag, comm, &r_request); /* receive indices */
    MPI_Isend(S[k], nS[k], MPI_INT, sk, tag, comm, &s_request); /* send indices */
    MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    p2k*=2;
    tag++;
  }

  /* Now we work locally */
  int **USR, *nUSR, **U, *nU;

  USR = (int **) malloc(steps*sizeof(int *));
  nUSR = (int *) malloc(steps*sizeof(int));
  U = (int **) malloc(steps*sizeof(int *));
  nU = (int *) malloc(steps*sizeof(int));

  for(k=0; k<steps; k++){
    nUSR[k] = modified_card_or(S[k], nS[k], R[k], nR[k]);
    USR[k] = (int *) malloc(nUSR[k]*sizeof(int));
    modified_set_or(S[k], nS[k], R[k], nR[k], USR[k]);
  }
  for(k=0; k<steps; k++){
    if(k==0){
      nU[k]=nUSR[k];
      U[k] = (int *) malloc(nU[k] * sizeof(int));
      memcpy( U[k], USR[k], nU[k]*sizeof(int));
    }
    else{
      nU[k] = modified_card_or(U[k-1], nU[k-1], USR[k], nUSR[k]);
      U[k] = (int *) malloc(nU[k]*sizeof(int *));
      modified_set_or(U[k-1], nU[k-1], USR[k], nUSR[k], U[k]);
    }
  }
  *com_count=nU[steps-1];
  *com_indices = (int *) malloc(*com_count * sizeof(int));
  memcpy(*com_indices, U[steps-1], *com_count * sizeof(int));
  /* ==================================================================== */

  for(k=0; k<steps; k++){
    subset2map(*com_indices, *com_count, S[k], nS[k]);
    subset2map(*com_indices, *com_count, R[k], nR[k]);
  }
  free(USR);
  free(U);

 return 0;
}

int butterfly_reshuffle_init(int *indices_in, int count_in, int *indices_out, int count_out, int **R, int *nR, int **S, int *nS, int **com_indices, int *com_count, int steps, MPI_Comm comm){

  /* initializes butterfly communication where the pixels distributions prior to and after it may be different.       *
   * The routine explicitly allows for some of the MPI processes to have no data either prior or after.               *
   * This is a recoding of the original MIDAPACK routine by Pierre Cargemel (called now butterfly_reduce_init         *
   * which should be equivalent to if the input and output distribution of the pixels coincide.       - rs 2022/06/10 */


  int i, k, p2k;
  int rank, size, rk, sk;
  int tag;
  MPI_Request s_request, r_request;
  int nbuf, *buf;
  int **I, *nI;
  int **J, *nJ;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  I = (int **) malloc(steps * sizeof(int*));
  nI = (int *) malloc(steps * sizeof(int));
  tag=0;
  p2k=size/2;

  for(k=0; k<steps; k++){		/* butterfly first pass : bottom up (fill tabs nI and I) */
    sk=(rank+size-p2k)%size;
    rk=(rank+p2k)%size;

    if(k==0){     					/* S^0 := A */
      nS[k] = count_out;
      if( nS[k]) {                                      /* do not allocate zero length objects */
         S[k] = (int *) malloc(nS[k] * sizeof(int));
         memcpy( S[k], indices_out, nS[k]*sizeof(int));
      }
    }
    else{      						/* S^k := S^{k-1} \cup R^{k-1} */
      nS[k] = modified_card_or(S[k-1], nS[k-1], I[steps-k], nI[steps-k]);
      if( nS[k]) {
	    S[k] = (int *) malloc(nS[k] * sizeof(int));
        modified_set_or(S[k-1], nS[k-1], I[steps-k], nI[steps-k], S[k]);
      }
    }

    MPI_Irecv(&nI[steps-k-1], 1, MPI_INT, rk, tag, comm, &r_request);	/* receive number of indices */
    MPI_Isend(&nS[k], 1, MPI_INT, sk, tag, comm, &s_request);		/* send number of indices */
    MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    if(nI[steps-k-1]) I[steps-k-1]= (int *) malloc(nI[steps-k-1] * sizeof(int));

    tag++;
    if( nI[steps-k-1]) MPI_Irecv(I[steps-k-1], nI[steps-k-1], MPI_INT, rk, tag, comm, &r_request);	/* receive indices */
    if( nS[k]) MPI_Isend(S[k], nS[k], MPI_INT, sk, tag, comm, &s_request);			        /* send indices */
    if( nI[steps-k-1]) MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    if( nS[k]) MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    p2k/=2;
    tag++;
  }

  J = (int **) malloc(steps * sizeof(int*));
  nJ = (int *) malloc(steps * sizeof(int));

  tag=0;
  p2k=1;
  for(k=0; k<steps; k++){		/* buuterfly second pass : top down (fill tabs nJ and J) */
    free(S[k]);
    sk=(rank+p2k)%size;
    rk=(rank+size-p2k)%size;
    if(k==0){
      nJ[k] = count_in;
      if( nJ[k]) {
         J[k] = (int *) malloc(nJ[k] * sizeof(int));
         memcpy( J[k], indices_in, nJ[k]*sizeof(int));
      }
    }
    else{
      nJ[k] = modified_card_or(J[k-1], nJ[k-1], R[k-1], nR[k-1]);
      if(nJ[k]) {
         J[k] = (int *) malloc(nJ[k] * sizeof(int));
         modified_set_or(J[k-1], nJ[k-1], R[k-1], nR[k-1], J[k]);  /* J^k=R^k-1 \cup J^k-1 */
      }
      free(R[k-1]);
    }
    if(k!=steps-1){
       MPI_Irecv(&nR[k], 1, MPI_INT, rk, tag, comm, &r_request);
       MPI_Isend(&nJ[k], 1, MPI_INT, sk, tag, comm, &s_request);
       MPI_Wait(&r_request, MPI_STATUS_IGNORE);
       MPI_Wait(&s_request, MPI_STATUS_IGNORE);

       if( nR[k]) R[k]= (int *) malloc( nR[k] * sizeof(int));
       tag++;

       if( nR[k]) MPI_Irecv(R[k], nR[k], MPI_INT, rk, tag, comm, &r_request);
       if( nJ[k]) MPI_Isend(J[k], nJ[k], MPI_INT, sk, tag, comm, &s_request);
       if( nR[k]) MPI_Wait(&r_request, MPI_STATUS_IGNORE);
       if( nJ[k]) MPI_Wait(&s_request, MPI_STATUS_IGNORE);
    }
    p2k*=2;
    tag++;
  }


  tag=0;
  p2k=1;
  for(k=0; k<steps; k++){		/* butterfly last pass : know that Sending tab is S = I \cap J, so send S and we'll get R */
    sk=(rank+p2k)%size;
    rk=(rank+size-p2k)%size;

    nS[k] = card_and(I[k], nI[k], J[k], nJ[k]);
    if( nS[k]) {
       S[k] = (int *) malloc(nJ[k] *sizeof(int));
       set_and( I[k], nI[k], J[k], nJ[k], S[k]);	/* S^k=I^k \cap J^k */
    }

    if(nI[k]) free(I[k]);
    if(nJ[k]) free(J[k]);

    MPI_Irecv(&nR[k], 1, MPI_INT, rk, tag, comm, &r_request);	/* receive size */
    MPI_Isend(&nS[k], 1, MPI_INT, sk, tag, comm, &s_request);	/* send size */
    MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    if( nR[k]) R[k]= (int *) malloc( nR[k] * sizeof(int));
    tag++;

    if( nR[k]) MPI_Irecv(R[k], nR[k], MPI_INT, rk, tag, comm, &r_request); /* receive indices */
    if( nS[k]) MPI_Isend(S[k], nS[k], MPI_INT, sk, tag, comm, &s_request); /* send indices */
    if( nR[k]) MPI_Wait(&r_request, MPI_STATUS_IGNORE);
    if( nS[k]) MPI_Wait(&s_request, MPI_STATUS_IGNORE);

    p2k*=2;
    tag++;
  }

  /* Now we work locally */
  int **USR, *nUSR, **U, *nU;

  USR = (int **) malloc(steps*sizeof(int *));
  nUSR = (int *) malloc(steps*sizeof(int));
  U = (int **) malloc(steps*sizeof(int *));
  nU = (int *) malloc(steps*sizeof(int));

  for(k=0; k<steps; k++){
    nUSR[k] = modified_card_or(S[k], nS[k], R[k], nR[k]);
    if( nUSR[k]) {
       USR[k] = (int *) malloc(nUSR[k]*sizeof(int));
       modified_set_or(S[k], nS[k], R[k], nR[k], USR[k]);
    }
  }
  for(k=0; k<steps; k++){
    if(k==0){
      nU[k]=nUSR[k];
      if( nU[k]) {
	    U[k] = (int *) malloc(nU[k] * sizeof(int));
        memcpy( U[k], USR[k], nU[k]*sizeof(int));
      }
    }
    else{
      nU[k] = modified_card_or(U[k-1], nU[k-1], USR[k], nUSR[k]);
      if( nU[k]) {
        U[k] = (int *) malloc(nU[k]*sizeof(int *));
        modified_set_or(U[k-1], nU[k-1], USR[k], nUSR[k], U[k]);
      }
    }
  }
  *com_count=nU[steps-1];
  if( *com_count) {
     *com_indices = (int *) malloc(*com_count * sizeof(int));
     memcpy(*com_indices, U[steps-1], *com_count * sizeof(int));   /* the full set of pixel indices dealt with by this proc at some point */
  }
  /* ==================================================================== */

  /* makes the indices in S and R relative to those stored in full set of pixels processed by this proc */

  if( *com_count) {
    for(k=0; k<steps; k++){
       if(nS[k]) subset2map(*com_indices, *com_count, S[k], nS[k]);
       if(nR[k]) subset2map(*com_indices, *com_count, R[k], nR[k]);
    }
  }
  free(USR);
  free(U);

 return 0;
}



int mirror_butterfly(double *values_local, int *indices_local, int size_local, double *values_received, int *indices_received, int *size_received, int flag_mirror_unmirror_size_indices_data, MPI_Comm world_comm)
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
    int nb_rank_Butterfly = pow(2,number_steps_Butterfly);

    int number_ranks_to_send = size - nb_rank_Butterfly;

    int *list_rank_excess_butterfly = (int *)malloc(number_ranks_to_send*sizeof(int));
    int *list_last_ranks_within_butterfly = (int *)malloc(number_ranks_to_send*sizeof(int));

    for (i=0; i<number_ranks_to_send; i++){
        list_rank_excess_butterfly[i] = nb_rank_Butterfly + i;
        list_last_ranks_within_butterfly[i] = nb_rank_Butterfly - i - 1;
    }

    // First, communicate the sizes, determine which processes will receive/send data
    int size_communicated = 0;
    switch(flag_mirror_unmirror_size_indices_data)
    {
        case MIRROR_SIZES: // Case 0 : MIRROR - send only size and indices to be exchanged, to have the excess processes over 2^k to send their data to the last 2^k processes before the Butterfly scheme
            mpi_send_data_from_list_rank(list_rank_excess_butterfly, list_last_ranks_within_butterfly, number_ranks_to_send, size_local, 1*sizeof(int), &size_communicated, world_comm);
            *size_received = size_communicated;
            break;

        case MIRROR_INDICES:
            size_communicated = *size_received;
            mpi_send_data_from_list_rank(list_rank_excess_butterfly, list_last_ranks_within_butterfly, number_ranks_to_send, indices_local, size_communicated*sizeof(int), indices_received, world_comm);
            break;

        case MIRROR_DATA:
            size_communicated = *size_received;
            mpi_send_data_from_list_rank(list_rank_excess_butterfly, list_last_ranks_within_butterfly, number_ranks_to_send, values_local, size_communicated*sizeof(double), values_received, world_comm);
            break;

        case UNMIRROR_SIZE:
            mpi_send_data_from_list_rank(list_last_ranks_within_butterfly, list_rank_excess_butterfly, number_ranks_to_send, size_local, 1*sizeof(int), &size_communicated, world_comm);
            *size_received = size_communicated;
            break;

        case UNMIRROR_INDICES:
            size_communicated = *size_received;
            mpi_send_data_from_list_rank(list_last_ranks_within_butterfly, list_rank_excess_butterfly, number_ranks_to_send, indices_local, size_communicated*sizeof(int), indices_received, world_comm);
            break;
        
        case UNMIRROR_DATA: // Case unmirror : the last processes within the 2^k butterfly processes send their data to the excess processes over the 2^k processes after the Butterfly scheme
            size_communicated = *size_received;
            mpi_send_data_from_list_rank(list_last_ranks_within_butterfly, list_rank_excess_butterfly, number_ranks_to_send, values_local, size_communicated*sizeof(double), values_received, world_comm);
            break;
    }

    free(list_rank_excess_butterfly);
    free(list_last_ranks_within_butterfly);
}


int construct_butterfly_struct(Butterfly_struct *Butterfly_obj, int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, int do_we_need_to_project_into_different_scheme, MPI_Comm worlcomm)
{
    /* Initialize the butterfly communication, assuming the mirroring of the data have been already done */
    int size;
    MPI_Comm_size(comm, &size);

    Butterfly_obj->classic_or_reshuffle_butterfly = flag_classic_or_reshuffle_butterfly;
    // 0 if classic butterfly, e.g. if pixel distributions are the same before/after ; 1 if reshuffle butterfly

    Butterfly_obj->steps = log_2(size);
    Butterfly_obj->S = (int **)malloc(Butterfly_obj->steps * sizeof(int *)); // allocate sending maps tab
    Butterfly_obj->R = (int **)malloc(Butterfly_obj->steps * sizeof(int *)); // allocate receiving maps tab
    Butterfly_obj->nS = (int *)malloc(Butterfly_obj->steps * sizeof(int));   // allocate sending map sizes tab
    Butterfly_obj->nR = (int *)malloc(Butterfly_obj->steps * sizeof(int));   // allocate receiving map size tab
    // butterfly_init(A->lindices + (A->nnz) * (A->trash_pix), A->lcount - (A->nnz) * (A->trash_pix), A->R, A->nR, A->S, A->nS, &(A->com_indices), &(A->com_count), A->steps, comm);
    
    // Construct the butterfly communicator
    MPI_Comm comm_butterfly;
    mpi_create_subset(pow(2,log_2(size)), world_comm, &comm_butterfly);

    switch(flag_classic_or_reshuffle_butterfly)
    {
        case 0: // Classic butterfly
            butterfly_reduce_init(indices_in, count_in, Butterfly_obj->R, Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS, &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count), Butterfly_obj->steps, comm_butterfly);

        case 1: // Reshuffle butterfly
            butterfly_reshuffle_init(indices_in, count_in, indices_out, count_out, Butterfly_obj->R, Butterfly_obj->nR, Butterfly_obj->S, Butterfly_obj->nS, &(Butterfly_obj->com_indices), &(Butterfly_obj->com_count), Butterfly_obj->steps, comm_butterfly);
    }

    // Butterfly_obj->do_we_need_to_project_into_different_scheme = do_we_need_to_project_into_different_scheme;
    return 0;
}


int modified_butterfly_reduce(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax, double *val, int steps, MPI_Comm comm){
  /* double st, t; */
  /* t=0.0; */

  /* this performs the butterfly all reduce as in the original routine of MIDAPACK butterfly() by P. Cargemel but allows for MPI *
   * process which have no data at the begining or at the end of the communication.                                              *
   * This routine needs to be preceded by a call to either butterfly_reduce_init() if the pixel distributions prior and after    *
   * communication are the same or butterfly_reshuffle_init() if they are not.                                   - rs 2022/06/09 */

  int k, p2k, tag;
  int rank, size, rk, sk;
  MPI_Request s_request, r_request;
  double *sbuf, *rbuf;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  sbuf = (double *) malloc(nSmax * sizeof(double));
  rbuf = (double *) malloc(nRmax * sizeof(double));
  tag=0;
  p2k=1;

  for(k=0; k<steps; k++){
    /* st=MPI_Wtime(); */
    rk=(rank+size-p2k)%size;
    if( nR[k]) MPI_Irecv(rbuf, nR[k], MPI_DOUBLE, rk, tag, comm, &r_request);
    sk=(rank+p2k)%size;
    if( nS[k]) {
       m2s(val, sbuf, S[k], nS[k]); /* fill the sending buffer */
       MPI_Isend(sbuf, nS[k], MPI_DOUBLE, sk, tag, comm, &s_request);
    }
    if( nR[k]) {
      MPI_Wait(&r_request, MPI_STATUS_IGNORE);
      s2m_sum(val, rbuf, R[k], nR[k]); /* sum receive buffer into values nR[k] floating sum */
    }
    p2k*=2;
    tag++;
    if( nS[k]) MPI_Wait(&s_request, MPI_STATUS_IGNORE);
    /* t=t+MPI_Wtime()-st; */
  }
  free(sbuf);
  free(rbuf);
  return 0;
}

int butterfly_reshuffle(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax, double *val, int steps, MPI_Comm comm){
  /* double st, t; */
  /* t=0.0; */

  /* this performs the butterfly reshuffle of all the pixels without modifying the values assigned to the pixels. In the case of the redundant pixels *
   * distribution in the input it is assumed that it is consistent between different MPI process. This is not checked by the routine! In such a case  *
   * the routine may not be optimal either. It should be however if the initial distribution is without any overlaps beyween the processes. The       *
   * allows for some of the mprocess to have no data either on the input or output.                                                                   *
   * This routine is a recoded MIDAPACK routine, butterfly(), by P. Cargamel.                                                                         *
   * This routine needs to be preceded by a call to butterfly_reshuffle_init().                                                       - rs 2022/06/09 */

  int k, p2k, tag;
  int rank, size, rk, sk;
  MPI_Request s_request, r_request;
  double *sbuf, *rbuf;

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  sbuf = (double *) malloc(nSmax * sizeof(double));
  rbuf = (double *) malloc(nRmax * sizeof(double));
  tag=0;
  p2k=1;

  for(k=0; k<steps; k++){
    /* st=MPI_Wtime(); */
    rk=(rank+size-p2k)%size;
    if( nR[k]) MPI_Irecv(rbuf, nR[k], MPI_DOUBLE, rk, tag, comm, &r_request);
    sk=(rank+p2k)%size;
    if( nS[k]) {
       m2s(val, sbuf, S[k], nS[k]); /* fill the sending buffer */
       MPI_Isend(sbuf, nS[k], MPI_DOUBLE, sk, tag, comm, &s_request);
    }
    if( nR[k]) {
      MPI_Wait(&r_request, MPI_STATUS_IGNORE);
      s2m(val, rbuf, R[k], nR[k]); /* sum receive buffer into values nR[k] floating sum */
    }
    p2k*=2;
    tag++;
    if( nS[k]) MPI_Wait(&s_request, MPI_STATUS_IGNORE);
    /* t=t+MPI_Wtime()-st; */
  }
  free(sbuf);
  free(rbuf);
  return 0;
}

// Already in alm.c
// void m2s(double *mapval, double *submapval, int *subset, int count){
//   int i;
//   for(i=0; i< count; i++){
//     submapval[i]=mapval[subset[i]];
//   }
// }

// Already in alm.c
// void s2m_sum(double *mapval, double *submapval, int *subset, int count){
//   int i;
//   for(i=0; i< count; i++){
//     mapval[subset[i]] += submapval[i];
//   }
// }

// void s2m_copy(double *mapval, double *submapval, int *subset, int count){

//   int i;

//   for(i=0; i< count; i++){
//     mapval[subset[i]] = submapval[i];
//   }
// }

int modified_set_or(int *A1, int n1, int *A2, int n2, int *A1orA2){

  /* added cases when either n1 or n2 are zero. One of which *has to* be nonzero. - rs 2022/06/09 */

  int i=0, j=0, k= 0;

  if( n1 && n2) {
    while( i<n1 || j<n2){
      if(A1[i] < A2[j]){
         if(i<n1){
           A1orA2[k]=A1[i];
           i++;
         }
         else{
           A1orA2[k]=A2[j];
           j++;
         }
      }
      else if(A1[i] > A2[j]){
         if(j<n2){
            A1orA2[k]=A2[j];
            j++;
         }
         else{
           A1orA2[k]=A1[i];
           i++;
         }
      }
      else{
         A1orA2[k]=A1[i];
         i++;
         j++;
      }
      k++;
    }
  } else {
    if( n1 == 0) {
      for( j=0; j<n2; j++) A1orA2[j]=A2[j];
      k=n2;
    } else {
      for( i=0; i<n1; i++) A1orA2[i]=A1[i];
      k=n1;
    }
  }

  return k;
}

// Already in als.c
// int set_and(int *A1, int n1, int *A2, int n2, int *A1andA2){   /* this one is only called if the result is a non-empty set so both n1 and n2 have to be non-zero */
//   int i=0, j=0, k= 0;
//   while( i<n1 && j<n2){
//     if(A1[i] < A2[j]){
//       i++;
//     }
//     else if(A1[i] > A2[j]){
//       j++;
//     }
//     else{
//       A1andA2[k]=A1[i];
//       k++;
//       i++;
//       j++;
//     }
//   }
//   return k;
// }

int modified_card_or(int *A1, int n1, int *A2, int n2){

  /* acounts for cases with either n1 or n2 or both equal to zero - rs 2022/06/09 */

  int i=0, j=0, k= 0;

  if (n1 && n2) {
    while( i<n1 || j<n2){
      if(A1[i] < A2[j]){
        if(i<n1){ i++; }
        else{ j++; }
      }
      else if(A1[i] > A2[j]){
        if(j<n2){ j++; }
        else{ i++; }
      }
      else{
        if(i<n1){ i++; }
        if(j<n2){ j++; }
      }
      k++;
    }
  } else {
    // k = ( n1 == 0) ? k=n2 : k=n1;
    if (n1 == 0){
        k = n2;
    }
    else{
        k=n1;
    }
  }
  return k;
}

// Already in als.c
// void subset2map(int *A, int nA, int *subA, int nsubA){
//   int i=0, j=0;
//   while( i<nA && j<nsubA){
//     if(A[i] < subA[j]){
//       i++;
//     }
//     else{
//       subA[j]=i;
//       i++;
//       j++;
//     }
//   }
// }

#endif
