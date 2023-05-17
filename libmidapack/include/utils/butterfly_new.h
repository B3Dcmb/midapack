/** @file   butterfly_new.h
    @brief <b> Declaration of routines for butterfly-like communication scheme. </b>
    @note Copyright (C) 2010 APC CNRS Université Paris Diderot @nThis program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. @nYou should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/gpl.html
    @author Pierre Cargemel
    @date April 2012*/

#include <mpi.h>

#define MIRROR_SIZE 0
#define MIRROR_INDICES 1
#define MIRROR_DATA 2
#define UNMIRROR_SIZE 3
#define UNMIRROR_INDICES 4
#define UNMIRROR_DATA 5

typedef struct Butterfly_struct{

    int		*com_indices, com_count; // communicated indices, and size
    int		steps;			         // number of steps in the butterfly scheme
    int		*nS, *nR;		         // number of indices (to send and to receive); size = steps
    int		**R, **S;		         // sending or receiving indices

    int classic_or_reshuffle_butterfly; 
    // Flag to indicate if classic or reshuffled butterfly is used, e.g. if we expect the (pixel) distributions prior and after communication to be the same or different
    // 0 for classic butterfly scheme with same pixel distributions prior/after ; 1 for after

    MPI_Comm comm_butterfly;
    // MPI communicator used by the buttefly scheme
} Butterfly_struct;

typedef struct Butterfly_struct_supplement{

    // Indices obtained from mirroring
    int *indices_mirror; // List of indices received in the mirroring step by the local MPI process  -> Will be used afterwards to resend back the relevant indices when unmirroring
    int size_from_mirror; // Size of the indices obtained from mirroring


    int new_size_local; // New size after mirroring when the redundant indices have been deleted
    int *ordered_indices; // Ordered indices used for butterfly scheme, obtained from both mirroring and any conversion of the (pixel) distribution
    // This MUST BE ORDERED INDICES
} Butterfly_struct_supplement;


typedef struct Butterfly_superstruct{
    Butterfly_struct *Butterfly_obj;
    
    Butterfly_struct_supplement *Butterfly_mirror_supp;
    Butterfly_struct_supplement *Butterfly_unmirror_supp;

} Butterfly_superstruct;

/* Content of butterfly_new.c */
int modified_set_or(int *A1, int n1, int *A2, int n2, int *A1orA2);
// int set_and(int *A1, int n1, int *A2, int n2, int *A1andA2);
int modified_card_or(int *A1, int n1, int *A2, int n2);
// void m2s(double *mapval, double *submapval, int *subset, int count);
// void subset2map(int *A, int nA, int *subA, int nsubA);
// void s2m_sum(double *mapval, double *submapval, int *subset, int count);
// void s2m_copy(double *mapval, double *submapval, int *subset, int count);

int modified_butterfly_reduce(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax, double *val, int steps, MPI_Comm comm);

int butterfly_reshuffle(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax, double *val, int steps, MPI_Comm comm);

int butterfly_reduce_init(int *indices, int count, int **R, int *nR, int **S, int *nS, int **com_indices, int *com_count, int steps, MPI_Comm comm);

int butterfly_reshuffle_init(int *indices_in, int count_in, int *indices_out, int count_out, int **R, int *nR, int **S, int *nS, int **com_indices, int *com_count, int steps, MPI_Comm comm);

/* Content of butterfly_wrapper.c */

int mirror_butterfly(double *values_local, int *indices_local, int size_local, double *values_received, int *indices_received, int *size_received, int flag_mirror_unmirror_size_indices_data, MPI_Comm worldcomm);

int construct_butterfly_struct(Butterfly_struct *Butterfly_obj, int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, MPI_Comm worldcomm);

int prepare_butterfly_communication(int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, Butterfly_superstruct *Butterfly_superstruct, MPI_Comm worldcomm);

int perform_butterfly_communication(double *values_to_communicate, int *indices_in, int count_in, double *values_out, int *indices_out, int count_out, Butterfly_superstruct *Butterfly_superstruct_obj, MPI_Comm worldcomm);


int free_butterfly_struct(Butterfly_struct *Butterfly_obj);
int free_butterfly_supplement(Butterfly_struct_supplement *Butterfly_supp);
int free_butterfly_superstruct(Butterfly_superstruct *Butterfly_superstruct_obj);