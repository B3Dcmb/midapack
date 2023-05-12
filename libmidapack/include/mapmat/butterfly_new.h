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
    int *indices_from_mirror; // List of indices received in the mirroring step by the local MPI process  -> Will be used afterwards to resend back the relevant indices when unmirroring
    int size_from_mirror; // Size of the indices obtained from mirroring

    // int do_we_need_to_project_into_different_scheme;
    // Flag to indicate if we need to project the input pixel distribution into a different scheme
    
    int *projector_values;
    // Projector for local values of map : allows to project the values from one (pixel) distribution to another
    // In the context of harmonic transformations, initialized for conversion between the nest pixel distribution of MAPPRAISER and ring pixel distribution of S2HAT
    
    int new_size_local; // New size after mirroring when the redundant indices have been deleted
    int *ordered_indices; // Indices ring used for butterfly scheme, obtained from both mirroring and any conversion of the (pixel) distribution
    // This MUST BE ORDERED INDICES
    // In the context of harmonic operations, it is used to store ordered indices, used for ring2nest and nest2ring transitions in which case it corresponds to ordered ring indices
} Butterfly_struct;


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

int mirror_butterfly(double *values_local, int *indices_local, int size_local, double *values_received, int *indices_received, int *size_received, int flag_mirror_unmirror_size_indices_data, MPI_Comm world_comm);
int construct_butterfly_struct(Butterfly_struct *Butterfly_obj, int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, int do_we_need_to_project_into_different_scheme, MPI_Comm worlcomm);
