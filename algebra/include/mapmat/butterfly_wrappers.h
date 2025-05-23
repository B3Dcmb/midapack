#ifndef BUTTERFLY_WRAPPERS_H
#define BUTTERFLY_WRAPPERS_H

#include <mpi.h>

#define MIRROR_SIZE 0
#define MIRROR_INDICES 1
#define MIRROR_DATA 2
#define UNMIRROR_SIZE 3
#define UNMIRROR_INDICES 4
#define UNMIRROR_DATA 5

typedef struct Butterfly_struct {
    int *com_indices, com_count; // communicated indices, and size
    int steps;                   // number of steps in the butterfly scheme
    int *nS, *nR; // number of indices (to send and to receive); size = steps
    int **R, **S; // sending or receiving indices

    int classic_or_reshuffle_butterfly;
    // Flag to indicate if classic or reshuffled butterfly is used, e.g. if we
    // expect the (pixel) distributions prior and after communication to be the
    // same or different 0 for classic butterfly scheme with same pixel
    // distributions prior/after ; 1 for after

    MPI_Comm comm_butterfly;
    // MPI communicator used by the butterfly scheme
} Butterfly_struct;

typedef struct Butterfly_struct_supplement {
    // Indices obtained from mirroring
    int *indices_mirror;
    // List of indices received in the mirroring step by the local MPI process
    // -> Will be used afterwards to resend back the relevant indices when
    // unmirroring
    int size_from_mirror; // Size of the indices obtained from mirroring

    int new_size_local; // New size after mirroring when the redundant indices
                        // have been deleted
    int *ordered_indices;
    // Ordered indices used for butterfly scheme, obtained from both mirroring
    // and any conversion of the (pixel) distribution This MUST BE ORDERED
    // INDICES
} Butterfly_struct_supplement;

typedef struct Butterfly_superstruct {
    Butterfly_struct Butterfly_obj;

    Butterfly_struct_supplement Butterfly_mirror_supp;
    Butterfly_struct_supplement Butterfly_unmirror_supp;
} Butterfly_superstruct;

int mirror_butterfly(double *values_local, int *indices_local, int size_local,
                     double *values_received, int *indices_received,
                     int *size_received,
                     int flag_mirror_unmirror_size_indices_data,
                     MPI_Comm worldcomm);

int construct_butterfly_struct(Butterfly_struct *Butterfly_obj, int *indices_in,
                               int count_in, int *indices_out, int count_out,
                               int flag_classic_or_reshuffle_butterfly,
                               MPI_Comm worldcomm);

int prepare_butterfly_communication(
    int *indices_in, int count_in, int *indices_out, int count_out,
    int flag_classic_or_reshuffle_butterfly,
    Butterfly_superstruct *Butterfly_superstruct, MPI_Comm worldcomm);

int perform_butterfly_communication(
    double *values_to_communicate, int *indices_in, int count_in,
    double *values_out, int *indices_out, int count_out,
    Butterfly_superstruct *Butterfly_superstruct_obj, MPI_Comm worldcomm);

int free_butterfly_struct(Butterfly_struct *Butterfly_obj, int rank);
int free_butterfly_supplement(Butterfly_struct_supplement *Butterfly_supp,
                              int rank);
int free_butterfly_superstruct(Butterfly_superstruct *Butterfly_superstruct_obj,
                               int rank);

#endif // BUTTERFLY_WRAPPERS_H
