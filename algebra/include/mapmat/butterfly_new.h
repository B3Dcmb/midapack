#ifndef BUTTERFLY_NEW_H
#define BUTTERFLY_NEW_H

#include <mpi.h>

int modified_set_or(int *A1, int n1, int *A2, int n2, int *A1orA2);
// int set_and(int *A1, int n1, int *A2, int n2, int *A1andA2);
int modified_card_or(int *A1, int n1, int *A2, int n2);
// void m2s(double *mapval, double *submapval, int *subset, int count);
// void subset2map(int *A, int nA, int *subA, int nsubA);
// void s2m_sum(double *mapval, double *submapval, int *subset, int count);
// void s2m_copy(double *mapval, double *submapval, int *subset, int count);

int modified_butterfly_reduce(int **R, int *nR, int nRmax, int **S, int *nS,
                              int nSmax, double *val, int steps, MPI_Comm comm);

int butterfly_reshuffle(int **R, int *nR, int nRmax, int **S, int *nS,
                        int nSmax, double *val, int steps, MPI_Comm comm);

int butterfly_reduce_init(int *indices, int count, int **R, int *nR, int **S,
                          int *nS, int **com_indices, int *com_count, int steps,
                          MPI_Comm comm);

int butterfly_reshuffle_init(int *indices_in, int count_in, int *indices_out,
                             int count_out, int **R, int *nR, int **S, int *nS,
                             int **com_indices, int *com_count, int steps,
                             MPI_Comm comm);

#endif // BUTTERFLY_NEW_H
