
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
// #include <mpi.h>
#include <unistd.h>
// #include "s2hat.h"
// #include <chealpix.h>

#include "midapack.h"


// int mirror_butterfly(double *values_local, int *indices_local, int size_local, double *values_received, int *indices_received, int *size_received, int flag_mirror_unmirror_size_indices_data, MPI_Comm worldcomm);

// int main_Butterfly_init(int argc, char** argv){
int main(int argc, char** argv){

    int rank, nprocs;
    MPI_Comm worldcomm;
    

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
    worldcomm = MPI_COMM_WORLD;

    int i;

    int size_local = rank + 10;
    int *indices_local;
    double *values_local;

    indices_local = (int *)malloc(size_local*sizeof(int));
    values_local = (double *)malloc(size_local*sizeof(double));

    for (i=0; i<size_local; i++){
        indices_local[i] = i+1 + rank*10;
        values_local[i] = rank+1 + i*0.01;
    }
    
    printf("%d --- ##### Data to send : \n", rank);
    printf("%d -- size : %d\n", rank, size_local);
    printf("%d -- indices_local : %d", rank, indices_local[0]);
    for (i=1; i<size_local; i++){
        printf("- %d -", indices_local[i]);
    }
    printf("\n");
    printf("%d -- values_local : %f", rank, values_local[0]);
    for (i=1; i<size_local; i++){
        printf("- %f -", values_local[i]);
    }
    printf("\n");
    fflush(stdout);
    
    
    int size_received;
    mirror_butterfly(NULL, NULL, size_local, NULL, NULL, &size_received, MIRROR_SIZE, worldcomm);
    int *indices_received = (int *)calloc(size_received,sizeof(int));
    printf("intermediate1 -- %d -- size_received : %d\n", rank, size_received);
    mirror_butterfly(NULL, indices_local, size_local, NULL, indices_received, &size_received, MIRROR_INDICES, worldcomm);
    if (size_received>0)
        printf("intermediate2 -- %d -- size_received : %d ; first indice %d \n", rank, size_received, indices_received[0]);
    double *values_received = (double *)calloc(size_received,sizeof(double));
    mirror_butterfly(values_local, NULL, size_local, values_received, NULL, &size_received, MIRROR_DATA, worldcomm);
    if (size_received>0)
        printf("intermediate3 -- %d -- size_received : %d ; first value %f \n", rank, size_received, values_received[0]);

    printf("%d --- ##### Data received : \n", rank);
    printf("%d -- size_received : %d\n", rank, size_received);
    printf("%d -- indices_received : %d", rank, indices_received[0]);
    for (i=1; i<size_received; i++){
        printf("- %d -", indices_received[i]);
    }
    printf("\n");
    printf("%d -- values_received : %f", rank, values_received[0]);
    for (i=1; i<size_received; i++){
        printf("- %f -", values_received[i]);
    }
    printf("\n");
    fflush(stdout);


    for (i=0; i<size_received; i++){
        indices_received[i] += 100;
        values_received[i] += 100;
    }

    printf("%d --- Starting send Data back !!!!!!!!!!!!!! \n", rank);
    int size_received_back;
    mirror_butterfly(NULL, NULL, size_received, NULL, NULL, &size_received_back, UNMIRROR_SIZE, worldcomm);
    printf("intermediate1U -- %d -- size_received_back : %d\n", rank, size_received_back);
    int *indices_received_back = (int *)calloc(size_received_back,sizeof(int));
    mirror_butterfly(NULL, indices_received, size_received, NULL, indices_received_back, &size_received_back, UNMIRROR_INDICES, worldcomm);
    if (size_received_back>0)
        printf("intermediate2U -- %d -- size_received_back : %d ; first indice %d \n", rank, size_received_back, indices_received_back[0]);
    double *values_received_back = (double *)calloc(size_received_back,sizeof(double));
    mirror_butterfly(values_received, NULL, size_received, values_received_back, NULL, &size_received_back, UNMIRROR_DATA, worldcomm);
    if (size_received_back>0)
        printf("intermediate3U -- %d -- size_received_back : %d ; first value %f \n", rank, size_received_back, values_received[0]);

    printf("%d --- ##### Data received_back : \n", rank);
    printf("%d -- size_back : %d\n", rank, size_received_back);
    printf("%d -- indices_received_back : %d", rank, indices_received_back[0]);
    for (i=1; i<size_received_back; i++){
        printf("- %d -", indices_received_back[i]);
    }
    printf("\n");
    printf("%d -- values_received_back : %f", rank, values_received_back[0]);
    for (i=1; i<size_received_back; i++){
        printf("- %f -", values_received_back[i]);
    }
    printf("\n");
    fflush(stdout);

    printf("%d --- ##### Data to send : \n", rank);
    printf("%d -- size : %d\n", rank, size_local);
    printf("%d -- indices_local : %d", rank, indices_local[0]);
    for (i=1; i<size_local; i++){
        printf("- %d -", indices_local[i]);
    }
    printf("\n");
    printf("%d -- values_local : %f", rank, values_local[0]);
    for (i=1; i<size_local; i++){
        printf("- %f -", values_local[i]);
    }
    printf("\n");
    fflush(stdout);

    printf("%d --- END TEST 0 \n", rank); fflush(stdout);
    free(indices_local);
    free(values_local);
    free(indices_received);
    printf("%d --- END TEST 1 \n", rank); fflush(stdout);
    free(values_received);
    free(indices_received_back);
    free(values_received_back);

    return 0;
}


// int construct_butterfly_struct(Butterfly_struct *Butterfly_obj, int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, MPI_Comm worldcomm);

// int prepare_butterfly_communication(int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, Butterfly_superstruct *Butterfly_superstruct, MPI_Comm worldcomm);

// int perform_butterfly_communication(double *values_to_communicate, int *indices_in, int count_in, double *values_out, int *indices_out, int count_out, Butterfly_superstruct *Butterfly_superstruct_obj, MPI_Comm worldcomm);


// int free_butterfly_struct(Butterfly_struct *Butterfly_obj);
// int free_butterfly_supplement(Butterfly_struct_supplement *Butterfly_supp);
// int free_butterfly_superstruct(Butterfly_superstruct *Butterfly_superstruct_obj);
