
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

int main_Butterfly_mirorring(int argc, char** argv){
// int main(int argc, char** argv){

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
        indices_local[i] = i+1 + rank*100;
        values_local[i] = rank + (i+1)*0.01;
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

// int butterfly_reshuffle(int **R, int *nR, int nRmax, int **S, int *nS, int nSmax, double *val, int steps, MPI_Comm comm);

// int free_butterfly_struct(Butterfly_struct *Butterfly_obj);



int main_Butterfly_scheme_tests_classic(int argc, char** argv){ // WORKS ONLY WITH POWER OF 2 NUMBER OF PROCESSES !!!!!
// int main(int argc, char** argv){

    int rank, nprocs;
    MPI_Comm worldcomm;
    

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
    worldcomm = MPI_COMM_WORLD;

    int i, k;

    int size_local = rank + 10;
    int *indices_local;
    double *values_local;

    indices_local = (int *)malloc(size_local*sizeof(int));
    values_local = (double *)malloc(size_local*sizeof(double));

    for (i=0; i<size_local; i++){
        indices_local[i] = i+1 + rank*100;
        values_local[i] = rank + (i+1)*0.01;
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
    
    // Butterfly_struct *Butterfly_obj_classic = (Butterfly_struct *)malloc(1*sizeof(Butterfly_struct));

    // int rank_to_receive = nprocs - rank - 1;
    // int rank_to_receive = fabs((rank+1)%nprocs);
    int rank_to_receive = rank;
    int size_to_receive = (rank_to_receive + 10);
    // int size_to_receive = (rank + 10)/2 + (rank_to_receive + 10)/2;
    int *indices_to_receive;
    double *values_to_receive;

    // if (rank == nprocs -2)
    //     size_to_receive = 0;

    indices_to_receive = (int *)malloc(size_to_receive*sizeof(int));
    values_to_receive = (double *)malloc(size_to_receive*sizeof(double));

    // int rank_min = min(rank, rank_to_receive);
    // int rank_max = max(rank, rank_to_receive);

    for (i=0; i<size_to_receive ; i++){
        indices_to_receive[i] = i+1 + rank*100;
    }
    // for (i=0; i<(rank_min + 10)/2 ; i++){
    //     indices_to_receive[i] = i+1 + rank_min*100;
    // }for (i=(rank_min + 10)/2; i<size_to_receive; i++){
    //     indices_to_receive[i] = i+1 + rank_max*100;
    // }

    int steps = log_2(nprocs);
    printf("%d --- Constructing butterfly_struct -- %d ; steps %d \n", rank, size_to_receive, steps); fflush(stdout);
    int flag_classic_or_reshuffle_butterfly = 0; // 1;
    // construct_butterfly_struct(Butterfly_obj_classic, indices_local, size_local, indices_to_receive, size_to_receive, flag_classic_or_reshuffle_butterfly, worldcomm);
    
    int **S = (int **)malloc(steps * sizeof(int *)); // allocate sending maps tab
    int **R = (int **)malloc(steps * sizeof(int *)); // allocate receiving maps tab
    int *nS = (int *)malloc(steps * sizeof(int));   // allocate sending map sizes tab
    int *nR = (int *)malloc(steps * sizeof(int));   // allocate receiving map size tab
    int *com_indices;
    int com_count;
    butterfly_init(indices_local, size_local, R, nR, S, nS, &com_indices, &com_count, steps, worldcomm);
    // printf("%d »»»»»» Results post init Butterfly !! %d \n", rank, Butterfly_obj_classic->com_count); fflush(stdout);
    // if (Butterfly_obj_classic->com_count){
    //     printf("%d »»»»»» com_indices : %d", rank, Butterfly_obj_classic->com_indices[0]);
    //     for (k=1; k<Butterfly_obj_classic->com_count; k++){
    //         printf("- %d -",Butterfly_obj_classic->com_indices[k]);
    //     }
    //     printf("\n"); fflush(stdout);
    // }

    int nRmax, nSmax;
    printf("%d --- Getting nRmax, nSmax \n", rank); fflush(stdout);
    if (size_to_receive>0){
        nRmax = nR[0];
        nSmax = nS[0];
        for (k = 1; k < steps; k++){
            if (nR[k] > nRmax)
                nRmax = nR[k];
            if (nS[k] > nSmax)
                nSmax = nS[k];
        }
    }

    printf("%d --- Defining com_val -- %d \n", rank, com_count); fflush(stdout);
    double *com_val = (double *)calloc(com_count, sizeof(double));
    // for (k = 0; k < Butterfly_obj_classic->com_count; k++)
    //     com_val[k] = 0.0;

    printf("%d --- Retrieving values \n", rank); fflush(stdout);
    m2m(values_local, indices_local, size_local, com_val, com_indices, com_count);

    butterfly_reduce(R, nR, nRmax, S, nS, nSmax, com_val, steps, worldcomm);
    printf("%d --- Performing butterfly : nRmax %d ; nSmax %d \n", rank, nRmax, nSmax); fflush(stdout);
    // switch(Butterfly_obj_classic->classic_or_reshuffle_butterfly)
    //     {
    //         case 0: // Classic butterfly
    //             modified_butterfly_reduce(Butterfly_obj_classic->R, Butterfly_obj_classic->nR, nRmax, Butterfly_obj_classic->S, Butterfly_obj_classic->nS, nSmax, com_val, Butterfly_obj_classic->steps, Butterfly_obj_classic->comm_butterfly);
    //             break;

    //         case 1: // Reshuffle butterfly
    //             butterfly_reshuffle(Butterfly_obj_classic->R, Butterfly_obj_classic->nR, nRmax, Butterfly_obj_classic->S, Butterfly_obj_classic->nS, nSmax, com_val, Butterfly_obj_classic->steps, Butterfly_obj_classic->comm_butterfly);
    //             break;
    //     }
    
    // printf("%d <<<<<< Results from Butterfly !! %d \n", rank, Butterfly_obj_classic->com_count); fflush(stdout);
    // if (Butterfly_obj_classic->com_count){
    //     printf("%d <<<<<< com_indices com_val : %d %f", rank, Butterfly_obj_classic->com_indices[0], com_val[0]);
    //     for (k=1; k<Butterfly_obj_classic->com_count; k++){
    //         printf("- %d %f -",Butterfly_obj_classic->com_indices[k], com_val[k]);
    //     }
    //     printf("\n"); fflush(stdout);
    // }

    printf("%d --- Performing butterfly2 \n", rank); fflush(stdout);
    // values_to_receive = (double *)malloc(size_to_receive*sizeof(double));
    m2m(com_val, com_indices, com_count, values_to_receive, indices_to_receive, size_to_receive);
    free(com_val);

    printf("%d --- ##### Data received : \n", rank); fflush(stdout);
    printf("%d -- size_received : %d\n", rank, size_to_receive);
    if (size_to_receive > 0){
        printf("%d -- indices_received : %d", rank, indices_to_receive[0]);
        for (i=1; i<size_to_receive; i++){
            printf("- %d -", indices_to_receive[i]);
        }
        printf("\n");
        printf("%d -- values_to_receive : %f", rank, values_to_receive[0]);
        for (i=1; i<size_to_receive; i++){
            printf("- %f -", values_to_receive[i]);
        }
        printf("\n");
    }
    fflush(stdout);

    // free_butterfly_struct(Butterfly_obj_classic, rank);
    printf("%d - Done ! \n", rank);

    return 0;
}

int main_Butterfly_scheme_tests_v0(int argc, char** argv){ // WORKS ONLY WITH POWER OF 2 NUMBER OF PROCESSES !!!!!
// int main(int argc, char** argv){

    int rank, nprocs;
    MPI_Comm worldcomm;
    

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
    worldcomm = MPI_COMM_WORLD;

    int i, k;

    int size_local = rank + 10;
    int *indices_local;
    double *values_local;

    indices_local = (int *)malloc(size_local*sizeof(int));
    values_local = (double *)malloc(size_local*sizeof(double));

    for (i=0; i<size_local; i++){
        indices_local[i] = i+1 + rank*100;
        values_local[i] = rank + (i+1)*0.01;
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
    
    Butterfly_struct *Butterfly_obj_classic = (Butterfly_struct *)malloc(1*sizeof(Butterfly_struct));

    // int rank_to_receive = nprocs - rank - 1;
    int rank_to_receive = fabs((rank+1)%nprocs);
    // int rank_to_receive = rank;
    // int size_to_receive = (rank_to_receive + 10)/2;
    int size_to_receive = (rank + 10)/2 + (rank_to_receive + 10)/2;
    int *indices_to_receive;
    double *values_to_receive;

    // if (rank == nprocs -2)
    //     size_to_receive = 0;

    indices_to_receive = (int *)malloc(size_to_receive*sizeof(int));
    values_to_receive = (double *)malloc(size_to_receive*sizeof(double));

    int rank_min = min(rank, rank_to_receive);
    int rank_max = max(rank, rank_to_receive);

    for (i=0; i<(rank_min + 10)/2 ; i++){
        indices_to_receive[i] = i+1 + rank_min*100;
    }for (i=(rank_min + 10)/2; i<size_to_receive; i++){
        indices_to_receive[i] = i+1 + rank_max*100;
    }

    printf("%d --- Constructing butterfly_struct -- %d \n", rank, size_to_receive); fflush(stdout);
    int flag_classic_or_reshuffle_butterfly = 1; // 1;
    construct_butterfly_struct(Butterfly_obj_classic, indices_local, size_local, indices_to_receive, size_to_receive, flag_classic_or_reshuffle_butterfly, worldcomm);
    
    printf("%d »»»»»» Results post init Butterfly !! %d \n", rank, Butterfly_obj_classic->com_count); fflush(stdout);
    if (Butterfly_obj_classic->com_count){
        printf("%d »»»»»» com_indices : %d", rank, Butterfly_obj_classic->com_indices[0]);
        for (k=1; k<Butterfly_obj_classic->com_count; k++){
            printf("- %d -",Butterfly_obj_classic->com_indices[k]);
        }
        printf("\n"); fflush(stdout);
    }

    int nRmax, nSmax;
    printf("%d --- Getting nRmax, nSmax \n", rank); fflush(stdout);
    if (size_to_receive>0){
        nRmax = Butterfly_obj_classic->nR[0];
        nSmax = Butterfly_obj_classic->nS[0];
        for (k = 1; k < Butterfly_obj_classic->steps; k++){
            if (Butterfly_obj_classic->nR[k] > nRmax)
                nRmax = Butterfly_obj_classic->nR[k];
            if (Butterfly_obj_classic->nS[k] > nSmax)
                nSmax = Butterfly_obj_classic->nS[k];
        }
    }

    printf("%d --- Defining com_val \n", rank); fflush(stdout);
    double *com_val = (double *)calloc(Butterfly_obj_classic->com_count, sizeof(double));
    // for (k = 0; k < Butterfly_obj_classic->com_count; k++)
    //     com_val[k] = 0.0;

    printf("%d --- Retrieving values \n", rank); fflush(stdout);
    m2m(values_local, indices_local, size_local, com_val, Butterfly_obj_classic->com_indices, Butterfly_obj_classic->com_count);

    printf("%d --- Performing butterfly : nRmax %d ; nSmax %d \n", rank, nRmax, nSmax); fflush(stdout);
    switch(Butterfly_obj_classic->classic_or_reshuffle_butterfly)
        {
            case 0: // Classic butterfly
                modified_butterfly_reduce(Butterfly_obj_classic->R, Butterfly_obj_classic->nR, nRmax, Butterfly_obj_classic->S, Butterfly_obj_classic->nS, nSmax, com_val, Butterfly_obj_classic->steps, Butterfly_obj_classic->comm_butterfly);
                break;

            case 1: // Reshuffle butterfly
                butterfly_reshuffle(Butterfly_obj_classic->R, Butterfly_obj_classic->nR, nRmax, Butterfly_obj_classic->S, Butterfly_obj_classic->nS, nSmax, com_val, Butterfly_obj_classic->steps, Butterfly_obj_classic->comm_butterfly);
                break;
        }
    
    printf("%d <<<<<< Results from Butterfly !! %d \n", rank, Butterfly_obj_classic->com_count); fflush(stdout);
    if (Butterfly_obj_classic->com_count){
        printf("%d <<<<<< com_indices com_val : %d %f", rank, Butterfly_obj_classic->com_indices[0], com_val[0]);
        for (k=1; k<Butterfly_obj_classic->com_count; k++){
            printf("- %d %f -",Butterfly_obj_classic->com_indices[k], com_val[k]);
        }
        printf("\n"); fflush(stdout);
    }

    printf("%d --- Performing butterfly2 \n", rank); fflush(stdout);
    values_to_receive = (double *)malloc(size_to_receive*sizeof(double));
    m2m(com_val, Butterfly_obj_classic->com_indices, Butterfly_obj_classic->com_count, values_to_receive, indices_to_receive, size_to_receive);
    free(com_val);

    printf("%d --- ##### Data received : \n", rank); fflush(stdout);
    printf("%d -- size_received : %d\n", rank, size_to_receive);
    if (size_to_receive > 0){
        printf("%d -- indices_received : %d", rank, indices_to_receive[0]);
        for (i=1; i<size_to_receive; i++){
            printf("- %d -", indices_to_receive[i]);
        }
        printf("\n");
        printf("%d -- values_to_receive : %f", rank, values_to_receive[0]);
        for (i=1; i<size_to_receive; i++){
            printf("- %f -", values_to_receive[i]);
        }
        printf("\n");
    }
    fflush(stdout);

    free_butterfly_struct(Butterfly_obj_classic, rank);
    printf("%d - Done ! \n", rank);

    return 0;
}


// int prepare_butterfly_communication(int *indices_in, int count_in, int *indices_out, int count_out, int flag_classic_or_reshuffle_butterfly, Butterfly_superstruct *Butterfly_superstruct, MPI_Comm worldcomm);

// int perform_butterfly_communication(double *values_to_communicate, int *indices_in, int count_in, double *values_out, int *indices_out, int count_out, Butterfly_superstruct *Butterfly_superstruct_obj, MPI_Comm worldcomm);

// int free_butterfly_superstruct(Butterfly_superstruct *Butterfly_superstruct_obj);


// int main_Butterfly_scheme_tests_vfinale(int argc, char** argv){
int main(int argc, char** argv){

    int rank, nprocs;
    MPI_Comm worldcomm;
    

    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
    worldcomm = MPI_COMM_WORLD;

    int i, k;

    int size_local = rank + 10;
    int *indices_local;
    double *values_local;

    indices_local = (int *)malloc(size_local*sizeof(int));
    values_local = (double *)malloc(size_local*sizeof(double));

    for (i=0; i<size_local; i++){
        indices_local[i] = i+1 + rank*100;
        values_local[i] = rank + (i+1)*0.01;
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
    
    // Butterfly_struct *Butterfly_obj_classic = (Butterfly_struct *)malloc(1*sizeof(Butterfly_struct));
    // Butterfly_superstruct *Butterfly_superstruct_obj = (Butterfly_superstruct *)malloc(1*sizeof(Butterfly_superstruct));
    Butterfly_superstruct Butterfly_superstruct_obj_;
    // int rank_to_receive = nprocs - rank - 1;
    // int rank_to_receive = fabs((rank-1)%nprocs);
    int rank_to_receive = (rank+1)%nprocs;
    // int rank_to_receive = rank;
    // int size_to_receive = (rank_to_receive + 10)/2;
    int size_to_receive = (rank_to_receive + 10);
    int *indices_to_receive;
    double *values_to_receive;

    // if (rank == nprocs-2){
    //     size_to_receive = 0;
    // }

    indices_to_receive = (int *)malloc(size_to_receive*sizeof(int));
    values_to_receive = (double *)malloc(size_to_receive*sizeof(double));

    for (i=0; i<size_to_receive; i++){
        indices_to_receive[i] = i+1 + rank_to_receive*100;
    }


    int flag_classic_or_reshuffle_butterfly = 1;
    // int flag_classic_or_reshuffle_butterfly = 0;
    // construct_butterfly_struct(Butterfly_obj_classic, indices_local, size_local, indices_to_receive, size_to_receive, flag_classic_or_reshuffle_butterfly, worldcomm);
    printf("%d --- Preparing butterfly communication \n", rank); fflush(stdout);
    prepare_butterfly_communication(indices_local, size_local, indices_to_receive, size_to_receive, flag_classic_or_reshuffle_butterfly, &Butterfly_superstruct_obj_, worldcomm);
    Butterfly_superstruct *Butterfly_superstruct_obj = &Butterfly_superstruct_obj_;

    int number_steps = log_2(nprocs);
    int nb_butterfly_ranks = pow_2(number_steps);

    Butterfly_struct *Butterfly_obj = &(Butterfly_superstruct_obj->Butterfly_obj);
    Butterfly_struct_supplement *Butterfly_mirror_supp = &(Butterfly_superstruct_obj->Butterfly_mirror_supp);
    Butterfly_struct_supplement *Butterfly_unmirror_supp = &(Butterfly_superstruct_obj->Butterfly_unmirror_supp);

    if (rank == nb_butterfly_ranks-1){

        printf("%d ##### Butterfly_supplement tests mirror !! nb_butterfly_ranks %d ; nb_steps %d ; new_size_local : %d \n", rank, nb_butterfly_ranks, number_steps, Butterfly_mirror_supp->new_size_local);
        printf("%d ##### order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
        for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
            printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
        }
        printf("\n");
        printf("%d ##### Butterfly_supplement tests !! nb_butterfly_ranks %d ; nb_steps %d ; new_size_local : %d \n", rank, nb_butterfly_ranks, number_steps, Butterfly_unmirror_supp->new_size_local);
        printf("%d ##### order_indices unmirror : %d -", rank, Butterfly_unmirror_supp->ordered_indices[0]);
        for (i=1; i<Butterfly_unmirror_supp->new_size_local; i++){
            printf("- %d -", Butterfly_unmirror_supp->ordered_indices[i]);
        }
        printf("\n");
        
        printf("%d ##### Prep Butterfly !! %d \n", rank, Butterfly_obj->com_count);
        if (Butterfly_obj->com_count){
            printf("%d ##### com_indices : %d -", rank, Butterfly_obj->com_indices[0]);
            for (k=1; k<Butterfly_obj->com_count; k++){
                printf("- %d -",Butterfly_obj->com_indices[k]);
            }
        }
        printf("\n");
    }


    printf("%d #### Test new_size_local %d \n", rank, Butterfly_mirror_supp->new_size_local); fflush(stdout);
    int *test_ordered_indices = (int *)malloc(Butterfly_mirror_supp->new_size_local*sizeof(int));
    test_ordered_indices = memcpy(test_ordered_indices, Butterfly_mirror_supp->ordered_indices, Butterfly_mirror_supp->new_size_local*sizeof(int));
    // for (i=0; i<Butterfly_mirror_supp->new_size_local; i++)
    //     test_ordered_indices[i] = Butterfly_mirror_supp->ordered_indices[i];
    
    // Butterfly_struct_supplement *Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
    // Butterfly_struct_supplement *Butterfly_unmirror_supp = Butterfly_superstruct_obj->Butterfly_unmirror_supp;

    if (rank == 1){
        if (Butterfly_mirror_supp->new_size_local){
            // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
            printf("%d #####0 Butterfly_supplement tests mirror !! nb_butterfly_ranks %d ; nb_steps %d ; new_size_local : %d \n", rank, nb_butterfly_ranks, number_steps, Butterfly_mirror_supp->new_size_local);
            printf("%d #####0 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
            for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
                printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
            }
            printf("\n"); fflush(stdout);
        }
    }

    printf("%d --- Performing butterfly communication \n", rank); fflush(stdout);
    perform_butterfly_communication(values_local, indices_local, size_local, values_to_receive, indices_to_receive, size_to_receive, Butterfly_superstruct_obj, worldcomm);

    if (rank == 1){
        if (Butterfly_mirror_supp->new_size_local){
            // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
            printf("%d #####3 Butterfly_supplement tests mirror !! nb_butterfly_ranks %d ; nb_steps %d ; new_size_local : %d \n", rank, nb_butterfly_ranks, number_steps, Butterfly_mirror_supp->new_size_local);
            printf("%d #####3 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
            for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
                printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
            }
            printf("\n"); fflush(stdout);
        }
    }
    printf("%d --- ##### Data received : \n", rank); fflush(stdout);
    printf("%d -- size_received : %d\n", rank, size_to_receive);
    if (size_to_receive > 0){
        printf("%d -- indices_received : %d", rank, indices_to_receive[0]);
        for (i=1; i<size_to_receive; i++){
            printf("- %d -", indices_to_receive[i]);
        }
        printf("\n");
        printf("%d -- values_to_receive : %f", rank, values_to_receive[0]);
        for (i=1; i<size_to_receive; i++){
            printf("- %f -", values_to_receive[i]); fflush(stdout);
        }
        printf("\n"); fflush(stdout);
    }
    fflush(stdout);

    // if (rank == nb_butterfly_ranks-1){

    printf("%d #####2 Butterfly_supplement tests mirror !! nb_butterfly_ranks %d ; nb_steps %d ; new_size_local : %d \n", rank, nb_butterfly_ranks, number_steps, Butterfly_mirror_supp->new_size_local);
    printf("%d #####2 order_indices mirror : %d -", rank, Butterfly_mirror_supp->ordered_indices[0]);
    for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
        printf("- %d -", Butterfly_mirror_supp->ordered_indices[i]);
    }
    printf("\n");
    printf("%d #####2 Butterfly_supplement tests !! nb_butterfly_ranks %d ; nb_steps %d ; new_size_local : %d \n", rank, nb_butterfly_ranks, number_steps, Butterfly_unmirror_supp->new_size_local);
    printf("%d #####2 order_indices unmirror : %d -", rank, Butterfly_unmirror_supp->ordered_indices[0]);
    for (i=1; i<Butterfly_unmirror_supp->new_size_local; i++){
        printf("- %d -", Butterfly_unmirror_supp->ordered_indices[i]);
    }
    printf("\n");
    // Butterfly_struct *Butterfly_obj = Butterfly_superstruct_obj->Butterfly_obj;
    // printf("%d #####2 Prep Butterfly !! %d \n", rank, Butterfly_obj->com_count);
    // if (Butterfly_obj->com_count){
    //     printf("%d #####2 com_indices : %d -", rank, Butterfly_obj->com_indices[0]);
    //     for (k=1; k<Butterfly_obj->com_count; k++){
    //         printf("- %d -",Butterfly_obj->com_indices[k]);
    //     }
    //     printf("\n");
    // }


    printf("%d - Free 1 \n", rank); fflush(stdout);
    // free_butterfly_supplement(Butterfly_superstruct_obj->Butterfly_mirror_supp, rank);
    printf("%d - Free 2 \n", rank); fflush(stdout);
    // free_butterfly_supplement(Butterfly_superstruct_obj->Butterfly_unmirror_supp, rank);
    printf("%d - Free 0 \n", rank); fflush(stdout);
    // free_butterfly_struct(Butterfly_superstruct_obj->Butterfly_obj, rank);
    printf("%d - Free 3 \n", rank); fflush(stdout);
    // if (Butterfly_mirror_supp->size_from_mirror)
    //     free(Butterfly_superstruct_obj->Butterfly_mirror_supp->indices_mirror);
    printf("%d - Free 4 \n", rank); fflush(stdout);

    if (Butterfly_mirror_supp->new_size_local){
        // Butterfly_mirror_supp = Butterfly_superstruct_obj->Butterfly_mirror_supp;
        printf("%d #####F Butterfly_supplement tests mirror !! nb_butterfly_ranks %d ; nb_steps %d ; new_size_local : %d \n", rank, nb_butterfly_ranks, number_steps, Butterfly_mirror_supp->new_size_local); fflush(stdout);
        printf("%d #####F order_indices mirror : %d %d -", rank, Butterfly_mirror_supp->ordered_indices[0], test_ordered_indices[0]); fflush(stdout);
        for (i=1; i<Butterfly_mirror_supp->new_size_local; i++){
            printf("- %d %d -", Butterfly_mirror_supp->ordered_indices[i], test_ordered_indices[i]);
        }
        printf("\n"); fflush(stdout);
        printf("%d #####F Butterfly_supplement tests indices_mirror !! nb_butterfly_ranks %d ; nb_steps %d ; size_from_mirror : %d \n", rank, nb_butterfly_ranks, number_steps, Butterfly_mirror_supp->size_from_mirror); fflush(stdout);
        if (Butterfly_mirror_supp->size_from_mirror){
            printf("%d #####F indices_mirror : %d -", rank, Butterfly_mirror_supp->indices_mirror[0]); fflush(stdout);
            for (i=1; i<Butterfly_mirror_supp->size_from_mirror; i++){
                printf("- %d -", Butterfly_mirror_supp->indices_mirror[i]);
            }
            printf("\n"); fflush(stdout);
        }
        // int *ordered_indices__ = Butterfly_mirror_supp->ordered_indices;
        // free(ordered_indices__);
    }

    printf("%d #### Test2 new_size_local %d \n", rank, Butterfly_mirror_supp->new_size_local); fflush(stdout);
    printf("%d --- ##### Free step \n", rank); fflush(stdout);
    free_butterfly_superstruct(Butterfly_superstruct_obj, rank);

    printf("%d - Free 5 \n", rank); fflush(stdout);
    // free(Butterfly_superstruct_obj->Butterfly_mirror_supp);
    
    free(indices_local);
    free(values_local);
    free(indices_to_receive);
    free(values_to_receive);
    
    printf("%d --- END TEST 1 \n", rank); fflush(stdout);
    free(test_ordered_indices);
    

    printf("%d - Done ! \n", rank);

    return 0;
}
