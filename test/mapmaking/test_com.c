// Midapack library
// Test the different communication schemes implmented in Midapack
// by performing the product At*V.
/** @file   test_com.c
    @author Hamza El Bouhargani
    @date   October 2018 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "midapack.h"


void usage(){
    printf("usage...\n");
}

//NERSC
extern const char *WORKDIR="/global/cscratch1/sd/elbouha/data2/";

int main(int argc, char *argv[])
{

  int64_t	M;       //Global number of rows
  int 		N, Nnz;  //of columns, of non-zeros values per column for the pointing matrix A
  int		m, n;  //local number of rows, of columns for the pointing matrix A
  int64_t	gif;			//global indice for the first local line
  int		i, j, k;
  Mat	A;			        //pointing matrix structure
  int 		*indices;
  double 	*values;
  int 		pointing_commflag ;	//option for the communication scheme for the pointing matrix
  double	*b, *Ag, *Ad, *pol_ang; 	 	//temporal domain vectors
  double *v; //map domain vector
  double	t0, t1, st, t;		 	//timer, start time
  int 		rank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("rank = %d, size = %d",rank,size);

//global data caracteristics
  int Nb_t_Intervals = 128;//1352;//128;//2;//256;//8;           //total number of stationnary intervals
  int t_Interval_length = 17899900;//1431992;//139992; //1431992;//2863984;//1431992;//pow(2,25);//pow(2,25);          //length for each stationnary interval
  int t_Interval_length_true = 17899900;//1431992;//139992;//1431992;//2863984;//1431992;//pow(2,20);
  Nnz=3;

//Number of loop we need to read the all t_Interval_length
  int t_Interval_loop = t_Interval_length/t_Interval_length_true ;
  printf("[rank %d] t_Interval_loop=%d\n", rank, t_Interval_loop );

//total length of the time domaine signal
  M = (int64_t) Nb_t_Intervals*t_Interval_length ;

  printf("[rank %d] M=%ld\n", rank, M);

//compute distribution indexes over the processes
  partition(&gif, &m, M, rank, size);

  double Nb_t_Intervals_loc_dble = Nb_t_Intervals/size;
  int Nb_t_Intervals_loc = ceil( Nb_t_Intervals*1.0/size );
  int nb_proc_shared_one_interval = max(1, size/Nb_t_Intervals ); //same as ceil( (size*1.0)/Nb_t_Intervals );
  printf("[rank %d] nb_proc_shared_one_interval=%d\n", rank, nb_proc_shared_one_interval );

  int t_Interval_length_loc = t_Interval_length/nb_proc_shared_one_interval; //just to check
//should be equal to min(m ,t_Interval_length)

  printf("[rank %d] size=%d \t m=%d \t Nb_t_Intervals=%d \t t_Interval_length=%d\n", rank, size, m, Nb_t_Intervals, t_Interval_length );
  printf("[rank %d] Nb_t_Intervals_loc=%d \t t_Interval_length_loc=%d\n", rank, Nb_t_Intervals_loc , t_Interval_length_loc);


//input data memory allocation
  indices  = (int *) malloc(Nnz*m * sizeof(int));     //for pointing matrix indices
  values  = (double *) malloc(Nnz*m * sizeof(double));//for pointing matrix values
  b   = (double *) malloc(m*sizeof(double));    //full raw data vector for the signal
  pol_ang = (double *) malloc(m * sizeof(double));

  //Read data from files:
//note: work only if the number of processes is a multiple of the number of stationary intervals

//Definition for the input data
  int part_id;      // stationnaly period id number
  int *point_data;  // scann strategy input data for the pointing matrix
  double *signal;   // signal input data
  double *polar; // linear polarization angles



//Init the pointing matrix values to unity
  for(i=0; i<(Nnz*m); i++)
    values[i] = 1.;
  int number_in_interval;

int flag_bigdata=1;

if (flag_bigdata==1 && nb_proc_shared_one_interval>1) {

if (rank==0) {
  printf("#######  ENTER BiG DATA MODE   ################\n");
}
  part_id = rank/nb_proc_shared_one_interval; //floor(rank/nb_proc_shared_one_interval);
  number_in_interval = rank%nb_proc_shared_one_interval;
  printf("[rank %d] part_id=%d\n", rank, part_id );  //interval id number
  printf("[rank %d] number_in_interval=%d\n", rank, number_in_interval );

  int t_Interval_loop_loc = ceil( t_Interval_loop*1.0/nb_proc_shared_one_interval);
  printf("[rank %d] t_Interval_loop_loc=%d\n", rank, t_Interval_loop_loc );
  printf("[rank %d] m=%d \t t_Interval_length_true*t_Interval_loop_loc=%d\n", rank, m, t_Interval_length_true*t_Interval_loop_loc );

  point_data  = (int *) malloc(Nnz*t_Interval_length_true*t_Interval_loop_loc * sizeof(int));
  signal      = (double *) malloc(t_Interval_length_true*t_Interval_loop_loc  * sizeof(double));
  polar = (double *) malloc(t_Interval_length_true*t_Interval_loop_loc * sizeof(double));

  for (i=0; i < t_Interval_loop_loc; ++i) {
    ioReadfile_pol(t_Interval_length, part_id, point_data+t_Interval_length_true*Nnz*i, signal+t_Interval_length_true*i, polar+t_Interval_length_true*i);
  }
//just keep the relevant part of the stationary interval for the local process
  int nb_proc_shared_one_subinterval = max(1, size/(Nb_t_Intervals*t_Interval_loop) );
  //same as ceil( (size*1.0)/Nb_t_Intervals );
  int number_in_subinterval = rank%nb_proc_shared_one_subinterval;

  int t_Interval_length_subinterval_loc = t_Interval_length_true/nb_proc_shared_one_subinterval;
//note: we must have this to be exactly an integer.

  printf("[rank %d] nb_proc_shared_one_subinterval=%d\n", rank, nb_proc_shared_one_subinterval );
  printf("[rank %d] number_in_subinterval=%d\n", rank, number_in_subinterval );


  for (i=0; i<(Nnz*m); i++)
    indices[i]=point_data[i+Nnz*number_in_subinterval*t_Interval_length_subinterval_loc];

  for(i=0; i<(m); i++){
    b[i] = signal[i+number_in_subinterval*t_Interval_length_subinterval_loc];
    pol_ang[i] = polar[i+ number_in_subinterval*t_Interval_length_subinterval_loc];
  }
  for(i=0; i<(Nnz*m);i+=3){
    values[i+1] = cos(2*pol_ang[(int)i/3]);
    values[i+2] = sin(2*pol_ang[(int)i/3]);
  }

  free(point_data);
  free(signal);
  free(polar);
  free(pol_ang);

}
else {

//for the case we share the stationary intervals in severals processes
if (nb_proc_shared_one_interval>1) {
  part_id = rank/nb_proc_shared_one_interval; //floor(rank/nb_proc_shared_one_interval);
  number_in_interval = rank%nb_proc_shared_one_interval;
  printf("[rank %d] part_id=%d\n", rank, part_id );  //interval id number
  printf("[rank %d] number_in_interval=%d\n", rank, number_in_interval );

  point_data = (int *) malloc(Nnz*t_Interval_length * sizeof(int));
  signal     = (double *) malloc(t_Interval_length * sizeof(double));

//read the all stationary interval
  for (i=0; i < t_Interval_loop; ++i) {
    ioReadfile(t_Interval_length, part_id, point_data+t_Interval_length_true*Nnz*i, signal+t_Interval_length_true*i);
  }
//  ioReadfile(t_Interval_length, part_id, point_data, signal);

//just keep the relevant part of the stationary interval for the local process
  for (i=0; i<(Nnz*m); i++)
    indices[i]=point_data[i+number_in_interval*t_Interval_length_loc];

  for(i=0; i<(m); i++)
    b[i] = signal[i+number_in_interval*t_Interval_length_loc];

  free(point_data);
  free(signal);

}
else { //for the case we dont need to share
//Read the relevants raw inputs data from files distributed by stationary period
//note: Work only for no sharing stationnary interval

  for (k=0; k < Nb_t_Intervals_loc; ++k) {
    point_data = indices + t_Interval_length*Nnz*k;
    signal = b + t_Interval_length*k;
    polar = pol_ang + t_Interval_length*k;
    part_id = Nb_t_Intervals_loc*rank + k;
  for (i=0; i < t_Interval_loop; ++i) {
    ioReadfile_pol(t_Interval_length, part_id, point_data+t_Interval_length_true*Nnz*i, signal+t_Interval_length_true*i, polar+t_Interval_length_true*i);
  }

//  ioReadfile(t_Interval_length, part_id, point_data, signal);

  }//end of the loop over the intervals
  // double alpha = 0;
  // srand(time(NULL));


  for (i=0; i<(Nnz*m); i+=3){
    // alpha = 2 * M_PI * rand()/((double)RAND_MAX);
    // printf("%d : pol_ang = %f, cos = %f, sin = %f\n",i/3, pol_ang[(int)i/3], cos(2*pol_ang[(int)i/3]), sin(2*pol_ang[(int)i/3]));
    values[i+1] = cos(2*pol_ang[(int)i/3]);//cos(alpha);//1 + rand()/((double)RAND_MAX);//cos(pol_ang[(int)i/3]);
    values[i+2] = sin(2*pol_ang[(int)i/3]);//sin(alpha);//1 + rand()/((double)RAND_MAX);//sin(pol_ang[(int)i/3]);
  }

}//end if

}//End if

/* ###################### TEST OF RING SCHEME #######################*/
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("\n\n########## TESTING RING COM SCHEME ###############\n\n");
  pointing_commflag=1;

//Pointing matrix init
  if(rank==0)
    printf("##### Initializing Pointing matrix ###############\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if(rank==0)
    printf("[rank %d] RING: Initializing pointing matrix time = %lf \n", rank, t1-t0);

  n = A.lcount; // number of local pixels
  v = (double *) malloc(n*sizeof(double));

//Pointing operation
  if(rank==0)
    printf("##### Product: pointing operation #################\n");
  MPI_Barrier(MPI_COMM_WORLD);
  st=MPI_Wtime();
  TrMatVecProd(&A, b, v, 0);		//  v = At*b
  MPI_Barrier(MPI_COMM_WORLD);
  t=MPI_Wtime();
  if(rank==0)
  printf("##### End pointing operation ####################\n");
  if (rank==0) {
    printf("[rank %d] RING: Total pointing operation time=%lf \n", rank, t-st);
  }

MatFree(&A); //free memory

/* ###################### TEST OF BUTTERFLY SCHEME #######################*/
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("\n\n########## TESTING BUTTERFLY COM SCHEME ###############\n\n");
  pointing_commflag=2;

//Pointing matrix init
  if(rank==0)
    printf("##### Initializing Pointing matrix ###############\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if(rank==0)
    printf("[rank %d] BUTTERFLY: Initializing pointing matrix time = %lf \n", rank, t1-t0);

  // n = A.lcount; // number of local pixels
  // v = (double *) malloc(n*sizeof(double));

//Pointing operation
  if(rank==0)
    printf("##### Product: pointing operation #################\n");
  MPI_Barrier(MPI_COMM_WORLD);
  st=MPI_Wtime();
  TrMatVecProd(&A, b, v, 0);		//  v = At*b
  MPI_Barrier(MPI_COMM_WORLD);
  t=MPI_Wtime();
  if(rank==0)
  printf("##### End pointing operation ####################\n");
  if (rank==0) {
    printf("[rank %d] BUTTERFLY: Total pointing operation time=%lf \n", rank, t-st);
  }

  MatFree(&A); //free memory

/* ###################### TEST OF NONBLOCKING SCHEME #######################*/
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
      printf("\n\n########## TESTING NONBLOCKING COM SCHEME ###############\n\n");
    pointing_commflag=3;

  //Pointing matrix init
    if(rank==0)
      printf("##### Initializing Pointing matrix ###############\n");
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
    MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    if(rank==0)
      printf("[rank %d] NONBLOCKING: Initializing pointing matrix time = %lf \n", rank, t1-t0);

    // n = A.lcount; // number of local pixels
    // v = (double *) malloc(n*sizeof(double));

  //Pointing operation
    if(rank==0)
      printf("##### Product: pointing operation #################\n");
    MPI_Barrier(MPI_COMM_WORLD);
    st=MPI_Wtime();
    TrMatVecProd(&A, b, v, 0);		//  v = At*b
    MPI_Barrier(MPI_COMM_WORLD);
    t=MPI_Wtime();
    if(rank==0)
    printf("##### End pointing operation ####################\n");
    if (rank==0) {
      printf("[rank %d] NONBLOCKING: Total pointing operation time=%lf \n", rank, t-st);
    }

    MatFree(&A); //free memory

/* ###################### TEST OF NOEMPTY SCHEME #######################*/
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("\n\n########## TESTING NOEMPTY COM SCHEME ###############\n\n");
  pointing_commflag=4;

//Pointing matrix init
  if(rank==0)
    printf("##### Initializing Pointing matrix ###############\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if(rank==0)
    printf("[rank %d] NOEMPTY: Initializing pointing matrix time = %lf \n", rank, t1-t0);

  // n = A.lcount; // number of local pixels
  // v = (double *) malloc(n*sizeof(double));

//Pointing operation
  if(rank==0)
    printf("##### Product: pointing operation #################\n");
  MPI_Barrier(MPI_COMM_WORLD);
  st=MPI_Wtime();
  TrMatVecProd(&A, b, v, 0);		//  v = At*b
  MPI_Barrier(MPI_COMM_WORLD);
  t=MPI_Wtime();
  if(rank==0)
  printf("##### End pointing operation ####################\n");
  if (rank==0) {
    printf("[rank %d] NOEMPTY: Total pointing operation time=%lf \n", rank, t-st);
  }

  MatFree(&A); //free memory

/* ###################### TEST OF ALLTOALLV SCHEME #######################*/
//   MPI_Barrier(MPI_COMM_WORLD);
//   if(rank==0)
//     printf("\n\n########## TESTING ALLTOALLV COM SCHEME ###############\n\n");
//   pointing_commflag=5;
//
// //Pointing matrix init
//   if(rank==0)
//     printf("##### Initializing Pointing matrix ###############\n");
//   MPI_Barrier(MPI_COMM_WORLD);
//   t0 = MPI_Wtime();
//   MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
//   MPI_Barrier(MPI_COMM_WORLD);
//   t1 = MPI_Wtime();
//   if(rank==0)
//     printf("[rank %d] ALLTOALLV: Initializing pointing matrix time = %lf \n", rank, t1-t0);
//
//   // n = A.lcount; // number of local pixels
//   // v = (double *) malloc(n*sizeof(double));
//
// //Pointing operation
//   if(rank==0)
//     printf("##### Product: pointing operation #################\n");
//   MPI_Barrier(MPI_COMM_WORLD);
//   st=MPI_Wtime();
//   TrMatVecProd(&A, b, v, 0);		//  v = At*b
//   MPI_Barrier(MPI_COMM_WORLD);
//   t=MPI_Wtime();
//   if(rank==0)
//   printf("##### End pointing operation ####################\n");
//   if (rank==0) {
//     printf("[rank %d] ALLTOALLV: Total pointing operation time=%lf \n", rank, t-st);
//   }
//
//   MatFree(&A); //free memory
//
// /* ###################### TEST OF ALLREDUCE SCHEME #######################*/
//   MPI_Barrier(MPI_COMM_WORLD);
//   if(rank==0)
//     printf("\n\n########## TESTING ALLREDUCE COM SCHEME ###############\n\n");
//   pointing_commflag=6;
//
// //Pointing matrix init
//   if(rank==0)
//     printf("##### Initializing Pointing matrix ###############\n");
//   MPI_Barrier(MPI_COMM_WORLD);
//   t0 = MPI_Wtime();
//   MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
//   MPI_Barrier(MPI_COMM_WORLD);
//   t1 = MPI_Wtime();
//   if(rank==0)
//     printf("[rank %d] ALLREDUCE: Initializing pointing matrix time = %lf \n", rank, t1-t0);
//
//   // n = A.lcount; // number of local pixels
//   // v = (double *) malloc(n*sizeof(double));
//
// //Pointing operation
//   if(rank==0)
//     printf("##### Product: pointing operation #################\n");
//   MPI_Barrier(MPI_COMM_WORLD);
//   st=MPI_Wtime();
//   TrMatVecProd(&A, b, v, 0);		//  v = At*b
//   MPI_Barrier(MPI_COMM_WORLD);
//   t=MPI_Wtime();
//   if(rank==0)
//   printf("##### End pointing operation ####################\n");
//   if (rank==0) {
//     printf("[rank %d] ALLREDUCE: Total pointing operation time=%lf \n", rank, t-st);
//   }
//
//   MatFree(&A); //free memory

/* ###################### TEST OF BUTTERFLY_BLOCKING_1 SCHEME #######################*/
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("\n\n########## TESTING BUTTERFLY_BLOCKING_1 COM SCHEME ###############\n\n");
  pointing_commflag=7;

//Pointing matrix init
  if(rank==0)
    printf("##### Initializing Pointing matrix ###############\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if(rank==0)
    printf("[rank %d] BUTTERFLY_BLOCKING_1: Initializing pointing matrix time = %lf \n", rank, t1-t0);

  // n = A.lcount; // number of local pixels
  // v = (double *) malloc(n*sizeof(double));

//Pointing operation
  if(rank==0)
    printf("##### Product: pointing operation #################\n");
  MPI_Barrier(MPI_COMM_WORLD);
  st=MPI_Wtime();
  TrMatVecProd(&A, b, v, 0);		//  v = At*b
  MPI_Barrier(MPI_COMM_WORLD);
  t=MPI_Wtime();
  if(rank==0)
  printf("##### End pointing operation ####################\n");
  if (rank==0) {
    printf("[rank %d] BUTTERFLY_BLOCKING_1: Total pointing operation time=%lf \n", rank, t-st);
  }

  MatFree(&A); //free memory

/* ###################### TEST OF BUTTERFLY_BLOCKING_2 SCHEME #######################*/
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("\n\n########## TESTING BUTTERFLY_BLOCKING_2 COM SCHEME ###############\n\n");
  pointing_commflag=8;

//Pointing matrix init
  if(rank==0)
    printf("##### Initializing Pointing matrix ###############\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if(rank==0)
    printf("[rank %d] BUTTERFLY_BLOCKING_2: Initializing pointing matrix time = %lf \n", rank, t1-t0);

  // n = A.lcount; // number of local pixels
  // v = (double *) malloc(n*sizeof(double));

//Pointing operation
  if(rank==0)
    printf("##### Product: pointing operation #################\n");
  MPI_Barrier(MPI_COMM_WORLD);
  st=MPI_Wtime();
  TrMatVecProd(&A, b, v, 0);		//  v = At*b
  MPI_Barrier(MPI_COMM_WORLD);
  t=MPI_Wtime();
  if(rank==0)
  printf("##### End pointing operation ####################\n");
  if (rank==0) {
    printf("[rank %d] BUTTERFLY_BLOCKING_2: Total pointing operation time=%lf \n", rank, t-st);
  }

MatFree(&A); //free memory

/* ###################### TEST OF NOEMPTYSTEPRING SCHEME #######################*/
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    printf("\n\n########## TESTING NOEMPTYSTEPRING COM SCHEME ###############\n\n");
  pointing_commflag=9;

//Pointing matrix init
  if(rank==0)
    printf("##### Initializing Pointing matrix ###############\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  MatInit( &A, m, Nnz, indices, values, pointing_commflag, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  if(rank==0)
    printf("[rank %d] NOEMPTYSTEPRING: Initializing pointing matrix time = %lf \n", rank, t1-t0);

  // n = A.lcount; // number of local pixels
  // v = (double *) malloc(n*sizeof(double));

//Pointing operation
  if(rank==0)
    printf("##### Product: pointing operation #################\n");
  MPI_Barrier(MPI_COMM_WORLD);
  st=MPI_Wtime();
  TrMatVecProd(&A, b, v, 0);		//  v = At*b
  MPI_Barrier(MPI_COMM_WORLD);
  t=MPI_Wtime();
  if(rank==0)
  printf("##### End pointing operation ####################\n");
  if (rank==0) {
    printf("[rank %d] NOEMPTYSTEPRING: Total pointing operation time=%lf \n", rank, t-st);
  }

  MatFree(&A); //free memory


  free(values);
  free(indices);
  free(v);
  free(b);
  MPI_Finalize();

  return 0;
 }


int partition(int64_t *gif, int *m, int64_t M, int rank, int size){
  int64_t r, k;
  k = M / size;
  r = M - k*size;
  if( rank < r){
    *gif = (k+1) * rank;
    *m = k+1;
  }
  else{
    *gif = r*(k+1) + k*(rank-r);
    *m = k;
    }
  return 0;
}
