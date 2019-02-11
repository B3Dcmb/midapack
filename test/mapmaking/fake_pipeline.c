// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// The routine reads data from binary files and writes the results in distributed binary files

/** @file   fake_pipeline.c
    @author Hamza El Bouhargani
    @date   February 2019 */


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

int main(int argc, char *argv[])
{

  int64_t	M;       //Global number of rows
  int 		N, Nnz;  //of columns, of non-zeros values per column for the pointing matrix A
  int		m, n;  //local number of rows, of columns for the pointing matrix A
  int64_t	gif;			//global indice for the first local line
  int		i, j, k;
  int           K;	                //maximum number of iteration for the PCG
  double 	tol;			//residual tolerence for the PCG
  Mat	A;			        //pointing matrix structure
  int 		*indices;
  double 	*values;
  int 		pointing_commflag ;	//option for the communication scheme for the pointing matrix
  double	*b, *Ag, *Ad, *wghts; 	 	//temporal domain vectors
  double	*x, *g, *d, *Ax_b;	//pixel domain vectors
  double *I, *Q, *U;
  double        alpha, beta, gamma, resold, resnew;
  double 	localreduce;
  double	st, t;		 	//timer, start time
  int 		output, timer, info;
  int 		rank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("rank = %d, size = %d",rank,size);

//communication scheme for the pointing matrix  (to move in .h)
  pointing_commflag=2; //2==BUTTERFLY - 1==RING

//global data caracteristics
  int Nb_t_Intervals = 1;//1352;//128;//2;//256;//8;           //total number of stationnary intervals
  int t_Interval_length = 8;//1749900;//17899900;//1431992;//139992; //1431992;//2863984;//1431992;//pow(2,25);//pow(2,25);          //length for each stationnary interval
  int t_Interval_length_true = 8;//1749900;//17899900;//1431992;//139992;//1431992;//2863984;//1431992;//pow(2,20);
  int LambdaBlock = pow(2,0);//pow(2,14)+1;  //lambda length for each stationnary interval
  Nnz=3;

// PCG parameters
  tol=pow(10,-6);
  K=500;

//total length of the time domaine signal
  M = (int64_t) Nb_t_Intervals*t_Interval_length ;

  printf("[rank %d] M=%ld\n", rank, M);

//compute distribution indexes over the processes
  partition(&gif, &m, M, rank, size);

//input data memory allocation
  indices  = (int *) malloc(Nnz*m * sizeof(int));     //for pointing matrix indices
  b   = (double *) malloc(m*sizeof(double));    //full raw data vector for the signal
  wghts = (double *) malloc(Nnz*m * sizeof(double));
  I = (double *) malloc(3*sizeof(double));
  Q = (double *) malloc(3*sizeof(double));
  U = (double *) malloc(3*sizeof(double));

// We generate fake data simulating the observation of three pixels 0, 1 and 2, where one pixel is badly conditioned
// Each pixel is observed exactly three times but the bad one is observed two times only
// Fake map
  I[0] = 200; I[2] = -300; I[1] = 178;
  Q[0] = 15; Q[2] = 10; Q[1] = -6;
  U[0] = -13; U[2] = 0; U[1] = 3;
// Fake pointing Data
for(i=0;i<m;i++){
  for(j=0;j<Nnz;j++){
    indices[i*Nnz+j] = Nnz*(i%3)+j;
  }
  wghts[i*Nnz] = 1;
  wghts[i*Nnz+1] = cos(2*(M_PI/3)*(int)i/3);
  wghts[i*Nnz+2] = sin(2*(M_PI/3)*(int)i/3);
  b[i] = I[(int)indices[i*Nnz]/Nnz] + Q[(int)indices[indices[i*Nnz]]/Nnz] * wghts[i*Nnz+1] + U[(int)indices[indices[i*Nnz]]/Nnz] * wghts[i*Nnz+2];
}

// Plot input data to check if everything is correct
for(i=0;i<m*Nnz;i++){
  printf("indices[%d] = %d, ",i,indices[i]);
}
printf("\n");
for(i=0;i<m;i++){
  for(j=0;j<Nnz;j++){
    printf("weights[%d] = %.2f, ",i*Nnz+j,wghts[i*Nnz+j]);
  }
  printf("\n");
}
for(i=0;i<m;i++){
  printf("b[%d] = %.2f\n",i,b[i]);
}


//Pointing matrix init
  MatInit( &A, m, Nnz, indices, wghts, pointing_commflag, MPI_COMM_WORLD);

// PCG begining vector input definition for the pixel domain map (MatInit gives A.lcount)
  x   = (double *) malloc(A.lcount*sizeof(double));

  for(j=0; j<A.lcount; j++){
    x[j] = 1000.;
  }



//Create piecewise Toeplitz matrix

//specifics parameters:
  int nb_blocks_tot = Nb_t_Intervals;
  int n_block_avg = M/nb_blocks_tot;  //should be equal to t_Intervals_length in this example
                                      //because we dont have flotting blocks
  int lambda_block_avg = LambdaBlock;

//flags for Toeplitz product strategy
  Flag flag_stgy;
  flag_stgy_init_auto(&flag_stgy);

//to print something on screen
  flag_stgy.flag_verbose=1;

//to define fixed bs:
//  flag_stgy.flag_bs = 1;
//  flag_stgy.fixed_bs = pow(2,17);

//define Toeplitz blocks list and structure for Nm1
  Block *tpltzblocks;
  Tpltz Nm1;

//dependants parameters:
  int64_t nrow = M;
  int mcol = 1;

  int64_t id0 = gif;
  int local_V_size = m;

  int Tsize = lambda_block_avg;
  double *T;  //toeplitz data storage
  T  = (double *) calloc(Tsize ,sizeof(double));

//For one identical block
  ioReadTpltzrandom( Tsize, T);


  int nb_blocks_loc;
  nb_blocks_loc = ceil( local_V_size*1.0/n_block_avg );

  double nb_blocks_loc_part =  (local_V_size*1.0)/(n_block_avg) ;



// check special cases to have exact number of local blocks
  if ((id0/n_block_avg + nb_blocks_loc) * n_block_avg < (id0+local_V_size))
    nb_blocks_loc=nb_blocks_loc+1;

  if (rank==0 | rank==1) {
    printf("M=%ld, m=%d \n", M, m);
    printf("gif = %ld \n", gif);
  }

  int nb_proc_shared_a_block = ceil( size*1.0/nb_blocks_tot );
  int nb_comm = (nb_proc_shared_a_block)-1 ;


//Block definition
  tpltzblocks = (Block *) malloc(nb_blocks_loc * sizeof(Block));
  defineBlocks_avg(tpltzblocks, T, nb_blocks_loc, n_block_avg, lambda_block_avg, id0 );
  defineTpltz_avg( &Nm1, nrow, 1, mcol, tpltzblocks, nb_blocks_loc, nb_blocks_tot, id0, local_V_size, flag_stgy, MPI_COMM_WORLD);

//print Toeplitz parameters for information
  if (rank==0 | rank==1) {
    printf("[rank %d] size=%d, nrow=%ld, local_V_size=%ld, id0=%ld \n", rank, size, nrow, local_V_size, id0);
    printf("[rank %d] nb_blocks_tot=%d, nb_blocks_loc=%d, n_block_avg=%d, lambda_block_avg=%d \n", rank, nb_blocks_tot, nb_blocks_loc, n_block_avg, lambda_block_avg);
    printf("[rank %d] nb_proc_shared_a_block=%d, nb_comm=%d \n", rank, nb_proc_shared_a_block, nb_comm);
  }


  MPI_Barrier(MPI_COMM_WORLD);
   if(rank==0)
 printf("##### Start PCG ####################\n");
  st=MPI_Wtime();
// Conjugate Gradient
  PCG_GLS_true( &A, Nm1, x, b, tol, K);


  MPI_Barrier(MPI_COMM_WORLD);
  t=MPI_Wtime();
   if(rank==0)
 printf("##### End PCG ####################\n");
  if (rank==0) {
    printf("[rank %d] Total PCG time=%lf \n", rank, t-st);
  }


//write output to binaries files:
  int mapsize= A.lcount;
  int map_id=rank;

  int *lstid;
  lstid = (int *) calloc(A.lcount, sizeof(int));
  for(i=0; i< A.lcount; i++){
    lstid[i] = A.lindices[i];
    printf("x[%d] = %.2f",i,x[i]);
  }

  MatFree(&A);                                                //free memory

  // free(indices);
  // free(values);

  free(b);
  free(x);
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
