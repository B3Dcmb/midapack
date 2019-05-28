// MAPPRAISER vdev
// Massively parallel iterative map-making code using the Midapack library v1.2b, Nov 2012
// The routine processes expanded pointing, signal and noise data arrays and produces local maps in binary files
// N.B: Future updates will produce the global map directly

/** @file   mappraiser.c
    @author Hamza El Bouhargani
    @date   May 2019 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "midapack.h"

int partition(int64_t *gif, int *m, int64_t M, int rank, int size);

void MLmap(MPI_Comm comm, int Nb_t_Intervals, int t_Interval_length, int Nnz, void *pix, void *pixweights, void *signal, int lambda, void *invcov0)
{
  int64_t	M;       //Global number of rows
  int		m;  //local number of rows of the pointing matrix A
  int64_t	gif;			//global indice for the first local line
  int		i, j, k;
  int          K;	  //maximum number of iteration for the PCG
  double 	tol;			//residual tolerence for the PCG
  Mat	A;			        //pointing matrix structure
  int 		*id0pix, *ll;
  int 		pointing_commflag;	//option for the communication scheme for the pointing matrix
  double	*x;	//pixel domain vectors
  double	st, t;		 	//timer, start time
  int 		rank, size;

  printf("\n############# MAPPRAISER : MidAPack PaRAllel Iterative Sky EstimatoR vDev, May 2019 ################\n");
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  printf("rank = %d, size = %d",rank,size);
  fflush(stdout);

//communication scheme for the pointing matrix  (to move in .h)
  pointing_commflag=1; //2==BUTTERFLY - 1==RING

//PCG parameters
  tol=pow(10,-6);
  K=500;

//total length of the time domaine signal
  M = (int64_t) Nb_t_Intervals*t_Interval_length;
  if(rank==0){
    printf("[rank %d] M=%ld\n", rank, M);
  }
  fflush(stdout);

//compute distribution indexes over the processes
  partition(&gif, &m, M, rank, size);

//Print information on data distribution
  int Nb_t_Intervals_loc = ceil( Nb_t_Intervals*1.0/size );
  int nb_proc_shared_one_interval = max(1, size/Nb_t_Intervals );
  if(rank==0)
    printf("[rank %d] nb_proc_shared_one_interval=%d\n", rank, nb_proc_shared_one_interval );
  int t_Interval_length_loc = t_Interval_length/nb_proc_shared_one_interval;
  if(rank==0){
    printf("[rank %d] size=%d \t m=%d \t Nb_t_Intervals=%d \t t_Interval_length=%d\n", rank, size, m, Nb_t_Intervals, t_Interval_length );
    printf("[rank %d] Nb_t_Intervals_loc=%d \t t_Interval_length_loc=%d\n", rank, Nb_t_Intervals_loc , t_Interval_length_loc);
  }
  fflush(stdout);

//Pointing matrix init
  st=MPI_Wtime();
  A.trash_pix =0;
  MatInit( &A, m, Nnz, pix, pixweights, pointing_commflag, comm);
  t=MPI_Wtime();
  if (rank==0) {
    printf("[rank %d] Initializing pointing matrix time=%lf \n", rank, t-st);
  }
  fflush(stdout);
  // printf("A.lcount = %d\n", A.lcount);

//Build pixel-to-time domain mapping
  st=MPI_Wtime();
  id0pix = (int *) malloc(A.lcount/(A.nnz) * sizeof(int)); //index of the last time sample pointing to each pixel
  ll = (int *) malloc(m * sizeof(int)); //linked list of time samples indexes

  //initialize the mapping arrays to -1
  for(i=0; i<m; i++){
    ll[i] = -1;
  }
  for(j=0; j<A.lcount/(A.nnz); j++){
    id0pix[j] = -1;
  }
  //build the linked list chain of time samples corresponding to each pixel
  for(i=0;i<m;i++){
    if(id0pix[A.indices[i*A.nnz]/(A.nnz)] == -1)
      id0pix[A.indices[i*A.nnz]/(A.nnz)] = i;
    else{
      ll[i] = id0pix[A.indices[i*A.nnz]/(A.nnz)];
      id0pix[A.indices[i*A.nnz]/(A.nnz)] = i;
    }
  }
  A.id0pix = id0pix;
  A.ll = ll;
  t=MPI_Wtime();
  if (rank==0) {
    printf("[rank %d] Total pixel-to-time domain mapping time=%lf \n", rank, t-st);
  }
  fflush(stdout);

//PCG beginning vector input definition for the pixel domain map (MatInit gives A.lcount)
  int *lhits;
  double *cond;
  x   = (double *) malloc(A.lcount*sizeof(double));
  cond = (double *) malloc((int)(A.lcount/3)*sizeof(double));
  lhits = (int *) malloc((int)(A.lcount/3) * sizeof(int));
  for(j=0; j<A.lcount; j++){
    x[j] = 0.;
    if(j%3 == 0){
      lhits[(int)(j/3)] = 0;
      cond[(int)(j/3)] = 0.;
    }
  }

//Create piecewise Toeplitz matrix
//specifics parameters:
  int nb_blocks_tot = Nb_t_Intervals;
  int n_block_avg = M/nb_blocks_tot;  //should be equal to t_Intervals_length in the current config
                                      //because we dont have flotting blocks
  int lambda_block_avg = lambda;

//flags for Toeplitz product strategy
  Flag flag_stgy;
  flag_stgy_init_auto(&flag_stgy);

//to print something on screen
  flag_stgy.flag_verbose=1;

//define Toeplitz blocks list and structure for Nm1
  Block *tpltzblocks;
  Tpltz Nm1;

//dependants parameters:
  int64_t nrow = M;
  int mcol = 1;

  int64_t id0 = gif;
  int local_V_size = m;

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
  defineBlocks_avg(tpltzblocks, invcov0, nb_blocks_loc, n_block_avg, lambda_block_avg, id0 );
  defineTpltz_avg( &Nm1, nrow, 1, mcol, tpltzblocks, nb_blocks_loc, nb_blocks_tot, id0, local_V_size, flag_stgy, comm);

//print Toeplitz parameters for information
  if (rank==0 | rank==1) {
    printf("[rank %d] size=%d, nrow=%ld, local_V_size=%ld, id0=%ld \n", rank, size, nrow, local_V_size, id0);
    printf("[rank %d] nb_blocks_tot=%d, nb_blocks_loc=%d, n_block_avg=%d, lambda_block_avg=%d \n", rank, nb_blocks_tot, nb_blocks_loc, n_block_avg, lambda_block_avg);
    printf("[rank %d] nb_proc_shared_a_block=%d, nb_comm=%d \n", rank, nb_proc_shared_a_block, nb_comm);
  }

  MPI_Barrier(comm);
   if(rank==0)
 printf("##### Start PCG ####################\n");
 fflush(stdout);

  st=MPI_Wtime();
// Conjugate Gradient
  PCG_GLS_true( &A, Nm1, x, signal, cond, lhits, tol, K);
  MPI_Barrier(comm);
  t=MPI_Wtime();
   if(rank==0)
 printf("##### End PCG ####################\n");
  if (rank==0) {
    printf("[rank %d] Total PCG time=%lf \n", rank, t-st);
  }
  fflush(stdout);

//write output to binaries files:
//This part will be updated later to create directly one global map in a fits file readable by healpy
  st=MPI_Wtime();
  int mapsize = A.lcount-(A.nnz)*(A.trash_pix);
  int map_id = rank;

  int *lstid;
  // int *lhits;
  lstid = (int *) calloc(mapsize, sizeof(int));
  for(i=0; i< mapsize; i++){
    lstid[i] = A.lindices[i+(A.nnz)*(A.trash_pix)];
  }
  ioWritebinfile( mapsize, map_id, lstid, x, cond, lhits);

//Write some parameters in txt file:
  //output file:
  FILE* file;
  char filenametxt [1024];
  sprintf(filenametxt,"/global/cscratch1/sd/elbouha/data_TOAST/output/mapout%01d.txt", map_id);
  file = fopen(filenametxt, "w");
  fprintf(file, "%d\n", size );
  fprintf(file, "%ld\n", gif );
  fprintf(file, "%d\n", m );
  fprintf(file, "%d\n", mapsize );
  fprintf(file, "%d\n", Nb_t_Intervals );
  fprintf(file, "%d\n", t_Interval_length );
  fprintf(file, "%d\n", LambdaBlock );
  fprintf(file, "%d\n", Nb_t_Intervals_loc );
  fprintf(file, "%d\n", rank );
  fprintf(file, "size idp m A.lcount Nb_t_Intervals t_Interval_length LambdaBlock Nb_t_Intervals_loc rank\n" );
  fclose(file);

  t=MPI_Wtime();
  if (rank==0) {
    printf("[rank %d] Write output files time=%lf \n", rank, t-st);
  }
  st=MPI_Wtime();

  MatFree(&A);                                                //free memory
  free(x);
  free(cond);
  free(lhits);
  free(tpltzblocks);
  free(id0pix);
  free(ll);
  MPI_Barrier(comm);
  t=MPI_Wtime();
  if (rank==0) {
    printf("[rank %d] Free memory time=%lf \n", rank, t-st);
  }
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
