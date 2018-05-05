#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "toeplitz.h"
#include "midapack.h"
#include <time.h>

//

char CHAR_RW='\0';  //global variable for write mode
extern int NFFT;


int main (int argc, char *argv[])
{

  // MPI parameters
  int rank, size;
  MPI_Init(&argc, &argv);                    //Initialise MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );    //Id rank of the process
  MPI_Comm_size( MPI_COMM_WORLD, &size );    //number of processes 

  int i,j,k; //some indexes
  srand (time (NULL));  //init seed

  //times variables
  double t1,t2;

  int n_thread;


#pragma omp parallel
{  n_thread = omp_get_num_threads(); }

/*
  if ((NB_OMPTHREADS <= n_thread) && (NB_OMPTHREADS != 0))
    omp_set_num_threads(NB_OMPTHREADS);

#pragma omp parallel
{  n_thread = omp_get_num_threads(); }

  printf("Using %d threads\n", n_thread);
*/


//parameters
  int   nrow = pow(2,28);//30); //28);//pow(2,25);//20);//you need to choose M so M/(nb_procs) = int to havec the same size everywhere

  int local_V_size = nrow/size;
  int id0 = rank*local_V_size;
  

  int nb_blocks_tot = pow(2,4);//16 //10;
  int n_block_avg = nrow/nb_blocks_tot;  //ca doit etre un multiple dans cet exemple, sinon la fin est sans block

  int lambda_block_avg = n_block_avg/pow(2,10);//pow(2,3);//3;
  Flag flag_stgy;
  flag_stgy_init_auto(&flag_stgy);

  int mcol = 1;  

  int Tsize = lambda_block_avg;
  double *T;  //toeplitz data storage
  T  = (double *) calloc(Tsize ,sizeof(double));

  createT(T, Tsize);

 int nb_blocks_loc;
// nb_blocks_loc =  local_V_size*1.0/n_block_avg  + (id0%n_block_avg != 0) + ((id0+local_V_size)%n_block_avg != 0);  //maximum number of blocks possible

  nb_blocks_loc = ceil( local_V_size*1.0/n_block_avg );
  double nb_blocks_loc_part =  (local_V_size*1.0)/(n_block_avg) ;
  int nb_proc_shared_a_block2 = 1./nb_blocks_loc_part;

// check special cases to have exact number of local blocks
  if ((id0/n_block_avg + nb_blocks_loc) * n_block_avg < (id0+local_V_size))
    nb_blocks_loc=nb_blocks_loc+1;

  int nb_proc_shared_a_block = ceil( size*1.0/nb_blocks_tot );
  int nb_comm = (nb_proc_shared_a_block)-1 ;

//Block definition 
  Block *tpltzblocks;
  tpltzblocks = (Block *) malloc(nb_blocks_loc * sizeof(Block));

  defineBlocks_avg(tpltzblocks, T, nb_blocks_loc, n_block_avg, lambda_block_avg, id0 );

  double *Vrank;
  Vrank = (double *) calloc(local_V_size ,sizeof(double));

  for(i=0;i<local_V_size;i++)
    Vrank[i] = rand()/((double) RAND_MAX);

//print parameters:

  if (rank==0) {
    printf("[rank %d] n_thread=%d \n", rank, n_thread);
    printf("[rank %d] size=%d, nrow=%d, local_V_size=%d, id0=%d \n", rank, size, nrow, local_V_size, id0);
    printf("[rank %d] nb_blocks_tot=%d, nb_blocks_loc=%d, n_block_avg=%d, lambda_block_avg=%d \n", rank, nb_blocks_tot, nb_blocks_loc, n_block_avg, lambda_block_avg);
    printf("[rank %d] nb_proc_shared_a_block=%d, nb_comm=%d \n", rank, nb_proc_shared_a_block, nb_comm);
    printf("[rank %d] nb_blocks_loc_part=%f, nb_proc_shared_a_block2=%d \n", rank, nb_blocks_loc_part, nb_proc_shared_a_block2);
    printf("[rank %d] \n", rank);

    if ((nrow%size != 0) || (nrow%nb_blocks_tot !=0))
      printf("[rank %d] residu check not good: nrowO/osize=%d,  nrow0/onb_blocks_tot=%d \n", rank, nrow%size,  nrow%nb_blocks_tot);
  }


//start computing product
  MPI_Barrier(MPI_COMM_WORLD);
  t1=MPI_Wtime();

//Toeplitz product
  mpi_stbmm(&Vrank, nrow, mcol, mcol, tpltzblocks, nb_blocks_loc, nb_blocks_tot, id0, local_V_size, flag_stgy, MPI_COMM_WORLD);

//  t2=  MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  t2=  MPI_Wtime();

  if (rank==0) {
    printf("[rank %d] Computation time : %f s.\n", rank, t2-t1);
    printf("[rank %d] Finish !\n", rank);
  }


  free(T);
  free(Vrank);
  MPI_Finalize();

  return 0;
}




//add fcts

int defineBlocks_avg(Block *tpltzblocks, double *T, int nb_blocks_loc, int n_block_avg, int lambda_block_avg, int id0 )
{

int i;


  for ( i=0; i<nb_blocks_loc; i++)
    tpltzblocks[i].n = n_block_avg;

  for ( i=0; i<nb_blocks_loc; i++)
    tpltzblocks[i].lambda = lambda_block_avg;

  tpltzblocks[0].idv = (id0/n_block_avg) * n_block_avg ;
  for(i=1;i<nb_blocks_loc;i++)
    tpltzblocks[i].idv = tpltzblocks[i-1].idv + tpltzblocks[i-1].n;

  for( i=0; i<nb_blocks_loc; i++) {
    tpltzblocks[i].T_block = (T);
  }

//  for(i=0; i<nb_blocks_loc; i++)
//    printf("tpltzblocks[%d].idv=%d\n", i, tpltzblocks[i].idv);

  return 0;
}




int createT(double *T, int Tsize)
{

  int i;
  srand (time (NULL));  //init seed


  if (1==1) {
  //input matrix definition of T
    for(i=0;i<Tsize;i++) 
      T[i]= rand()/((double) RAND_MAX);
    

  }
  else if (1==2) {

  //input matrix definition of T
    for(i=0;i<Tsize;i++) {
      if (i%3 == 0) {
        T[i]=10.;}
      else if (i%3 == 1) {
        T[i]=2.;}
      else if (i%3 == 2) {
        T[i]=3.;}
      else {
        T[i]=0.;//rand()/((double) RAND_MAX);
     }}

  }
  else {

#include "invtt_params.h"

  double *invtt;

  T = invtt;
//  createinvtt(invtt);

  printf("toto=%d\n", toto);

exit(1);

  } //end if`


  return 0;
}

