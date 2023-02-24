#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "toeplitz.h"
#include <time.h>

char CHAR_RW='\0';  //global variable for write mode

int main (int argc, char *argv[])
{
  int test=0;  //index of the choosen test
  int nb_test=1;  //number of example tests

  if (argc>1 && atoi(argv[1])>0) {
    test=atoi(argv[1]);

    if(test>nb_test) {
      printf("test %d not defined, use a test number between 1 and %d.\n", test, nb_test);
      printf("execfile [testnumber] [w]\n" );
      exit(1);}

  if(argc>2 && strcmp(argv[2],"w")==0) {
     CHAR_RW=argv[2][0];
     printf("Write mode...\n");
  }}
  else {
    printf( "execfile [testnumber] [w]\n" );
    exit(1);
  }


  int i,j;  //some indexes
  int *n, m, *lambda, nb_blocks; //row dimension and column dimension of the toeplitz matrix, number of blocks
  int *lv; //first indexes interval between two following blocks
  int m_rowwise;

  // MPI parameters
  int rank, size;
  MPI_Init(&argc, &argv);                    //Initialise MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );    //Id rank of the process
  MPI_Comm_size( MPI_COMM_WORLD, &size );    //number of processes 
  int *displs = (int*) calloc(size,sizeof(int));  //indexes to scatter over the processes
  int *nranks = (int*) calloc(size,sizeof(int));  //size of each process 


//test case 1: (nblocks = 4, m = 2, n[i] = 100 and lambda[i] = 3) ;
//data dispatched fairly between the processes)
//This is showing a basic example of how to scatter the data from two vector (as a matrix n
//by 2) fairly between the processes. The last process contains the bigger data size due to
//the remainder of the euclidean division. here we use m=2 for two "global" columns in a
//column-wise order.
  if(test==1){
    if(rank==0) {
    printf("============================\n");
    printf("test 1 hugeNew\n");
    printf("----------------------------\n");
    }

//parameters
//    if(rank==0)
//      printf("Parameters setup...\n");

    nb_blocks = 100; //number of Toeplitz blocks on the diagonal of the matrix
    m = 1; //global column dimension
    m_rowwise = 1; //only m_rowwise==1 work //global column dimension
    n      = (int *) calloc(nb_blocks, sizeof(int));
    lambda = (int *) calloc(nb_blocks, sizeof(int));
    lv = (int *) calloc(nb_blocks, sizeof(int));

    for(i=0;i<nb_blocks;i++)
      n[i] = 1000; //row dimension of each toeplitz block

    for (i=0;i<nb_blocks;i++) 
      lv[i]=n[i];

    //half bandwith size for each toeplitz block (meaning no-zeros elements of the first row)
    for (i=0;i<nb_blocks;i++)
      lambda[i] = 3;


    int nsample=0; //global row dimension of the matrix
    for(j=0;j<nb_blocks;j++)
      nsample += lv[j];

 //   if (rank==0)
 //     printf("nsample=%d\n", nsample);

    int nrank  = (nsample*m)/size;

    for(i=0;i<size;i++){
      displs[i]=i*nrank*m_rowwise; //indexes to scatter over the processes
      nranks[i]=nrank; //size of each process
    }

    displs[size-1] = displs[size-2]+nranks[size-2]; //last process index
    nranks[size-1] = nsample*m-displs[size-1]; //size of the last process

    test_toeplitz(n, m, m_rowwise, lambda, nb_blocks, displs, nranks, lv);
    }  //End test==1



  free(displs);
  free(nranks);
  MPI_Finalize();

  return 0;
}

// Perform the multiplication of a Toeplitz matrix by a matrix with MPI using the
//defined parameters. The input matrix V is construct on the process 0 and then 
//scattered to all the processes for the computation. We compare the obtain result
//to the exact solution and print the relative euclidiean norm residual. For this
//example, V and T are filled with constant values. The right result of TV is in the files.
int test_toeplitz(int *n, int m, int m_rowwise, int *lambda, int nb_blocks, int *displs, int *nranks, int *lv)
{

  double *T;  //input toeplitz
  double *V;  //input matrix
  double *TV;  //output matrix
 
  //MPI parameters
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 

  int i,j,k,l,p;  // some indexes
  srand (time (NULL));  //init seed

  int ntot = 0;  //global n size
  int lambdatot = 0;  //T length
  int lvtot = 0;
  int mbloc=1;  //row size of each blocks (constant for this example)

  for(i=0;i<nb_blocks;i++) {
    lambdatot += lambda[i];  //total number of lambdas
    ntot += n[i];  //total row size of the defined blocks
    lvtot += lv[i];}  //total row size of the Toeplitz block diagonal matrix

  int *idv = (int *) calloc(nb_blocks ,sizeof(int));  //indexes for the blocks

  int flag_scatter=0;

  idv[0]=0;
  for(i=1;i<nb_blocks;i++) 
    idv[i] = idv[i-1] + lv[i-1];

//allocation of T
  T  = (double *) calloc(lambdatot ,sizeof(double));

  if(rank==0) {
//allocation of V and TV
    
  //input matrix definition of T
  for(i=0;i<lambdatot;i++) {  // half band Toeplitz blocks 
    if (i%3 == 0) {
      T[i]=10;}
    else if (i%3 == 1) {
      T[i]=2;}
    else if (i%3 == 2) {
      T[i]=3;}
    else {
      T[i]=i%10;//rand()/((double) RAND_MAX);
   }}

/*
 if (flag_scatter==1) {
  V  = (double *) calloc(lvtot*m ,sizeof(double));
  //definitiion of the input matrix V which will be next scattered over the processes
  //This is just for the case you want first create the all vector in the rank zero
  for(i=0;i<lvtot*m;i++) 
    V[i] = i%20+1;//rand()/((double) RAND_MAX);
  }//end (flag_scatter==1)
*/
  }//end rank==0

/*
  if (rank==0) {
    printf("Scatter over the processes :\n");
    for(i=0;i<size;i++)
      printf("Rank %d : index %d ; size %d\n", i, displs[i], nranks[i]);
    printf("Scatter over the blocks :\n");
    for(i=0;i<nb_blocks;i++)
      printf("idv[%d]=%d ; n[%d]=%d ; lv[%d]=%d ; lambda[%d]=%d\n", i, idv[i], i, n[i], i, lv[i], i, lambda[i]);
  }//end rank==0
*/

//spread T over the processes
  MPI_Bcast(T, lambdatot, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //times variables
  double t1,t2;

//Product computation
  int nrow = lvtot;  //global row size
  int id0 = displs[rank];  //first index of the process
  int local_V_size=nranks[rank];  //data size for the process

  //allocations of the local V
  double *Vrank = (double *) calloc(local_V_size, sizeof(double));  //local V input for the process

  if (flag_scatter==1) {
  //scatter V over the processes 
  MPI_Scatterv(V, nranks, displs, MPI_DOUBLE, Vrank, local_V_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank==0)
    free(V);
  }
  else {
  //direct definition for the local input matrix Vrank 
  for(i=0;i<local_V_size;i++)
    Vrank[i] = (i+id0)%20+1;//rand()/((double) RAND_MAX);
  }

  if (rank==0)
    printf("Product computation...\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t1=MPI_Wtime();
  //perform local product using midas DA lib
  mpi_stbmm(&Vrank, n, m, nrow, m_rowwise, T, nb_blocks, nb_blocks, lambda, idv, id0, local_V_size, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  t2=  MPI_Wtime();



  if (rank==0)
    printf("Result...\n");

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank==0) {
    printf("[rank %d] Barrier - Computation time : %f s.\n", rank, t2-t1);
    printf("[rank %d] Finish !\n", rank);
  }
 
  free(T);
  free(Vrank);
 
 
  return 0;
}


