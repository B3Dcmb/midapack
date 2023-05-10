//
// File: testexple_mpi_gstbmm.c, version 1.1b, July 2012  
//
// Provide using examples for the symmetric Toeplitz block-diagonal matrix-matrix product
// routine including defined gaps locations (mpi_gstbmm).
// 
// Author:  Frederic Dauvergne (APC) 
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "toeplitz.h"
#include "midapack.h"
#include <time.h>

// Test cases examples are showning using differents parameters.
// Test number as first argument :
// 1 : nblocks = 7, m = 1, n[i] = 500 and lambda = [10 51 159 73 102 55 49]; 
// 2 : nblocks = 4, m = 2, n[i] = 100 and lambda[i] = 3;

char CHAR_RW='\0';  //global variable for write mode

int main (int argc, char *argv[])
{
  int test=0;  //index of the choosen test
  int nb_test=2;  //number of example tests

  if (argc>1 && atoi(argv[1])>0) {
    test=atoi(argv[1]);

    if(test>nb_test) {
      printf("test %d not defined, use a test number between 1 and %d.\n", test, nb_test);
      printf("test_mpi_gstbmm_example [testnumber] [w]\n" );
      exit(1);}

  if(argc>2 && strcmp(argv[2],"w")==0) {
     CHAR_RW=argv[2][0];
     printf("Write mode...\n");
  }}
  else {
    printf( "test_mpi_gstbmm_example [testnumber] [w]\n" );
    exit(1);
  }


  int i,j,k;  //some indexes
  int *n, m, *lambda, nb_blocks; //row dimension and column dimension of the toeplitz matrix, number of blocks
  int *lv; //first indexes interval between two following blocks
  int m_rowwise;

  // MPI parameters
  int rank, size;
  MPI_Init(&argc, &argv);                    //Initialise MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );    //Id rank of the process
  MPI_Comm_size( MPI_COMM_WORLD, &size );    //number of processes 
  int *disp_ranks_row = (int*) calloc(size,sizeof(int));  //indexes to scatter over the processes
  int *nranks = (int*) calloc(size,sizeof(int));  //size of each process 

  int ngap;
  int64_t *id0gap;
  int *lgap;

//test case 1: (nblocks = 2, m = 1, n[i] = 10, lambda[i] = 3 and lv[i]=20 ;
//data dispatched fairly between the processes)
//This is showing a basic example of how to scatter the data fairly between the processes.
//The last process contains the bigger data size due to the remainder of the euclidean division.
//here we have defined two no contigious flotting blocks for the Toeplitz block diagonal matrix.
  if(test==1){
    if(rank==0) {
    printf("============================\n");
    printf("test 1\n");
    printf("----------------------------\n");
    }

//parameters
    if(rank==0)
      printf("Parameters setup...\n");

    nb_blocks = 2; //number of Toeplitz blocks on the diagonal of the matrix
    m = 2; //global column dimension
    m_rowwise=1;
    n      = (int *) calloc(nb_blocks, sizeof(int));
    lambda = (int *) calloc(nb_blocks, sizeof(int));
    lv = (int *) calloc(nb_blocks, sizeof(int));  //define the space between the firsts indexes of
                                                  // two following blocks

    for(i=0;i<nb_blocks;i++)
      n[i] = 10; //row dimension of each toeplitz block

    for (i=0;i<nb_blocks;i++)
      lv[i]=2*n[i];

    //half bandwith size for each toeplitz block (meaning no-zeros elements of the first row)
    for (i=0;i<nb_blocks;i++)
      lambda[i] = 3;


    int lvtot=0; //global row dimension of the matrix
    for(j=0;j<nb_blocks;j++)
      lvtot += lv[j];

    if (rank==0)
      printf("lvtot=%d\n", lvtot);

    int nrank  = (lvtot*m)/size;

    for(i=0;i<size;i++){
      disp_ranks_row[i]=i*nrank; //indexes to scatter over the processes
      nranks[i]=nrank; //size of each process
    }

    disp_ranks_row[size-1] = disp_ranks_row[size-2]+nranks[size-2]; //last process index
    nranks[size-1] = lvtot*m-disp_ranks_row[size-1]; //size of the last process


  //gaps description
    ngap=2;
    id0gap = (int64_t *) calloc(ngap, sizeof(int64_t));
    lgap = (int *) calloc(ngap, sizeof(int));

    id0gap[0]=1;
    lgap[0]=4;
    id0gap[1]=6;
    lgap[1]=3;


    test_toeplitz(n, m, m_rowwise, lambda, nb_blocks, disp_ranks_row, nranks, lv, id0gap, lgap, ngap);
    }  //End test==1



//test case 1: (nblocks = 2, m = 1, n[i] = 10, lambda[i] = 3 and lv[i]=20 ;
//data dispatched fairly between the processes)
//This is showing a basic example of how to scatter the data fairly between the processes.
//The last process contains the bigger data size due to the remainder of the euclidean division.
//here we have defined two no contigious flotting blocks for the Toeplitz block diagonal matrix.
  if(test==2){
    if(rank==0) {
    printf("============================\n");
    printf("test 2\n");
    printf("----------------------------\n");
    }


//parameters
    if(rank==0)
      printf("Parameters setup...\n");

    nb_blocks = 2; //number of Toeplitz blocks on the diagonal of the matrix
    m = 1; //global column dimension
    m_rowwise=2;
    n      = (int *) calloc(nb_blocks, sizeof(int));
    lambda = (int *) calloc(nb_blocks, sizeof(int));
    lv = (int *) calloc(nb_blocks, sizeof(int));  //define the space between the firsts indexes of
                                                  // two following blocks

    for(i=0;i<nb_blocks;i++)
      n[i] = 20; //row dimension of each toeplitz block

    for (i=0;i<nb_blocks;i++)
      lv[i]=2*n[i];

    //half bandwith size for each toeplitz block (meaning no-zeros elements of the first row)
    for (i=0;i<nb_blocks;i++)
      lambda[i] = 3;


    int lvtot=0; //global row dimension of the matrix
    for(j=0;j<nb_blocks;j++)
      lvtot += lv[j];

    if (rank==0)
      printf("lvtot=%d\n", lvtot);

    int nrank  = (lvtot*m)/size;


    for(i=0;i<size;i++){
      disp_ranks_row[i]=i*nrank; //indexes to scatter over the processes
      nranks[i]=nrank; //size of each process
    }

    disp_ranks_row[size-1] = disp_ranks_row[size-2]+nranks[size-2]; //last process index
    nranks[size-1] = lvtot*m-disp_ranks_row[size-1]; //size of the last process


  //gaps description
    ngap=nb_blocks;
    id0gap = (int64_t *) calloc(ngap, sizeof(int64_t));
    lgap = (int *) calloc(ngap, sizeof(int));


  for (k=0;k<ngap;k++){
    id0gap[k]=k*lv[k]+n[k]/4;
    lgap[k]=n[k]/2;//5;2;
  }

    test_toeplitz(n, m, m_rowwise, lambda, nb_blocks, disp_ranks_row, nranks, lv, id0gap, lgap, ngap);

    }  //End test==2


  free(disp_ranks_row);
  free(nranks);
  MPI_Finalize();

  return 0;
}



// Perform the multiplication of a Toeplitz matrix by a matrix with MPI using the
//defined parameters. The input matrix V is construct on the process 0 and then 
//scattered to all the processes for the computation. We compare the obtain result
//to the exact solution and print the relative euclidiean norm residual. For this
//example, V and T are filled with constant values. The right result of TV is in the files.
int test_toeplitz(int *n, int m, int m_rowwise, int *lambda, int nb_blocks, int *disp_ranks_row, int *nranks, int *lv, int64_t *id0gap, int *lgap, int ngap)
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

  int *idv = (int *) calloc(nb_blocks, sizeof(int));  //indexes for the blocks

  int *disp_ranks;
  disp_ranks = (int *) calloc(size, sizeof(int));

  for(i=0;i<size;i++)
    disp_ranks[i] = disp_ranks_row[i]*m_rowwise ;


  int *ranks_size;
  ranks_size = (int *) calloc(size, sizeof(int));

  for(i=0;i<size;i++)
    ranks_size[i] = nranks[i]*m_rowwise ;


  int flag_scatter=0;

  idv[0]=0;
  for(i=1;i<nb_blocks;i++)
    idv[i] = idv[i-1] + lv[i-1];



//allocation of T
  T  = (double *) calloc(lambdatot ,sizeof(double));


//Blocks struture definition and allocation
  Block *tpltzblocks;
  tpltzblocks = (Block *) malloc(nb_blocks * sizeof(Block));


  for ( i=0; i<nb_blocks; i++)
    tpltzblocks[i].n = n[i];

  for ( i=0; i<nb_blocks; i++)
    tpltzblocks[i].idv = idv[i];

  for ( i=0; i<nb_blocks; i++)
    tpltzblocks[i].lambda = lambda[i];


  int lambdaShft=0;
  for( lambdaShft=i=0; i<nb_blocks; i++) {
    tpltzblocks[i].T_block = (T+lambdaShft);
    lambdaShft += lambda[i];
  }


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

 if (flag_scatter==1) {
  V  = (double *) calloc(lvtot*m*m_rowwise ,sizeof(double));
  //definitiion of the input matrix V which will be next scattered over the processes
  //This is just for the case you want first create the all vector in the rank zero

  for(k=0;k<size;k++) {
  for(j=0;j<m_rowwise;j++)
  for(i=0;i<nranks[k];i++)
    V[i+j*nranks[k]+disp_ranks[k]] = (i+disp_ranks_row[k])%20+1;//rand()/((double) RAND_MAX);
  }

  }//end (flag_scatter==1)
  }//end rank==0


  if (rank==0) {
    printf("Global parameters: lvtot=%d ; ntot=%d ; m=%d ; m_rowwise=%d \n", lvtot, ntot, m, m_rowwise);
    printf("Data distribution over the processes :\n");
    for(i=0;i<size;i++)
      printf("Rank %d : index %d ; total rank size %d\n", i, disp_ranks_row[i], nranks[i]*m_rowwise);
    printf("Scatter over the blocks :\n");
    for(i=0;i<nb_blocks;i++)
      printf("idv[%d]=%d ; n[%d]=%d ; lv[%d]=%d ; lambda[%d]=%d\n", i, idv[i], i, n[i], i, lv[i], i, lambda[i]);
  }//end rank==0

//spread T over the processes
  MPI_Bcast(T, lambdatot, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //times variables
  double t1,t2;

//Product computation
  int nrow = lvtot;  //global row size
  int id0 = disp_ranks_row[rank];  //first index of the process
  int local_V_size=nranks[rank];  //data size for the process

  //allocations of the local V
  double *Vrank = (double *) calloc(local_V_size*m_rowwise, sizeof(double));  //local V input for the process


  if (flag_scatter==1) {
  //scatter V over the processes 
  MPI_Scatterv(V, ranks_size, disp_ranks, MPI_DOUBLE, Vrank, local_V_size*m_rowwise, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  if (rank==0)
    free(V);
  }
  else {
  //direct definition for the local input matrix Vrank 
  for(j=0;j<m_rowwise;j++)
  for(i=0;i<local_V_size;i++)
    Vrank[i+j*local_V_size] = (i+id0)%20+1;//rand()/((double) RAND_MAX);

  }

  if (rank==0 && VERBOSE)
    printf("Init flag...\n");

  Flag flag_stgy;

//  flag_stgy_init_defined(&flag_stgy);
  flag_stgy_init_auto(&flag_stgy);

  if (rank==0 && VERBOSE)
    print_flag_stgy_init(flag_stgy);

  if (rank==0)
    printf("Product computation...\n");

  MPI_Barrier(MPI_COMM_WORLD);
  t1=MPI_Wtime();

  //perform local product using midas DA lib
  mpi_gstbmm(&Vrank, nrow, m, m_rowwise, tpltzblocks, nb_blocks, nb_blocks, id0, local_V_size, id0gap, lgap, ngap, flag_stgy, MPI_COMM_WORLD);

  t2=  MPI_Wtime();


  if (rank==0)
    printf("Result...\n");

  MPI_Barrier(MPI_COMM_WORLD);

  printf("[rank %d] Computation time : %f s.\n", rank, t2-t1);
  printf("[rank %d] Finish !\n", rank);


//write to file if "write mode on" 
  if (CHAR_RW=='w') {

  TV  = (double *) calloc(lvtot*m*m_rowwise, sizeof(double));

  //receive data from each process
  MPI_Gatherv(Vrank, ranks_size[rank], MPI_DOUBLE, TV, ranks_size, disp_ranks, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //output file:
  FILE* file;
  char filename [1024];
  sprintf(filename,"output-test.txt");
  file = fopen(filename, "w");


  if(rank==0){
    fprintf(file,"Test result:\n");
    fprintf(file,"ntot=%d, nrow=%d, m=%d, m_rowwise=%d\n", ntot, nrow, m, m_rowwise);
    fprintf(file,"TV=\n");



  for(k=0;k<size;k++)
  for(i=0;i<nranks[k];i++) {
  for(j=0;j<m_rowwise;j++) {
    fprintf(file,"%lf", TV[i+j*nranks[k]+disp_ranks[k]]);
    if (j!=m_rowwise) fprintf(file,"\t");
  }
  fprintf(file,"\n");
  }


  free(TV);
  fclose(file);
  }

    printf("Write done...\n");
  }//end  CHAR_RW=='w'


  free(T);
  free(Vrank);


  return 0;
}




