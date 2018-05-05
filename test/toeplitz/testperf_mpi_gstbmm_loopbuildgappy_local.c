#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "toeplitz.h"
#include <time.h>

// Using examples for the symmetric Toeplitz block-diagonal matrix-matrix product routine (mpi_stbmm).
// Test cases examples are showning using differents parameters.
// Test number as first argument
// -1 (default) : all test will be run. 
// 1 : nblocks = 7, m = 1, n[i] = 500 and lambda = [10 51 159 73 102 55 49]; 
// 2 : nblocks = 4, m = 2, n[i] = 100 and lambda[i] = 3;

char CHAR_RW='\0';  //global variable for write mode
int IDLOOP;
int LGAP;

int PRINTOUT;
int FIRSTPRINT;
FILE* FILEOUT;
FILE* FILEOUT_MOY;
int IDLOOP_NB;

int main (int argc, char *argv[])
{
  int test=0;  //index of the choosen test
  int nb_test=1;  //number of example tests

  if (argc>1 && atoi(argv[1])>0) {
    test=atoi(argv[1]);

    if(test>nb_test) {
      printf("test %d not defined, use a test number between 1 and %d.\n", test, nb_test);
      printf("test_mpi_stmm_example [testnumber] [w]\n" );
      exit(1);}

  if(argc>2 && strcmp(argv[2],"w")==0) {
     CHAR_RW=argv[2][0];
     printf("Write mode...\n");
  }}
  else {
    printf( "test_mpi_gstbmm_example [testnumber] [w]\n" );
    exit(1);
  }


  int i,j;  //some indexes
  int *n, m, *lambda, nb_blocks; //row dimension and column dimension of the toeplitz matrix, number of blocks
  int *lv; //first indexes interval between two following blocks

  // MPI parameters
  int rank, size;
  MPI_Init(&argc, &argv);                    //Initialise MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );    //Id rank of the process
  MPI_Comm_size( MPI_COMM_WORLD, &size );    //number of processes 
  int size0=1;
  int rank0=0;
  int *displs = (int*) calloc(size0,sizeof(int));  //indexes to scatter over the processes
  int *nranks = (int*) calloc(size0,sizeof(int));  //size of each process 

  int idloop;


  FILE* fileout;
  char filenameout [1024];
  PRINTOUT=rank;
  sprintf(filenameout,"output-test-%d.txt", PRINTOUT);
  fileout = fopen(filenameout, "w");
 
  FILEOUT = fileout;

  FILE* fileout_moy;
  char filenameout_moy [1024];
  sprintf(filenameout_moy,"output-test-moy.txt");
  fileout_moy = fopen(filenameout_moy, "w");

  FILEOUT_MOY = fileout_moy;


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
    printf("idloop \t lgap \t cputime (s)\n");
    }

 //PRINTOUT=1;
 FIRSTPRINT=1;

//parameters
//    if(rank==0)
//      printf("Parameters setup...\n");

//  for (idloop=0;idloop<333;idloop++) {
  IDLOOP_NB=30;
  for (idloop=1;idloop<IDLOOP_NB;idloop++) {

    if(rank==0) {
 //   printf("----------------------------\n");
 //   printf("loop %d\n", idloop);
 //   printf("----------------------------\n");
    }

    IDLOOP = idloop;

    nb_blocks = 10; //10; //number of Toeplitz blocks on the diagonal of the matrix
    m = 1; //global column dimension
    n      = (int *) calloc(nb_blocks, sizeof(int));
    lambda = (int *) calloc(nb_blocks, sizeof(int));
    lv = (int *) calloc(nb_blocks, sizeof(int));  //define the space between the firsts indexes of
                                                  // two following blocks

    for(i=0;i<nb_blocks;i++)
      n[i] = 1000; //row dimension of each toeplitz block

    for (i=0;i<nb_blocks;i++) 
      lv[i]=n[i];//2*n[i];

    //half bandwith size for each toeplitz block (meaning no-zeros elements of the first row)
    for (i=0;i<nb_blocks;i++)
      lambda[i] = 3;


    int nsample=0; //global row dimension of the matrix
    for(j=0;j<nb_blocks;j++)
      nsample += lv[j];

//    if (rank==0)
//      printf("nsample=%d\n", nsample);

    int nrank  = (nsample*m)/size0;
//    int nrank  = (nsample*m);


    for(i=0;i<size0;i++){
      displs[i]=i*nrank; //indexes to scatter over the processes
      nranks[i]=nrank; //size of each process
    }

    //displs[size-1] = displs[size-2]+nranks[size-2]; //last process index
    //nranks[size-1] = nsample*m-displs[size-1]; //size of the last process


  //gaps description
    int ngap;
    int *id0gap;
    int *lgap;

    ngap=nb_blocks;
    id0gap = (int *) calloc(ngap, sizeof(int));
    lgap = (int *) calloc(ngap, sizeof(int));

  for (i=0;i<ngap;i++){
    id0gap[i]=i*lv[i]+n[i]/2 - (idloop)*10*lambda[0]/2;
    lgap[i]=(idloop)*10*lambda[0];
  }

  LGAP = lgap[0];

    test_toeplitz(n, m, lambda, nb_blocks, displs, nranks, lv, id0gap, lgap, ngap);

 }//end for
    }  //End test==1


  if(test==2){
    if(rank==0) {
    printf("============================\n");
    printf("test 2\n");
    printf("----------------------------\n");
    printf("idloop \t lgap \t cputime (s)\n");
    }

 //PRINTOUT=1;
 FIRSTPRINT=1;

  IDLOOP_NB=30;
  for (idloop=1;idloop<IDLOOP_NB;idloop++) {

    if(rank==0) {
 //   printf("----------------------------\n");
 //   printf("loop %d\n", idloop);
 //   printf("----------------------------\n");
    }

    IDLOOP = idloop;

    nb_blocks = 10; //10; //number of Toeplitz blocks on the diagonal of the matrix
    m = 1; //global column dimension
    n      = (int *) calloc(nb_blocks, sizeof(int));
    lambda = (int *) calloc(nb_blocks, sizeof(int));
    lv = (int *) calloc(nb_blocks, sizeof(int));  //define the space between the firsts indexes of
                                                  // two following blocks

    for(i=0;i<nb_blocks;i++)
      n[i] = 1000; //row dimension of each toeplitz block

    for (i=0;i<nb_blocks;i++)
      lv[i]=n[i];//2*n[i];

    //half bandwith size for each toeplitz block (meaning no-zeros elements of the first row)
    for (i=0;i<nb_blocks;i++)
      lambda[i] = 30;


    int nsample=0; //global row dimension of the matrix
    for(j=0;j<nb_blocks;j++)
      nsample += lv[j];

//    if (rank==0)
//      printf("nsample=%d\n", nsample);

    int nrank  = (nsample*m)/size0;
//    int nrank  = (nsample*m);


    for(i=0;i<size0;i++){
      displs[i]=i*nrank; //indexes to scatter over the processes
      nranks[i]=nrank; //size of each process
    }

    //displs[size-1] = displs[size-2]+nranks[size-2]; //last process index
    //nranks[size-1] = nsample*m-displs[size-1]; //size of the last process


  //gaps description
    int ngap;
    int *id0gap;
    int *lgap;

    ngap=nb_blocks;
    id0gap = (int *) calloc(ngap, sizeof(int));
    lgap = (int *) calloc(ngap, sizeof(int));

  for (i=0;i<ngap;i++){
    id0gap[i]=i*lv[i]+n[i]/2 - (idloop)*lambda[0]/2;
    lgap[i]=(idloop)*lambda[0];
  }

  LGAP = lgap[0];

    test_toeplitz(n, m, lambda, nb_blocks, displs, nranks, lv, id0gap, lgap, ngap);

 }//end for
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
int test_toeplitz(int *n, int m, int *lambda, int nb_blocks, int *displs, int *nranks, int *lv, int *id0gap, int *lgap, int ngap)
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

 // if(rank==0) {
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

 int flag_IwantVlocal=1;
 if (flag_scatter==1 || flag_IwantVlocal==1) {
  V  = (double *) calloc(lvtot*m ,sizeof(double));
  //definitiion of the input matrix V which will be next scattered over the processes
  //This is just for the case you want first create the all vector in the rank zero
  for(i=0;i<lvtot*m;i++)
    V[i] = i%20+1;//rand()/((double) RAND_MAX);
  }//end (flag_scatter==1)
  //}//end rank==0


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
//  MPI_Bcast(T, lambdatot, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //times variables
  double t1,t2;

//Product computation
  int nrow = lvtot;  //global row size
  int id0 = displs[rank];  //first index of the process
  int local_V_size=nranks[rank];  //data size for the process

  //allocations of the local V
//  double *Vrank = (double *) calloc(local_V_size, sizeof(double));  //local V input for the process
/*
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
*/
//  if (rank==0)
//    printf("Product computation...\n");
///  MPI_Barrier(MPI_COMM_WORLD);
  t1=MPI_Wtime();
  //perform local product using midas DA lib
//  mpi_gstbmm(&Vrank, n, m, nrow, T, nb_blocks, nb_blocks, lambda, idv, id0, local_V_size, id0gap, lgap, ngap, MPI_COMM_WORLD);
  int m_rowwise=1;
  mpi_gstbmm(&V, n, m, nrow, m_rowwise, T, nb_blocks, nb_blocks, lambda, idv, 0, lvtot*m, id0gap, lgap, ngap, MPI_COMM_WORLD);
//  mpi_gstbmm(&Vrank, n, m, nrow, T, nb_blocks, nb_blocks, lambda, idv, id0, local_V_size, id0gap, lgap, ngap, MPI_COMM_WORLD);
// sleep(1);

  t2=  MPI_Wtime();


//  if (rank==0)
//    printf("Result...\n");

///  MPI_Barrier(MPI_COMM_WORLD);

//  printf("[rank %d] Computation time : %f s.\n", rank, t2-t1);

//  printf("%d \t %d \t %f\n", IDLOOP, LGAP, t2-t1);
  if (FIRSTPRINT==1) { 
    fprintf(FILEOUT,"idloop \t lgap \t cputime (s)\n");
  }

  fprintf(FILEOUT, "%d \t %d \t %f\n", IDLOOP, LGAP, t2-t1);


//  printf("[rank %d] Finish !\n", rank);


  int flag_reduceforstats=1;
  int *displs_moy = (int*) calloc(size,sizeof(int));
  int *nranks_moy = (int*) calloc(size,sizeof(int));

  double *cputime_moy;
  cputime_moy  = (double *) calloc(IDLOOP_NB ,sizeof(double));

  double *cputime_var_nobiais;
  cputime_var_nobiais  = (double *) calloc(IDLOOP_NB ,sizeof(double));

  double *cputime_allrank;
  cputime_allrank  = (double *) calloc(size ,sizeof(double));
  double *cputime_loc;
  cputime_loc  = (double *) calloc(1 ,sizeof(double));
  cputime_loc[0] = t2-t1;



 if (flag_reduceforstats==1) {

  for(i=0;i<size;i++){
    displs_moy[i]=i*1; //indexes to scatter over the processes
    nranks_moy[i]=1; //size of each process
  }

  MPI_Gatherv(cputime_loc, nranks_moy[rank], MPI_DOUBLE, cputime_allrank, nranks_moy, displs_moy, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

if (rank==0) {
//  for(i=0;i<size;i++){
//    printf("IDLOOP=%d, cputime_allrank[%d]=%lf\n", IDLOOP, i, cputime_allrank[i]);
//  }

  cputime_moy[IDLOOP]=0;
  for(i=0;i<size;i++){
    cputime_moy[IDLOOP] += cputime_allrank[i];
  }
  cputime_moy[IDLOOP] = cputime_moy[IDLOOP]/(size);


  cputime_var_nobiais[IDLOOP]=0;
  for(i=0;i<size;i++){
    cputime_var_nobiais[IDLOOP] += (cputime_allrank[i]-cputime_moy[IDLOOP])*(cputime_allrank[i]-cputime_moy[IDLOOP]);
  }
//  cputime_var_nobiais[IDLOOP] = cputime_var_nobiais[IDLOOP]/(size-1) - (cputime_moy[IDLOOP]*cputime_moy[IDLOOP])*size/(size-1);  //size-1 for nobiais estimator
//    fprintf(FILEOUT_MOY,"%lf\n",cputime_var_nobiais[IDLOOP]);

  cputime_var_nobiais[IDLOOP] = cputime_var_nobiais[IDLOOP]/(size-1);
  cputime_var_nobiais[IDLOOP] = sqrt( cputime_var_nobiais[IDLOOP] );


  if (FIRSTPRINT==1) {
    fprintf(FILEOUT_MOY,"idloop \t lgap \t cputime_moy (s) \t cputime_var_nobiais\n");
  }

//  fprintf(FILEOUT_MOY, "IDLOOP=%d, cputime_moy[IDLOOP]=%lf\n", IDLOOP, cputime_moy[IDLOOP]);
  fprintf(FILEOUT_MOY, "%d \t %d \t %lf \t %lf\n", IDLOOP, LGAP, cputime_moy[IDLOOP], cputime_var_nobiais[IDLOOP]);


}//End Rank
}//end flag_reduceforstats==1

  FIRSTPRINT=0;


  free(T);
//  free(Vrank);


  return 0;
}



