#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "toeplitz.h"
//#include "midapack.h"
#include <time.h>

// Using examples for the symmetric Toeplitz matrix-matrix product routine (mpi_stmm).
// Test cases examples are showning using differents parameters.
//Test number as first argument :
// -1 (default) : all tests will be run
// 1 : n = 5000, m = 1 and lambda = 50; data dispatched fairly between the processes
// 2 : n = 5000, m = 2 and lambda = 50; data dispatched fairly between the processes
// 3 : n = 5000, m = 5 and lambda = 50; data dispatched with variables sizes

char CHAR_RW='\0';  //global variable for write mode

int main (int argc, char *argv[])
{
  int test=0;  //index of the choosen test
  int nb_test=3;  //number of example tests

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
    printf( "test_mpi_stmm_example [testnumber] [w]\n" );
    exit(1);
  }



  int i,j; //some indexes
  int n; //row dimension of the toeplitz matrix 
  int lambda; //half bandwith size
  int m; //column dimension
 
  //MPI parameters
  int rank, size;
  MPI_Init(&argc, &argv);                    //Initialise MPI 
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );    //Id rank of the process    
  MPI_Comm_size( MPI_COMM_WORLD, &size );    //Number of processes     
  int *displs = (int*) calloc(size,sizeof(int)); //indexes to scatter over the processes
  int *nranks = (int*) calloc(size,sizeof(int)); //size of each process



//test case 1: (n = 5000, m = 1 and lambda = 50; data dispatched fairly between the processes)
//This is showing a basic example of how to scatter the data from a vector fairly between
//the processes. The last process contains the bigger data size due to the remainder of the
//euclidean division. 
  if(test==1) {
    if(rank==0) {
    printf("============================\n");
    printf("test 1\n");
    printf("----------------------------\n");
    }

    //Parameters
    n = 5000;//10;//5000; //row dimension of the toeplitz matrix
    lambda = 50;//3;//50; //half bandwith size 
    m = 1; //column dimension

    int nrank  = n*m/size;
    for(i=0;i<size-1;i++){
      displs[i]=i*nrank; //indexes to scatter over the processes
      nranks[i]=nrank;   //size of each process
    }
    displs[size-1] = displs[size-2]+nranks[size-2]; //last process index
    nranks[size-1] = m*n-displs[size-1]; //size of the last process

    if(rank==0) {
      printf("Parameters setup...\n");
      printf("n = %d, m = %d and lambda = %d\n", n, m, lambda);}
    test_toeplitz(n, m, lambda, displs, nranks);
    }


//test case 2: (n = 5000, m = 2 and lambda = 50; data dispatched fairly between the processes)
//This is showing a basic example of how to scatter the data from two vector (as a matrix n
//by 2) fairly between the processes. The last process contains the bigger data size due to
//the remainder of the euclidean division. Just the m parameter change from the test case 1.
  if(test==2) {
    if(rank==0) {
    printf("============================\n");
    printf("test 2\n");
    printf("----------------------------\n");
    }

    //Parameters
    n = 5000; //row dimension of the toeplitz matrix
    lambda = 50; //half bandwith size
    m = 2; //column dimension
    int nrank  = n*m/size;
    for(i=0;i<size-1;i++){
      displs[i]=i*nrank; //indexes to scatter over the processes
      nranks[i]=nrank; //size of each process
    }
    displs[size-1] = displs[size-2]+nranks[size-2]; //last process index
    nranks[size-1] = m*n-displs[size-1]; //size of the last process

    if(rank==0) {
      printf("Parameters setup...\n");
      printf("n = %d, m = %d and lambda = %d\n", n, m, lambda); }
    test_toeplitz(n, m, lambda, displs, nranks);
    }


//test case 3: (n = 5000, m = 5 and lambda = 50; data dispatched with variables sizes)
//This is showing an advance example of how to scatter the data from five vectors
//(as a matrix n by 5) dispatched with variables sizes. The last process contains the
//bigger data size due to the remainder of the euclidean division. The first process
//take lambda values, the second 3*n-lambda, the third n (all the fourth column) and
//the data left are dispatched fairly between the others processes. The last process
//still contains the bigger data size due to the remainder of the euclidean division.
  if(test==3) {
    if(rank==0) {
    printf("============================\n");
    printf("test 3\n");
    printf("----------------------------\n");
    }

    //Parameters
    n = 5000; //row dimension of the toeplitz matrix
    lambda = 50; //half bandwith size
    m = 5; //column dimension

    if (size<5){
      printf("Error : this test uses 5 MPI processes at least and you used %d only.\n", size);
      return -1;}

//Indexes to scatter over the processes
    int l0 = lambda;  //size of the first process
    int l1 = 2*lambda; //size of the last process
    int lc = 3*n; //cumul size of first and second processes 
    int r  = n*m-l0-(lc-l0)-n-l1; //size of left data
    displs[0] = 0;  //index for the first process
    nranks[0] = l0; //number of values for the first process
    displs[1] = displs[0] + nranks[0];  //index for the second process
    nranks[1] = lc-l0; //number of values for the second process
    displs[2] = lc; //index of the third process (first element of the fourth column)
    nranks[2] = n; //number of values for the third process

    //the data left are dispatched fairly between the others processes
    int nr = size-4;
    for (i=3; i<size-1; i++){
      displs[i] = displs[i-1]+nranks[i-1]; //indexes to scatter over the processes
      nranks[i] = r/nr;} //size of each process

    //last process
    displs[size-1] = displs[size-2]+nranks[size-2]; //index for the last process 
    nranks[size-1] = m*n-displs[size-1];  //size of the last process

    if(rank==0) {
      printf("Parameters setup...\n");
      printf("n = %d, m = %d and lambda = %d\n", n, m, lambda); }
    test_toeplitz(n, m, lambda, displs, nranks);
    }

 
  free(displs);
  free(nranks);
  MPI_Finalize();

  return 0;
}

 
// Perform the multiplication of a Toeplitz matrix by a matrix with MPI using the
//defined parameters. The input matrix V is construct on the process 0 and then 
//scattered to all the processes for the computation. We compare the obtain result
//to the exact solution and print the relative euclidiean norm residual. For this
//example, V is filled with 1 everywhere and T composed by a constant value
//on the bandwith so that the sum over the all bandwith is equal to 1.
int test_toeplitz(int n, int m, int lambda, int* displs, int* nranks)
{

  double *T;  //input toeplitz
  double *V, *Vrank;  //input matrix
  double *TV;//output matrix
  int i,j; //some indexes

  //MPI parameters
  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );    
  MPI_Comm_size( MPI_COMM_WORLD, &size );     
  
  int nrank  = n*m/size;

//allocation of T
  T  = (double *) calloc(lambda,sizeof(double));
 
  if(rank==0) {
    for(i=0;i<size;i++)
      printf("Rank %d : index %d ; size %d\n", i, displs[i], nranks[i]);

    //input matrix definition for T  
    for(i=0;i<lambda;i++) // half bandwith Toeplitz matrix
      T[i]= 1./(2*lambda-1);  
  }

//spread T over the processes  
  MPI_Bcast(T, lambda, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(rank==0) {
//allocation for V and TV
    V  = (double *) calloc(n*m   ,sizeof(double));
    TV = (double *) calloc(n*m   ,sizeof(double));

    //input matrix definition of V    
    for(j=0;j<m;j++) 
      for(i=0;i<n;i++)
    	V[j*n+i] = 1.;  

  }
 
  //Compute max(nranks) for imformation
  int maxsize = nranks[0];
  for (i=1; i<size; i++)
    if (nranks[i]>maxsize)
      maxsize=nranks[i];


//Product computation
  //allocations of the local V
  Vrank = (double *) calloc(nranks[rank], sizeof(double));


  //times variables
  double t1,t2;

  //scatter V over the processes
  MPI_Scatterv(V, nranks, displs, MPI_DOUBLE, Vrank, nranks[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);


  if (rank==0)
    printf("Product computation...\n");
  MPI_Barrier(MPI_COMM_WORLD);

  t1=MPI_Wtime();
  //perform local product using midas DA lib
  mpi_stmm_old(&Vrank, n, m, displs[rank], nranks[rank], T, lambda, MPI_COMM_WORLD);
  t2=  MPI_Wtime();

  //receive data from each process
  MPI_Gatherv(Vrank, nranks[rank], MPI_DOUBLE, TV, nranks, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

  if (rank==0)
    printf("Results...\n");

  MPI_Barrier(MPI_COMM_WORLD);
  printf("[rank %d] Computation time : %f s.\n", rank, t2-t1);

//compute the relative residual error for each process
  if (rank==0) {
    for (j=0; j<size; j++){  //loop on the process
      int nb_false=0;

      double resloc;
      double res2=0;
      double norm2relres=0;
      double sol;
      double norm2sol=0;
      double res=0;


      int ibis;
      for (i=displs[j]; i<displs[j]+nranks[j]; i++){ //loop on each element of the process
        ibis = i%(n); //row index on the actual column

        if (ibis<lambda) {
          sol = (lambda+ibis)*1./(2*lambda-1) ; //exact solution for the lambda first elements
          resloc = fabs(TV[i]- sol);        }   //residual error
        else if (ibis>(n-1-lambda)) {
          sol = (lambda+(n-1)-ibis)*1./(2*lambda-1) ; //exact solution for the lambda last elements
          resloc = fabs(TV[i]- sol);        } //residual error
        else {
          sol = 1. ; //exact solution for the middle elements
          resloc = fabs(TV[i]- sol);        }  //residual error

        if(resloc>1e-8){ //number of false values to count the local error
          nb_false += 1; }

      res2 = res2 + resloc*resloc ;
      norm2sol = norm2sol + sol*sol;
      }//end of the loop processes

      res2 = sqrt(res2); //euclidean norm of the residual error
      norm2sol =  sqrt(norm2sol); //euclidean norm of the solution
      norm2relres = res2/norm2sol; //euclidean norm of the relative residual error
 
      if (fabs(norm2relres)<1e-8) {
	printf("Rank %d : Success !\n", j);  }
      else {
	printf("Rank %d : failed...\n", j);
	printf("Rank %d : Number of false values : %d over %d\n",j, nb_false, nranks[j]); }

      printf("Rank %d : norm2 residual is %e\n", j, norm2relres);
  }//end of column loop
 }//end of rank==0


//write to file if "write mode on" 
  if (CHAR_RW=='w') {

  //fich:
  FILE* file;
  char filename [1024];
  sprintf(filename,"output-test.txt");
  file = fopen(filename, "w");


  if(rank==0){
    fprintf(file,"Test result:\n");
    fprintf(file,"n=%d, m=%d\n", n, m);
    fprintf(file,"TV=\n");

    for(i=0;i<(n*m);i++)
      fprintf(file,"%lf\n", TV[i]);

  fclose(file);
  }

    printf("Write done...\n");
  }//end  CHAR_RW=='w'


    free(T);
    if (rank==0){
      free(Vrank);
      free(TV);
      free(V);
    }
   
    return 0;
  }


