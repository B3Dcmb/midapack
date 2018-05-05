#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "toeplitz.h"
#include <time.h>

// Using examples for the symmetric Toeplitz block-diagonal matrix-matrix product routine (mpi_stbmm).
// Test cases examples are showning using differents parameters.
// Test number as first argument
// -1 (default) : all test will be run. 
// 1 : nblocks = 7, m = 1, n[i] = 500 and lambda = [10 51 159 73 102 55 49] ; 
// 2 : nblocks = 4, m = 2, n[i] = 100 and lambda[i] = 3 ;
// 3:  nblocks = 2, m = 1, n[i] = 10, lambda[i] = 3 and lv[i]=20 

char CHAR_RW='\0';  //global variable for write mode

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
    printf( "test_mpi_stbmm_example [testnumber] [w]\n" );
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
  int *displs = (int*) calloc(size,sizeof(int));  //indexes to scatter over the processes
  int *nranks = (int*) calloc(size,sizeof(int));  //size of each process 

  
  int toto;            /* 16-bit integer          */
  //unsigned int toto2;
  int64_t var_int64_t;
  int64_t var_int64_t_id0;
  double var_dble;
  double toto3;
  int toto4;

  toto=1000*1000*1000*2;
  //toto2=1000*1000*1000*4;
  toto3=1000LL*1000LL*1000LL*8LL;
  toto4 = toto3/4;
  printf("toto=%d\n", toto);
  //printf("toto2=%u\n", toto2);
  printf("toto3=%lld\n", toto3);
  printf("toto4=%lld\n", toto4);


  var_dble = 1234567891011;

//Affecte la valeur 4294967296 dans a
long int a = 4294967296LL;
//Affiche cette valeur
printf("aaaaaa=%lld\n", a);

long int b = 1000000000000LL;
long int c = 1000000000000LL;

  var_int64_t = 8000000000LL;
  var_int64_t_id0 = 8000000003LL;

printf("var_int64_t_id0=%lld\n", var_int64_t_id0);
printf("var_int64_t=%lld\n", var_int64_t);
printf("modulo=%lld\n", var_int64_t_id0%var_int64_t);
printf("modulo=%lld\n", var_int64_t%var_int64_t_id0);

//Affiche cette valeur
printf("b=%lld\n", b);
printf("c=%lld\n", c);

printf("b=%ld\n", b);
printf("c=%ld\n", c);

//int b = 2000000000;
//int c = 2000000000;
long int d;
d = b+c;
printf("d=%ld\n", d);


return 0;

}
