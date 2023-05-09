/** @file   test_diags_mat_0.c
    - create a matrix with cyclic diaganals, A.
    - create a vector b.
    - run conjugate gradient for the system A'*A*x=A'*b
    @author Pierre Cargemel
    @date   Septemeber 2012*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "midapack.h"

extern char *optarg;

void usage(){ 
    printf("\n*#  Usage: mpirun -n [nb procs] rand_mapmat -m [nb local rows]");
    printf("\n*#  -m [int] set the number of local rows. Global number of rows will be m times nb procs");
    printf("\n*#  -N [int] set the maxmum number of columns (column indices will belong [0, N-1])");
    printf("\n*#  -Z [int] set the number non-zero values per column");
    printf("\n*#  -f [int] special flag for the communication scheme");
    printf("\n*#  -K [int] maximum iterate");
    printf("\n*#  -r [int] residual tolerance\n");
}


int main(int argc, char *argv[]){
  int		M, N, Nnz;	//global rows and columns, non-zero per column
  int		m, n;		//local rows and columns
  int		err, i, j, k, r;	
  Mat		A;		//matrix struct
  int 		*indices;	//list of column indices
  double 	*values;	//list of non-zero values
  int 		flag;		//communication flag, sort method
  double	*x, *b;		//solution, right hand side 
  double	lr, norm;	 
  char 		ch, *endp;
  int 		K;
  double 	rtol;
  int 		rank, size;	
  
  /*mpi init*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  M=1;
  N=1;
  flag=1;
  Nnz=1;
  K=1;
  rtol=0.01;
  if(argc==0){
      usage();
  }

  /*set options*/
  while((ch = getopt( argc, argv, "m:N:Z:f:K:r:" )) != EOF){
    switch(ch) {
      case 'N':
 	N = strtol(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && N <= 0 )){
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'm':
	m = strtol(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && m <= 0 )){
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'Z':
        Nnz = strtol(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && Nnz <= 0 )){
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'f':
 	flag = strtol(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (flag < 0 || flag > 5))){
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'K':
 	K = strtol(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && K<=0)){
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'r':
        rtol = strtod(optarg, &endp);
	if ((errno == ERANGE) || (errno != 0 && rtol <= 0)){
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
    }
  }
  M=m*size;							//global rows
 
  values  = (double *) malloc(m*Nnz * sizeof(double));		//set non-zero values
  for(j=0; j<m; j++){					//
    for(k=0; k<Nnz; k++)					//
      values[j*Nnz+k]=1.0/(k+1);						//
  }     		

  indices = (int *) malloc(m*Nnz * sizeof(int));		//set columns indices
  for(j=0; j<m; j++){					//
    for(k=0; k<Nnz; k++)					//
      indices[j*Nnz+k] = ((rank*m)+j+k)%N;
  }     		
  for(i=0; i<m*Nnz; i++)                      //assemblage
 //   printf("(%d,%d,%lf) ",rank,indices[i],values[i]);
 // printf("\n ");
 						//
  MatInit(&A, m, Nnz, indices, values, flag, MPI_COMM_WORLD); 	//init matrix


  x = (double *) malloc(A.lcount*sizeof(double));		//allocate pixel domain vector 
  for(j=0; j<A.lcount; j++)						//set values
    x[j] = 0.0;							//
  b = (double *) malloc(m*sizeof(double));    			//allocate time domain vector
  for(i=0; i<m; i++)						//set values
    b[i] = 1.0;							//

  MatInfo(&A, 0, "Diags");					//save matrix info in Rand_info.txt
  
  CGNE(A, x, b, rtol, K);
   
  MatFree(&A);							//free memory  
  free(b);							//	
  free(x);							//
  MPI_Finalize();
  return 0;
}


