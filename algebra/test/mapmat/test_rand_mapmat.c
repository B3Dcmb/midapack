/** @file   test_rand_mapmat.c
  
    @author Pierre Cargemel
    @date   December 2011 */

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
    printf("\n*#  Usage: mpirun -n [number of processes] test_mat [option[value]] ");
    printf("\n*#  -M [int] set the global number of rows. Rows are uniformly distributed over processes");
    printf("\n*#  -N [int] set the global number of columns. ");
    printf("\n*#  -Z [int] set the number non-zero values per column (default = 1)");
    printf("\n*#  -f [int] special flag for the communication scheme (default = 1)");
    printf("\n*#  -o [filename] specify an output filename");
    printf("\n*#  -t use timers and print execution to benchmark"); 
    printf("\n*#  -i print matrix info"); 
}

int check(int err){
  if(err==0)
    printf("COMPLETED\n");
  return err;
}

int main(int argc, char *argv[]){
  int		M, N, Nnz;		//global rows and columns, non-zero per column
  int		m, n;			//local rows and columns
  int 		gif;			//global index of the first local column
  int		err, i, j, k;	
  Mat	A;			//matrix struct
  int 		*indices;
  double 	*values;
  int 		flag, sort;		//communication flag, sort method
  double	*y, *y_f, *x, *x_n;	//time domain vector, pixel domain vector 
  double	t[5], st;       	//timer table, start time
  char 		ch;
  char		*filename;		//output filename
  int 		output, timer, info;         
  int 		rank, size;	
  
  /*mpi init*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
 
  output=0;
  timer=0;
  info=0; 
  err=0;
  Nnz=1;
  if(argc==0){
      usage();
  }

  /*set options*/
  while((ch = getopt( argc, argv, "M:N:Z:f:s:o:ti" )) != EOF){
    switch(ch) {
      case 'N':
 	N = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (N <= 0 || N >= 2147483648))){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'M':
	M = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (M <= 0 || M >= 2147483648))){
          err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'Z':
        Nnz = strtol(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && (Nnz <= 0 || Nnz >= 2147483648))){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'f':
 	flag = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (flag < 0 || flag > 5))){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 's':
 	sort = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (sort < 0 || sort > 5))){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'o':
        filename = strdup(optarg);
        output=1;
      break;
      case 'i':
        info=1;
      break;
      case 't':
        timer=1;
      break;
    }
  }

  partition(&gif, &m, M, rank, size);                     

  err = MatAllocate(&A, m, Nnz, MPI_COMM_WORLD);	//create matrix

  indices  = (int *) malloc(Nnz * sizeof(int));			//assembly/
  values  = (double *) malloc(Nnz * sizeof(double));
  for(j=0; j<Nnz; j++)
    values[j]=1.0;
  srand(gif*Nnz); //init seed
  for(i=0; i<m; i++){                      
    for(j=0; j<Nnz; j++)
      indices[j] = rand()%N;  
    MatSetIndices(&A, Nnz, i*Nnz, 1, indices);
    MatSetValues(&A, Nnz, i*Nnz, 1, values);
  }

  if(output==1)
    MatSave(&A, filename);				//matrix save

  k=0;
  MPI_Barrier(MPI_COMM_WORLD);                                
  st = MPI_Wtime();
  err = MatLocalShape(&A, sort, 0);
  t[k] = MPI_Wtime()-st;

  x_n = (double *) malloc(A.lcount*sizeof(double));		// allocate pixel domain vector 
  x   = (double *) malloc(A.lcount*sizeof(double));		// allocate pixel domain vector 
  y   = (double *) malloc(m*sizeof(double));    			// allocate time domain vector
  y_f = (double *) malloc(m*sizeof(double));    			// allocate time domain vector
  for(i=0; i<m; i++)
    y[i] = 1.0;
  
  k++;
  MPI_Barrier(MPI_COMM_WORLD);                                
  st = MPI_Wtime();
  err = TrMatVecProd_Naive(&A, y, x_n, 0);	//trmatvec_naive
  t[k] = MPI_Wtime()-st;

  k++;
  MPI_Barrier(MPI_COMM_WORLD);  
  st = MPI_Wtime();
  err = MatComShape(&A, flag);         		//com_shape/
  t[k] = MPI_Wtime()-st;

  if(info==1)
    MatInfo(&A, 0, "Rand");

  k++;
  MPI_Barrier(MPI_COMM_WORLD);                 
  st = MPI_Wtime();
  err = TrMatVecProd(&A, y, x, 0);		//trmatvec
  t[k] = MPI_Wtime()-st;

  for(j=0; j<A.lcount; j++){			//little check
    if(x_n[j]!=x[j])
      err++;
    x[j]=1/x[j];
  }
  //if(err !=0)
    //MPI_Abort(MPI_COMM_WORLD, 1);
 
  k++;
  MPI_Barrier(MPI_COMM_WORLD);                 
  st = MPI_Wtime();
  err = MatVecProd(&A, x, y_f, 0);//matvec
  t[k] = MPI_Wtime()-st;

//  for(i=0; i<m; i++){			//little check
 //   if(abs(y[i]-y_f[i])<0.01){
  //    printf("(%lf %lf",y[i],y_f[i]);
  //    err++;
  //  }
 // }
  //if(err !=0)
    //MPI_Abort(MPI_COMM_WORLD, 1);
  
  if(rank==0 && timer==1){
    printf("\nFunction\ttime(sec)");
    printf("\nLocalShape\t%lf\nTrMatVec_Naive\t%lf\nComShape\t%lf\nTrMatVec\t%lf\nMatVec\t%lf\n", t[0], t[1], t[2], t[3], t[4]);
  }
 
  MatFree(&A);                                                //free memory  
  free(y);
  free(y_f);
  free(x);
  free(x_n);
  free(indices);
  free(values);
  MPI_Finalize();
  
  return check(err);
}


