/** @file   test_spI_mapmat.c
  
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
#include <midapack.h>


extern char *optarg;

void usage(){ 
    printf("\n*#  Usage: mpirun -n [number of processes] test_spI_mapmat [option[value]] ");
    printf("\n*#  -M [int] set the global number of rows. Rows are uniformly distributed over processes");
    printf("\n*#  -Z [int] set the number non-zero values per column (default = 1)");
    printf("\n*#  -f [int] special flag for the communication scheme (default = 1)");
    printf("\n*#  -s [int] special flag for the sorting method");
    printf("\n*#  -p [int] special flag for mutli-thread function");
    printf("\n*#  -o [filename] specify an output filename");
    printf("\n*#  -t use timers and print execution to benchmark"); 
    printf("\n*#  -i print matrix info"); 
}

int check(int err){
  int rank;
  int gerr;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Reduce(&err, &gerr, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(gerr==0){
    if(rank==0)
    printf("COMPLETED\n");
  }
  else{
    if(rank==0)
      printf("FAILED err %d\n", gerr);
//    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  return gerr;
}

int main(int argc, char *argv[]){
  int		M, N, Nnz;		//global rows and columns, non-zero per column
  int		m, n;			//local rows and columns
  int 		gif;			//global index of the first local column
  int		err, i, j, k;	
  Mat	A;			//matrix struct
  int 		*indices;
  double 	*values;
  int 		flag, sort, omp ;		//communication flag, sort method
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
 
  flag=1;
  sort=1;
  omp=0;
  output=0;
  timer=0;
  info=0; 
  err=0;
  Nnz=1;
  if(argc==0){
      usage();
  }

  /*set options*/
  while((ch = getopt( argc, argv, "M:Z:f:s:p:o:ti" )) != EOF){
    switch(ch) {
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
      case 'p':
 	omp = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (omp < 0 || omp > 3))){
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


  N=M;

  partition(&gif, &m, M, rank, size);                     

  MatAllocate(&A, m, Nnz, MPI_COMM_WORLD);	//create matrix

  indices  = (int *) malloc(m * sizeof(int));			//assembly/
  values  = (double *) malloc(m * sizeof(double));


  for(i=0; i<m; i++){
    values[i]=1.0;
    indices[i] = gif+i;
  }  
  MatSetIndices(&A, m, 0, Nnz, indices);
  MatSetValues(&A, m, 0, Nnz, values);

  srand(gif*(Nnz-1)); //init seed
  for(j=1; j<Nnz; j++){
    for(i=0; i<m; i++){                      
      values[i] = 0.0;
      indices[i] = rand()%N;  
    }
    MatSetIndices(&A, m, j, Nnz, indices);
    MatSetValues(&A, m, j, Nnz, values);
  }
  free(indices);
  free(values);
  

  if(output==1)
    MatSave(&A, filename);				//matrix save

  k=0;
  MPI_Barrier(MPI_COMM_WORLD);                                
  st = MPI_Wtime();
  MatLocalShape(&A, sort, omp);
  t[k] = MPI_Wtime()-st;

  x   = (double *) malloc(A.lcount*sizeof(double));		// allocate pixel domain vector 
  for(j=0; j<A.lcount; j++)
    x[j] = 1.0;
  y   = (double *) malloc(m*sizeof(double));    			// allocate time domain vector
  for(i=0; i<m; i++)
    y[i] = 1.0;
  
  k++;
  MPI_Barrier(MPI_COMM_WORLD);                                
  st = MPI_Wtime();
  TrMatVecProd_Naive(&A, y, x, omp);	//trmatvec_naive
  t[k] = MPI_Wtime()-st;
  for(j=0; j<A.lcount; j++){			//little check
    if(x[j]!=1.0){
      err++;
    }
  }
  check(err);
  for(j=0; j<A.lcount; j++)			//little check
    x[j]=1.0;
 
  k++;
  MPI_Barrier(MPI_COMM_WORLD);  
  st = MPI_Wtime();
  err = MatComShape(&A, flag);         		//com_shape/
  t[k] = MPI_Wtime()-st;

  if(info==1)
    MatInfo(&A, 0, "SpI");

  k++;
  MPI_Barrier(MPI_COMM_WORLD);                 
  st = MPI_Wtime();
  err += TrMatVecProd(&A, y, x, 0);		//trmatvec
  t[k] = MPI_Wtime()-st;
  /*if(rank==0 || rank ==3 || rank==6 || rank ==7){
    printf("\nid Rank %d :", rank);
    for(j=0; j<A.lcount; j++){			//little check
      printf(" %d ",  A.lindices[j]);
    }
    printf("\nid Rank %d :", rank);
    for(j=0; j<A.com_count; j++){			//little check
      printf(" %d ",  A.com_indices[j]);
    }
  }*/
    
  for(j=0; j<A.lcount; j++){			//little check
    if(x[j]!=1.0){
      printf(" (%d: %d %lf) ", rank, A.lindices[j], x[j]);
      err++;
    }
  }
  check(err);
  for(j=0; j<A.lcount; j++)			//little check
    x[j]=1.0;

  k++;
  MPI_Barrier(MPI_COMM_WORLD);                 
  st = MPI_Wtime();
  MatVecProd(&A, x, y, omp);//matvec
  t[k] = MPI_Wtime()-st;

  for(i=0; i<m; i++){			//little check
    if(y[i]!=1.0){
      printf(" (%d: %d %lf) ", rank, i, y[i]);
     err++;
    }
  }
  check(err);
  
  if(rank==0 && timer==1){
    printf("\nFunction\ttime(sec)");
    printf("\nLocalShape\t%lf\nTrMatVec_Naive\t%lf\nComShape\t%lf\nTrMatVec\t%lf\nMatVec\t%lf\n", t[0], t[1], t[2], t[3], t[4]);
  }
 
  MatFree(&A);                                                //free memory  
  free(y);
  free(x);
  MPI_Finalize();
  
  return 0;
}


