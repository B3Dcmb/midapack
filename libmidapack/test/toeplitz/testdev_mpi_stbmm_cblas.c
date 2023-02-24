#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "toeplitz.h"
#include <cblas.h>
#include <time.h>

// Perform unit testing of the MPI Bloc-diagonal routine
// Test number as first argument
// -1 (default) : all test will be run. 
// 0 : nblocs = 7, m = 1, n[i] = 500 and lambda = [10 51 159 73 102 55 49]; 
// 1 : nblocs = 1, m = 1, n = 500 and lambda = 51; 
// 2 : nblocs = 1, m = 6, n = 500 and lambda = 51; 
// 3 : nblocs = 7, m = 6, n[i] = 500 and lambda = [10 51 159 73 102 55 49]; 


 extern int NFFT;

int main (int argc, char *argv[])
{


  int test = -1;
  if (argc>1)
    test = atoi(argv[1]);

  int i,j;

  int rank, size;
  
 
  // MPI parameters
  MPI_Init(&argc, &argv);                    // Initialise MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );    // Récupère l'id du proc
  MPI_Comm_size( MPI_COMM_WORLD, &size );    // Récupère le nombre de 
  int *displs = (int*) calloc(size,sizeof(int));
  int *nranks = (int*) calloc(size,sizeof(int)); 


  int *n, m, *lambda, nb_blocs, nb_col_blocs;
  int *mblocks;



  if(test==0 || test < 0){
    if(rank==0)
      printf("test 0\n");


//    NFFT=2;

    nb_blocs = 4;
    m = 1;
    n      = (int *) calloc(nb_blocs, sizeof(int));
    lambda = (int *) calloc(nb_blocs, sizeof(int));
    mblocks = (int *) calloc(nb_blocs, sizeof(int));

    for(i=0;i<nb_blocs;i++)
      n[i] = 10;

    for (i=0;i<nb_blocs;i++)
      mblocks[i] = 1;

    lambda[0] = 3;
    lambda[1] = 3;//4;
    lambda[2] = 3;
    lambda[3] = 3;

    int nsample=0;
    for(j=0;j<nb_blocs;j++)
      nsample += (n[j]*mblocks[j]);

    if (rank==0)
      printf("nsample=%d\n", nsample);

    int nrank  = (nsample*m)/size;

    for(i=0;i<size;i++){
      displs[i]=i*nrank;
      nranks[i]=nrank;
    }

    displs[size-1] = displs[size-2]+nranks[size-2];
    nranks[size-1] = nsample*m-displs[size-1];

    if (rank==0) {
    for(i=0;i<size;i++){
      printf("displs[i=%d]=%d\n", i, displs[i]);
      printf("nranks[i=%d]=%d\n", i, nranks[i]);
    }}
    test_toeplitz(n, m, lambda, nb_blocs, displs, nranks, mblocks);
    }


  if(test==1 || test < 0){
    if(rank==0)
      printf("test 1\n");
    nb_blocs = 2;
    m = 2;
    n      = (int *) calloc(nb_blocs, sizeof(int));
    lambda = (int *) calloc(nb_blocs, sizeof(int));
    mblocks = (int *) calloc(nb_blocs, sizeof(int));

    for(i=0;i<nb_blocs;i++)
      n[i] = 10;

    for (i=0;i<nb_blocs;i++)
      mblocks[i] = 1;

    lambda[0] = 3;
    lambda[1] = 3;//4;

    int nsample=0;
    for(j=0;j<nb_blocs;j++)
      nsample += (n[j]*mblocks[j]);

    if (rank==0)
      printf("nsample=%d\n", nsample);

    int nrank  = (nsample*m)/size;

    for(i=0;i<size;i++){
      displs[i]=i*nrank;
      nranks[i]=nrank;
    }

    displs[size-1] = displs[size-2]+nranks[size-2];
    nranks[size-1] = nsample*m-displs[size-1];

    if (rank==0) {
    for(i=0;i<size;i++){
      printf("displs[i=%d]=%d\n", i, displs[i]);
      printf("nranks[i=%d]=%d\n", i, nranks[i]);
    }}
    test_toeplitz(n, m, lambda, nb_blocs, displs, nranks, mblocks);
    }



  if(test==2 || test < 0){
    if(rank==0)
      printf("test 1\n");
    nb_blocs = 3;
    m = 1;
    n      = (int *) calloc(nb_blocs, sizeof(int));
    lambda = (int *) calloc(nb_blocs, sizeof(int));
    mblocks = (int *) calloc(nb_blocs, sizeof(int));

    for(i=0;i<nb_blocs;i++)
      n[i] = 10;

    for (i=0;i<nb_blocs;i++)
      mblocks[i] = 1;

    lambda[0] = 3;
    lambda[1] = 3;//4;
    lambda[2] = 3;

    int nsample=0;
    for(j=0;j<nb_blocs;j++)
      nsample += (n[j]*mblocks[j]);

    if (rank==0)
      printf("nsample=%d\n", nsample);

    int nrank  = (nsample*m)/size;

    for(i=0;i<size;i++){
      displs[i]=i*nrank;
      nranks[i]=nrank;
    }

    displs[size-1] = displs[size-2]+nranks[size-2];
    nranks[size-1] = nsample*m-displs[size-1];

    if (rank==0) {
    for(i=0;i<size;i++){
      printf("displs[i=%d]=%d\n", i, displs[i]);
      printf("nranks[i=%d]=%d\n", i, nranks[i]);
    }}
    test_toeplitz(n, m, lambda, nb_blocs, displs, nranks, mblocks);
    }






  if(test==3 || test < 0){
    if(rank==0)
      printf("test 1\n");
    nb_blocs = 2;
    m = 2;
    n      = (int *) calloc(nb_blocs, sizeof(int));
    lambda = (int *) calloc(nb_blocs, sizeof(int));
    mblocks = (int *) calloc(nb_blocs, sizeof(int));

    for(i=0;i<nb_blocs;i++)
      n[i] = 10;

    for (i=0;i<nb_blocs;i++)
      mblocks[i] = 1;

    lambda[0] = 3;
    lambda[1] = 3;//4;

    int nsample=0;
    for(j=0;j<nb_blocs;j++)
      nsample += (n[j]*mblocks[j]);

    if (rank==0)
      printf("nsample=%d\n", nsample);

    int nrank  = (nsample*m)/size;

    for(i=0;i<size;i++){
      displs[i]=i*nrank;
      nranks[i]=nrank;
    }

    displs[size-1] = displs[size-2]+nranks[size-2];
    nranks[size-1] = nsample*m-displs[size-1];

    if (rank==0) {
    for(i=0;i<size;i++){
      printf("displs[i=%d]=%d\n", i, displs[i]);
      printf("nranks[i=%d]=%d\n", i, nranks[i]);
    }}
    test_toeplitz(n, m, lambda, nb_blocs, displs, nranks, mblocks);
    }




  if(test==10 || test < 0){
    if(rank==0) {
    printf("============================\n");
    printf("test 1\n");
    printf("----------------------------\n");
    }

// parameters
    if(rank==0)
      printf("Parameters setup...\n");

    nb_blocs = 7; //number of Toeplitz blocks on the diagonal of the matrix
    m = 2; //global column dimension
    n      = (int *) calloc(nb_blocs, sizeof(int));
    lambda = (int *) calloc(nb_blocs, sizeof(int));
    mblocks = (int *) calloc(nb_blocs, sizeof(int));

    for(i=0;i<nb_blocs;i++)
      n[i] = 500;   //row dimension of each toeplitz block

    for (i=0;i<nb_blocs;i++)
      mblocks[i] = 1;  //column dimension of each toeplitz block

    //half bandwith size for each toeplitz block
    lambda[0] = 10;
    lambda[1] = 51;
    lambda[2] = 159;
    lambda[3] = 73;
    lambda[4] = 102;
    lambda[5] = 55;
    lambda[6] = 49;

    int nsample=0; //global row dimension of the matrix
    for(j=0;j<nb_blocs;j++)
      nsample += (n[j]*mblocks[j]);

    if (rank==0) {
      printf("global n = %d\n", nsample);
      printf("global m = %d\n", m);
    }

    int nrank  = (nsample*m)/size;
    for(i=0;i<size-1;i++){
      displs[i]=i*nrank;  //indexes to scatter over the processes
      nranks[i]=nrank;  //size of each process
    }
    displs[size-1] = displs[size-2]+nranks[size-2]; //last process index
    nranks[size-1] = nsample*m-displs[size-1]; //size of the last process

    test_toeplitz(n, m, lambda, nb_blocs, displs, nranks, mblocks);
  } //End test==1



  free(displs);
  free(nranks);
  MPI_Finalize();

  return 0;
}



// Compare middle level routine output to Cblas
int test_toeplitz(int *n, int m, int *lambda, int nb_blocs, int *displs, int *nranks, int *mblocks)
{
  time_t start,end;
  double dif;

  double *T, *T2;  // input toeplitz
  double *V;       // input matrix
  double *TV, *TV2;// output matrix
  
  int rank,size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );    // Récupère l'id du proc
  MPI_Comm_size( MPI_COMM_WORLD, &size );    // Récupère le nombre de 


  FILE *file; 
  file = stdout;

  int i,j,k,l,p; // some index
  int ntot = 0;
  int lambdatot = 0;
  int lvtot = 0;

  for(i=0;i<nb_blocs;i++) {
    lambdatot += lambda[i];
    ntot += n[i];
    lvtot += mblocks[i]*n[i];
   }

  int m_rowwise=2;

  int *dispn = (int *) calloc(nb_blocs, sizeof(int));
  int *displ = (int *) calloc(nb_blocs, sizeof(int));

  dispn[0]=0;
  displ[0]=0;
  for(i=1;i<nb_blocs;i++) {
    dispn[i] = dispn[i-1] + n[i-1]*mblocks[i-1];
    displ[i] = displ[i-1] + lambda[i-1]; }

  // alloc arrays
 
// define mbloc as the max of the mblocks for convenience 
  int mbloc = mblocks[0];
  for(i=1;i<nb_blocs;i++) 
    if (mblocks[i]>mbloc)
     mbloc=mblocks[i];



 
  T  = (double *) calloc(lambdatot,sizeof(double));
  if(rank==0) {
    T2 = (double *) calloc(ntot*ntot,sizeof(double));
    V  = (double *) calloc(mbloc*ntot*m*m_rowwise   ,sizeof(double));  
    TV = (double *) calloc(mbloc*ntot*m*m_rowwise   ,sizeof(double));  
    TV2= (double *) calloc(mbloc*ntot*m*m_rowwise   ,sizeof(double));  
    
    
    // init arrays
//    for(i=0;i<lambdatot;i++) // Toeplitz matrix (band only)
//      T[i]=rand()/((double) RAND_MAX); 
    

  for(i=0;i<lambdatot;i++) {// Toeplitz matrix (band only)
    if (i%3 == 0) {
//      T[i]= rand()/((double) RAND_MAX);}
      T[i]=10;}//rand()/((double) RAND_MAX);
    else if (i%3 == 1) {
//      T[i]= rand()/((double) RAND_MAX);}
      T[i]=2;}
    else if (i%3 == 2) {
//      T[i]=rand()/((double) RAND_MAX);}
      T[i]=3;}

    else {
      T[i]= 0.;//rand()/((double) RAND_MAX);}
   }}

    
    for (j=0;j<ntot;j++) // Full Toeplitz matrix needed for cblas computation
      for(i=0;i<ntot;i++)
	T2[j*ntot+i] = 0.0;

    

  printf("ntot=%d\n", ntot);
/*
    for(i=0;i<ntot*m*mbloc;i++) // input matrix
      V[i] = (double) ((i)%10 +1);


    for(i=0;i<10;i++) // input matrix
      V[i] = (double) (1);
*/
    for(i=0;i<lambdatot;i++) // input matrix
      printf("T[i=%d]=%f\n", i, T[i]);


    for(i=0;i<ntot*m*mbloc;i++) // input matrix
      V[i] = (i)%10+1;//rand()/((double) RAND_MAX);



//    for(i=0;i<ntot*m*mbloc;i++) // input matrix
//      printf("V[i]=%f\n", V[i]);


   }//End Rank0



   int *lv;
   lv = (int *) calloc(nb_blocs,sizeof(int));

    for(i=0;i<nb_blocs;i++) // input matrix
      lv[i] = mblocks[i]*n[i];

 if (rank==0) {
    for(i=0;i<nb_blocs;i++) {
      printf("mblocks[i=%d]=%d\n", i, mblocks[i]);
      printf("n[i=%d]=%d\n", i, n[i]);
      printf("lv[i=%d]=%d\n", i, lv[i]);
    }}


  MPI_Bcast(T, lambdatot, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    
  int maxsize = nranks[0];
  for (i=1; i<size; i++)
    if (nranks[i]>maxsize)
      maxsize=nranks[i];
  
  // MPI test
  double *Vrank = (double *) calloc(maxsize*m_rowwise, sizeof(double));
  double *TVrank = (double *) calloc(maxsize*m_rowwise,sizeof(double));
  
  int nrow = ntot;
  int id0 = displs[rank];
  int local_V_size=nranks[rank];

  
  MPI_Scatterv(V, nranks, displs, MPI_DOUBLE, Vrank, maxsize, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

    for(i=0;i<nranks[rank];i++) // input matrix
      fprintf(file, "[%d] Vrank[%d]=%f\n", rank, i, Vrank[i]);

      fprintf(file, "=========\n");

    for(i=0;i<nb_blocs;i++) 
      fprintf(file, "[%d] n[%d]=%d\n", rank, i, n[i]);
    fprintf(file, "[%d] m=%d, nrow=%d, nb_blocs=%d\n", rank, m, nrow, nb_blocs);
    fprintf(file, "[%d] id0=%d, local_V_size=%d\n", rank, id0, local_V_size);

    for(i=0;i<nb_blocs;i++)
      fprintf(file, "[%d] dispn[%d]=%d\n", rank, i, dispn[i]);


    for(i=0;i<nranks[rank];i++) // input matrix
      Vrank[i+local_V_size]= Vrank[i];


  mpi_stbmm(&Vrank, n, m, nrow, m_rowwise, T, nb_blocs, nb_blocs, lambda, dispn, id0, local_V_size, MPI_COMM_WORLD);

    fprintf(file, "==out=======\n");
  for(i=0;i<nranks[rank];i++) // output matrix
    fprintf(file, "[%d] Vrank[%d]=%f\n", rank, i, Vrank[i]);

  // receive data from each proc
  MPI_Gatherv(Vrank, nranks[rank], MPI_DOUBLE, TV, nranks, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
 
   
  // perform product using cblas routine
  // cblas_dgemm (CblasColMajor, CblasNoTrans,CblasNoTrans, n, m, n, 1, T2, n, V, n, 1, TV2, n);
  if(rank==0) {
    double res = 0.0;
    int nb_false = 0;
    for(l=0;l<m;l++) {
      for(k=0;k<nb_blocs;k++) {
	res = 0.0;
	nb_false = 0;

/*
	for (j=0;j<n[k];j++){ // Full Toeplitz matrix needed for cblas computation
	  for(i=0;i<n[k];i++)
	    T2[j*n[k]+i] = 0;
	  for (i=0;i<lambda[k];i++){
	    if (j-i>=0) 
	      T2[j*n[k]+j-i] = T[i+displ[k]]; 
	    if (j+i<n[k])
	    T2[j*n[k]+j+i] = T[i+displ[k]]; }}
*/	

  for(i=0;i<ntot*m;i++)
    TV2[i] = 0.0;


  build_full_Toeplitz(n[k], T+displ[k], lambda[k], T2);

/*
  print_full_Toeplitz(n[k], T2);
  fprintf(file, "---\n");
  print_full_matrix(n[k], m_rowwise, (V+dispn[k]));
  fprintf(file, "---\n");
  stmm_cblas(n[k], m_rowwise, T2, (V+dispn[k]), TV2);
  fprintf(file, "---\n");
  print_full_matrix(n[k], m_rowwise, TV2);
*/

  for(i=0;i<ntot*m;i++)
    TV2[i] = 0.0;


// Make the computation for only the first column
   cblas_dgemm (CblasColMajor, CblasNoTrans,CblasNoTrans, n[k], mblocks[k]*m_rowwise, n[k], 1, T2, n[k], &(*(V+dispn[k]+lvtot*l)), n[k], 1, TV2, n[k]);


  fprintf(file, "---\n");

  print_full_matrix(n[k], m_rowwise, TV2);


	// compute difference
//	for(i=0;i<n[k]*mblocks[k];i++)
//	  res += fabs(TV[i+dispn[k]+ntot*l]-TV2[i]);

        for(i=0;i<n[k]*mblocks[k]*m_rowwise;i++)
          res += fabs(TV[i+dispn[k]+lvtot*l]-TV2[i]);


  fprintf(file, "=========\n");
  for (j=0;j<m_rowwise;j++) {
  for(i=0;i<n[k];i++) // output matrix
    fprintf(file, "[%d] TV[%d]=%f\n", rank, i+j*n[k], TV[i+j*n[k]]);
  fprintf(file, "------\n");
  }

//        for (i=0;i<n[k]*mblocks[k];i++)
//          if(fabs(TV[i+dispn[k]+ntot*l]-TV2[i])>1e-5) {

	for (i=0;i<n[k]*mblocks[k]*m_rowwise;i++)
	  if(fabs(TV[i+dispn[k]+lvtot*l]-TV2[i])>1e-5) {
	    printf("bloc %d : TV[%d,0]=%e \t TV2[%d,0]=%e \t diff=%e\n", k, i+dispn[k], TV[i+dispn[k]], i, TV2[i], TV[i]-TV2[i]);
	    nb_false += 1;}
	res = res/(n[k]*mblocks[k]*m_rowwise); 
	
	if (fabs(res)<1e-8)
	  printf("Success ! (difference is about %e)\n", fabs(res));
	else {
	  printf("Test failed %d ... (difference is about %e)\n", k, fabs(res));
	  printf("Number of false values : %d/%d\n",nb_false, n[k]*mblocks[k]);}
      } } }
  
  if(rank==0){
    free(T2);
    free(V);
    free(TV);
    free(TV2);
   }
  free(T);
  free(Vrank);
  free(TVrank);  
  
  return 0;
}


