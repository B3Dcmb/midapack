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

    nb_blocs = 2;
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
//    lambda[2] = 3;
//    lambda[3] = 3;

    int nsample=40;
//    for(j=0;j<nb_blocs;j++)
//      nsample += (n[j]*mblocks[j]);

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


  ntot=40;
 
  T  = (double *) calloc(lambdatot,sizeof(double));
  if(rank==0) {
    T2 = (double *) calloc(ntot*ntot,sizeof(double));
    V  = (double *) calloc(mbloc*ntot*m   ,sizeof(double));  
    TV = (double *) calloc(mbloc*ntot*m   ,sizeof(double));  
    TV2= (double *) calloc(mbloc*ntot*m   ,sizeof(double));  
    
    
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

    

  printf("ntot=%d, mbloc=%d\n", ntot, mbloc);
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


/*
    for(i=0;i<ntot*m*mbloc;i++) // input matrix
      printf("V[i]=%f\n", V[i]);
*/

   }//End Rank0



   int *lv;
   lv = (int *) calloc(nb_blocs,sizeof(int));

    for(i=0;i<nb_blocs;i++) // input matrix
      lv[i] = mblocks[i]*n[i];

// if (rank==0) {
//    for(i=0;i<nb_blocs;i++) {
//      printf("mblocks[i=%d]=%d\n", i, mblocks[i]);
//      printf("n[i=%d]=%d\n", i, n[i]);
//      printf("lv[i=%d]=%d\n", i, lv[i]);
//    }
//}


  MPI_Bcast(T, lambdatot, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    
  int maxsize = nranks[0];
  for (i=1; i<size; i++)
    if (nranks[i]>maxsize)
      maxsize=nranks[i];
  
  // MPI test
  double *Vrank = (double *) calloc(maxsize, sizeof(double));
  double *TVrank = (double *) calloc(maxsize,sizeof(double));
  
  int nrow = ntot;
  int id0 = displs[rank];
  int local_V_size=nranks[rank];

  
  MPI_Scatterv(V, nranks, displs, MPI_DOUBLE, Vrank, maxsize, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
/*
    for(i=0;i<nranks[rank];i++) // input matrix
      fprintf(file, "[%d] Vrank[%d]=%f\n", rank, i, Vrank[i]);

      fprintf(file, "=========\n");

    for(i=0;i<nb_blocs;i++) 
      fprintf(file, "[%d] n[%d]=%d\n", rank, i, n[i]);
    fprintf(file, "[%d] m=%d, nrow=%d, nb_blocs=%d\n", rank, m, nrow, nb_blocs);
    fprintf(file, "[%d] id0=%d, local_V_size=%d\n", rank, id0, local_V_size);

    for(i=0;i<nb_blocs;i++)
      fprintf(file, "[%d] dispn[%d]=%d\n", rank, i, dispn[i]);
*/



  int ngap;
  int *id0gap;
  int *lgap;

  ngap=2;
  id0gap = (int *) calloc(ngap, sizeof(int));
  lgap = (int *) calloc(ngap, sizeof(int));

/*
  id0gap[0]=3;
  lgap[0]=5;
  id0gap[1]=10;
  lgap[1]=2;
*/
/*
  id0gap[0]=10;
  lgap[0]=5;
*/
/*
  id0gap[0]=1;
  lgap[0]=3;
  id0gap[1]=5;
  lgap[1]=3;
  id0gap[2]=8;
  lgap[2]=14;
*/

/*
  id0gap[0]=3;
  lgap[0]=4;
  id0gap[1]=9;
  lgap[1]=14;
*/

  id0gap[0]=1;
  lgap[0]=4;
  id0gap[1]=6;
  lgap[1]=3;



//  for(i=0;i<nb_blocs;i++)
//    fprintf(file, "[%d] idv[%d]=%d\n", rank, i, idv[i]);
//  for(i=0;i<nb_blocs;i++)
//    fprintf(file, "[%d] n[%d]=%d ; nnew[%d]=%d\n", rank, i, n[i], i, nnew[i]);

  int nb_blocs_local=nb_blocs;
  int nb_blocs_all=nb_blocs;
  int *idv = (int *) calloc(nb_blocs, sizeof(int));

//  for(i=0;i<nb_blocs;i++)
//    idv[i]=dispn[i];

  idv[0]=0;
  idv[1]=20;

  fprintf(file, "=========\n");
  fprintf(file, "[%d] id0=%d, local_V_size=%d\n", rank, id0, local_V_size);
  fprintf(file, "[%d] nrow=%d ; m=%d\n", rank, nrow, m);
  fprintf(file, "[%d] nb_blocs_local=%d ; nb_blocs_all=%d\n", rank, nb_blocs_local, nb_blocs_all);
  for(i=0;i<nb_blocs;i++)
    fprintf(file, "[%d] n[%d]=%d \n", rank, i, n[i]);
  for(i=0;i<nb_blocs;i++)
    fprintf(file, "[%d] idv[%d]=%d\n", rank, i, idv[i]);
  for(i=0;i<nb_blocs;i++)
    fprintf(file, "[%d] lambda[%d]=%d\n", rank, i, lambda[i]);
  fprintf(file, "---------\n");
    fprintf(file, "[%d] ngap=%d\n", rank, ngap);
  for(i=0;i<ngap;i++)
    fprintf(file, "[%d] id0gap[%d]=%d ; lgap[%d]=%d\n", rank, i, id0gap[i],  i, lgap[i]);



  int flag_param_small=0;
  int param_small_fixed=0;

  int nb_blocs_gappy;
  double *Tgappy;
  int *ngappy;
  int *lambdagappy;
  int *idvgappy;


//  get_gappy_blocks_params(n, m, nrow, T, nb_blocs_local, nb_blocs_all, lambda, idv, id0gap, lgap, ngap, &nb_blocs_gappy, Tgappy, idvgappy, ngappy, lambdagappy, flag_param_small, param_small_fixed);

  int idp=id0;


  if(rank==0){
 for (i=0; i<nrow; i++)
   printf("V=%f\n", V[i]);
  }


 
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




