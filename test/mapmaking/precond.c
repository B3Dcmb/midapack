// Midapack library
// mapmaking code example using the Midapack library - release 1.2b, Nov 2012
// Compute the diagonal and block-diagonal Jacobi preconditioners for the PCG

/** @file   precond.c
    @author Frederic Dauvergne
    @date   November 2012
    @Last_update October 2018 by Hamza El Bouhargani*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include "midapack.h"

int precondblockjacobilike(Mat A, Mat *BJ, int *lhits)
{
  int           i, j, k ;                       // some indexes
  int           m, n, rank, size;
  int *indices_new;
  double *vpixBlock, *vpixBlock_loc;
  double det, invdet;
  int pointing_commflag = 1;

  m = A.m;
  n = A.lcount;
  MPI_Comm_rank(A.comm, &rank);                 //
  MPI_Comm_size(A.comm, &size);

  indices_new = (int *) malloc((A.nnz)*n * sizeof(int));
  vpixBlock_loc = (double *) malloc(n * sizeof(double));
  vpixBlock = (double *) malloc(n*(A.nnz) * sizeof(double));

  //Init vpixBlock
  for(i=0;i<n*(A.nnz);i++)
    vpixBlock[i] = 0.;

  //Compute local AtA blocks (local sum over same pix blocks)
  for(i=0;i<A.m;i++){
    for(j=0;j<A.nnz;j++){
      for(k=0;k<A.nnz;k++){
        vpixBlock[(A.nnz)*A.indices[i*(A.nnz)+j]+k] += A.values[i*(A.nnz)+j]*A.values[i*(A.nnz)+k];
      }
    }
  }

  //communicate with the other processes to have the global reduce
    for(i=0;i<n*(A.nnz);i+=(A.nnz)*(A.nnz)){
    for(j=0;j<(A.nnz);j++){
      vpixBlock_loc[(i/(A.nnz))+j] = vpixBlock[i+j];
    }
  }
  commScheme(&A, vpixBlock_loc, 2);
  for(i=0;i<n*(A.nnz);i+=(A.nnz)*(A.nnz)){
    for(j=0;j<(A.nnz);j++){
      vpixBlock[i+j] = vpixBlock_loc[(i/(A.nnz))+j];
    }
  }

  for(i=3;i<n*(A.nnz);i+=(A.nnz)*(A.nnz)){
    for(j=0;j<(A.nnz);j++)
      vpixBlock_loc[(i-3)/(A.nnz)+j] = vpixBlock[i+j];
  }
  commScheme(&A, vpixBlock_loc, 2);
  for(i=3;i<n*(A.nnz);i+=(A.nnz)*(A.nnz)){
    for(j=0;j<(A.nnz);j++)
      vpixBlock[i+j] = vpixBlock_loc[(i-3)/(A.nnz)+j];
  }

  for(i=6;i<n*(A.nnz);i+=(A.nnz)*(A.nnz)){
    for(j=0;j<(A.nnz);j++)
      vpixBlock_loc[(i-6)/(A.nnz)+j] = vpixBlock[i+j];
  }
  commScheme(&A, vpixBlock_loc, 2);
  for(i=6;i<n*(A.nnz);i+=(A.nnz)*(A.nnz)){
    for(j=0;j<(A.nnz);j++)
      vpixBlock[i+j] = vpixBlock_loc[(i-6)/(A.nnz)+j];
  }

  //Compute the inverse of the global AtA blocks (beware this part is only valid for nnz = 3)
  for(i=0;i<n*(A.nnz);i+=(A.nnz)*(A.nnz)){
    lhits[(int)i/9] = (int)vpixBlock[i];
    //init 3x3 block
    double block[3][3];
    for(j=0;j<3;j++){
      for(k=0;k<3;k++){
        block[j][k] = vpixBlock[i+(j*3)+k];
      }
    }
    //Compute det
    det = block[0][0] * (block[1][1] * block[2][2] - block[2][1] * block[1][2]) -
             block[0][1] * (block[1][0] * block[2][2] - block[1][2] * block[2][0]) +
             block[0][2] * (block[1][0] * block[2][1] - block[1][1] * block[2][0]);

    if(fabs(det) > 1e-7){
    invdet = 1 / det;

    //Compute the inverse coeffs
    vpixBlock[i] = (block[1][1] * block[2][2] - block[2][1] * block[1][2]) * invdet;
    vpixBlock[i+1] = (block[0][2] * block[2][1] - block[0][1] * block[2][2]) * invdet;
    vpixBlock[i+2] = (block[0][1] * block[1][2] - block[0][2] * block[1][1]) * invdet;
    vpixBlock[i+3] = (block[1][2] * block[2][0] - block[1][0] * block[2][2]) * invdet;
    vpixBlock[i+4] = (block[0][0] * block[2][2] - block[0][2] * block[2][0]) * invdet;
    vpixBlock[i+5] = (block[1][0] * block[0][2] - block[0][0] * block[1][2]) * invdet;
    vpixBlock[i+6] = (block[1][0] * block[2][1] - block[2][0] * block[1][1]) * invdet;
    vpixBlock[i+7] = (block[2][0] * block[0][1] - block[0][0] * block[2][1]) * invdet;
    vpixBlock[i+8] = (block[0][0] * block[1][1] - block[1][0] * block[0][1]) * invdet;
  }
  else{// Cancel the resolution of pixels corresponding to singular blocks
    vpixBlock[i] = 0.;
    vpixBlock[i+4] = 0.;
    vpixBlock[i+8] = 0.;
    vpixBlock[i+1] = 0.;
    vpixBlock[i+2] = 0.;
    vpixBlock[i+3] = 0.;
    vpixBlock[i+5] = 0.;
    vpixBlock[i+6] = 0.;
    vpixBlock[i+7] = 0.;
  }
  }

  //Define Block-Jacobi preconditioner indices
  for(i=0;i<n;i++){
    for(j=0;j<(A.nnz);j++){
        indices_new[i*(A.nnz)+j] = A.lindices[(A.nnz)*((int)i/(A.nnz))+j];
    }
  }
  // Init Block-Jacobi preconditioner
  MatInit(BJ, n, A.nnz, indices_new, vpixBlock, pointing_commflag, MPI_COMM_WORLD);


  return 0;

}
int precondjacobilike_avg(Mat A, Tpltz Nm1, double *c)
{

  int           i, j, k ;                       // some indexes
  int           m, n, rank, size;
  double        localreduce;                    //reduce buffer
  double        st, t;                          //timers

  m=A.m;                                        //number of local time samples
  n=A.lcount;                                   //number of local pixels
  MPI_Comm_rank(A.comm, &rank);                 //
  MPI_Comm_size(A.comm, &size);                 //

  double diagNm1;


//Compute diag( AtA )
  DiagAtA(&A, c, 2);

//multiply by the diagonal Toeplitz
  diagNm1 = Nm1.tpltzblocks[0].T_block[0];

  printf("diagNm1 = %f \n", i, diagNm1 );
  for(j=0; j<n; j++)
    c[j] = diagNm1 * c[j];

// compute c inverse vector
  for(j=0; j<n; j++)
    c[j] = 1./c[j] ;


  return 0;
}




int precondjacobilike(Mat A, Tpltz Nm1, int *lhits, double *cond, double *vpixDiag)
{

  int           i, j, k ;                       // some indexes
  int           m, n, rank, size;
  double        localreduce;                    //reduce buffer
  double        st, t;                          //timers

  MPI_Comm_rank(A.comm, &rank);                 //
  MPI_Comm_size(A.comm, &size);                 //

  m=A.m;                                        //number of local time samples
  n=A.lcount;                                   //number of local pixels

//Compute local diag( At diag(N^1) A )
  getlocDiagN(&A, Nm1, vpixDiag);

/*
  for(i=0; i<10; i++)
    printf("rank=%d, vpixDiag[%d]=%lf \n", rank, i, vpixDiag[i]) ;

    printf("rank=%d, vpixDiag[n-1]=%lf \n", rank, vpixDiag[n-1]);
*/

//communicate with the other processes to have the global reduce
  commScheme(&A, vpixDiag, 2);
  for(i=0;i<50;i++){
    printf("global AtA block: vpixDiag[%d]=%f\n",i,vpixDiag[i]);
  }
// compute the inverse vector
  for(i=0; i<n; i++){
    if(i%3 == 0){
      lhits[(int)i/3] = (int)vpixDiag[i];
      cond[(int)i/3] = vpixDiag[i+1] + vpixDiag[i+2];
    }
    vpixDiag[i] = 1./vpixDiag[i] ;
  }
  return 0;
}



//do the local diag( At diag(Nm1) A ) with as output a vector in the pixel domain
int getlocDiagN(Mat *A, Tpltz Nm1, double *vpixDiag)
{
  int           i, j, k ;                       // some indexes
  int           m;

  m=Nm1.local_V_size;  //number of local time samples
//  int nnz=(A->nnz);

  //Define the indices for each process
  int idv0, idvn;  //indice of the first and the last block of V for each processes
  int *nnew;
  nnew = (int*) calloc(Nm1.nb_blocks_loc, sizeof(int));
  int64_t idpnew;
  int local_V_size_new;
//get idv0 and idvn
  get_overlapping_blocks_params( Nm1.nb_blocks_loc, Nm1.tpltzblocks, Nm1.local_V_size, Nm1.nrow, Nm1.idp, &idpnew, &local_V_size_new, nnew, &idv0, &idvn);
 // double *vpixDiag;
 // vpixDiag = (double *) malloc(A->lcount *sizeof(double));

  int istart, il, istartn;
  for(i=0; i < A->lcount; i++)
    vpixDiag[i]=0.0;//0.0;

  int vShft=idpnew-Nm1.idp; //=Nm1.tpltzblocks[idv0].idv-Nm1.idp in principle
/*
  printf("Nm1.idp=%d, idpnew=%d, vShft=%d\n", Nm1.idp, idpnew, vShft);
  printf("idv0=%d, idvn=%d\n", idv0, idvn);
  printf("Nm1.nb_blocks_loc=%d, Nm1.local_V_size=%d\n", Nm1.nb_blocks_loc, Nm1.local_V_size);

  for(i=0; i < Nm1.nb_blocks_loc; i++)
    printf("Nm1.tpltzblocks[%d].idv=%d\n", i, Nm1.tpltzblocks[i].idv);
*/

//go until the first piecewise stationary period
    for (i= 0; i<vShft; i++) {
      for (j=0; j<(A->nnz); j++)
        vpixDiag[A->indices[i*(A->nnz)+j]]+=(A->values[i*(A->nnz)+j]*A->values[i*(A->nnz)+j]);
    }

//temporary buffer for one diag value of Nm1
  int diagNm1;
//loop on the blocks
  for(k=idv0; k<(idv0+Nm1.nb_blocks_loc); k++) {
  if (nnew[idv0]>0) {  //if nnew==0, this is a wrong defined block

    if (k+1<idv0+Nm1.nb_blocks_loc)   //if there is a next block, compute his next first indice
      istartn= Nm1.tpltzblocks[k+1].idv-Nm1.idp ;
    else
      istartn = 0;

    istart = max( 0, Nm1.tpltzblocks[k].idv-Nm1.idp);

    //if block cut from the left:
    if (k==idv0)
      il     = min( Nm1.tpltzblocks[k].n, Nm1.tpltzblocks[k].idv + Nm1.tpltzblocks[k].n - Nm1.idp );
    //if block cut from the right:
    if (k==idv0+Nm1.nb_blocks_loc-1)
      il     = min(il , (Nm1.idp + Nm1.local_V_size) - Nm1.tpltzblocks[k].idv );
    //if block alone in the middle, and cut from both sides
    if (Nm1.nb_blocks_loc==1)
      il     = min(il , Nm1.local_V_size);

    //get the diagonal value of the Toeplitz
    diagNm1 = Nm1.tpltzblocks[k].T_block[0];

/*
    printf("istart=%d, il=%d, istartn=%d\n", istart, il, istartn);
    printf("Nm1.tpltzblocks[k].idv=%d, Nm1.tpltzblocks[k].n=%d, Nm1.idp=%d\n", Nm1.tpltzblocks[k].idv, Nm1.tpltzblocks[k].n, Nm1.idp);
*/
//a piecewise stationary period
    for (i= istart; i<istart+il; i++) {
      for (j=0; j<(A->nnz); j++)
        vpixDiag[A->indices[i*(A->nnz)+j]]+=(A->values[i*(A->nnz)+j]*A->values[i*(A->nnz)+j])*diagNm1;
    }

//continue until the next period if exist
    for (i= istart+il; i<istartn; i++) {
      for (j=0; j<(A->nnz); j++)
        vpixDiag[A->indices[i*(A->nnz)+j]]+=(A->values[i*(A->nnz)+j]*A->values[i*(A->nnz)+j]);
    }


  }}//end of the loop over the blocks


  return 0;
}


//communication scheme in the pixel domain for the vector vpixDiag
//extract from a Madmap routine
int commScheme(Mat *A, double *vpixDiag, int pflag){
  int i, j, k;
  int nSmax, nRmax;
  double *lvalues;

#if W_MPI
  lvalues = (double *) malloc(A->lcount *sizeof(double));    /*<allocate and set to 0.0 local values*/
  memcpy(lvalues, vpixDiag, (A->lcount) *sizeof(double)); /*<copy local values into result values*/

  nRmax=0;
  nSmax=0;

  if(A->flag  == BUTTERFLY){                                  /*<branch butterfly*/
    //memcpy(out_values, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    for(k=0; k< A->steps; k++)                                  /*compute max communication buffer size*/
      if(A->nR[k] > nRmax)
        nRmax = A->nR[k];
    for(k=0; k< A->steps; k++)
      if(A->nS[k] > nSmax)
        nSmax = A->nS[k];

    double *com_val;
    com_val=(double *) malloc( A->com_count *sizeof(double));
    for(i=0; i < A->com_count; i++){
      com_val[i]=0.0;
    }
//already done    memcpy(vpixDiag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    m2m(lvalues, A->lindices, A->lcount, com_val, A->com_indices, A->com_count);
    butterfly_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
    m2m(com_val, A->com_indices, A->com_count, vpixDiag, A->lindices, A->lcount);
    free(com_val);
  }
  else if(A->flag == RING){
//already done    memcpy(vpixDiag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    for(k=1; k< A->steps; k++)                                  /*compute max communication buffer size*/
      if(A->nR[k] > nRmax)
        nRmax = A->nR[k];

    nSmax = nRmax;
    ring_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, vpixDiag, A->steps, A->comm);
  }
  else if(A->flag == NONBLOCKING){
//already done    memcpy(vpixDiag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    ring_nonblocking_reduce(A->R, A->nR, A->S, A->nS, lvalues, vpixDiag, A->steps, A->comm);
  }
  else if(A->flag == NOEMPTY){
//already done    memcpy(vpixDiag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    int ne=0;
    for(k=1; k< A->steps; k++)
      if(A->nR[k]!=0)
        ne++;

    ring_noempty_reduce(A->R, A->nR, ne, A->S, A->nS, ne, lvalues, vpixDiag, A->steps, A->comm);
  }
  else{
    return 1;
  }
#endif
  free(lvalues);
  return 0;
}







/** @brief Compute Diag(A' diag(Nm1) A).
    @param out_values local output array of doubles*/

int DiagAtA(Mat *A, double *diag, int pflag){
  int i, j, k;
  int nSmax, nRmax;
  double *lvalues;

  lvalues = (double *) malloc(A->lcount *sizeof(double));    /*<allocate and set to 0.0 local va
lues*/
  for(i=0; i < A->lcount; i++)
    lvalues[i]=0.0;

//Naive computation with a full defined diag(Nm1):
  for(i=0; i< A->m; i++)
    for (j=0; j< A->nnz; j++)                                   /*<dot products */
      lvalues[A->indices[i*(A->nnz)+j]]+=(A->values[i*(A->nnz)+j]*A->values[i*(A->nnz)+j]) ;//*vdiagNm1[i];



#if W_MPI
  nRmax=0;  nSmax=0;

  if(A->flag  == BUTTERFLY){                                  /*<branch butterfly*/
    //memcpy(out_values, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    for(k=0; k< A->steps; k++)                                  /*compute max communication buffer size*/
      if(A->nR[k] > nRmax)
        nRmax = A->nR[k];
    for(k=0; k< A->steps; k++)
      if(A->nS[k] > nSmax)
        nSmax = A->nS[k];

    double *com_val;
    com_val=(double *) malloc( A->com_count *sizeof(double));
    for(i=0; i < A->com_count; i++){
      com_val[i]=0.0;
    }
    memcpy(diag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    m2m(lvalues, A->lindices, A->lcount, com_val, A->com_indices, A->com_count);
    butterfly_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
    m2m(com_val, A->com_indices, A->com_count, diag, A->lindices, A->lcount);
    free(com_val);
  }
  else if(A->flag == RING){
    memcpy(diag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    for(k=1; k< A->steps; k++)                                  /*compute max communication buffer size*/
      if(A->nR[k] > nRmax)
        nRmax = A->nR[k];

    nSmax = nRmax;
    ring_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, diag, A->steps, A->comm);
  }
  else if(A->flag == NONBLOCKING){
    memcpy(diag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    ring_nonblocking_reduce(A->R, A->nR, A->S, A->nS, lvalues, diag, A->steps, A->comm);
  }
  else if(A->flag == NOEMPTY){
    memcpy(diag, lvalues, (A->lcount) *sizeof(double)); /*<copy local values into result values*/
    int ne=0;
    for(k=1; k< A->steps; k++)
      if(A->nR[k]!=0)
        ne++;

    ring_noempty_reduce(A->R, A->nR, ne, A->S, A->nS, ne, lvalues, diag, A->steps, A->comm);
  }
  else{
    return 1;
  }
#endif
  free(lvalues);
  return 0;
}




int get_pixshare_pond(Mat A, double *pixpond )
{

  int           i, j, k ;                       // some indexes
  int           m, n, rank, size;
  double        localreduce;                    //reduce buffer
  double        st, t;                          //timers

//  double        *eyesdble;

  MPI_Comm_rank(A.comm, &rank);                 //
  MPI_Comm_size(A.comm, &size);                 //

  m=A.m;                                        //number of local time samples
  n=A.lcount;                                   //number of local pixels

//  eyesdble = (double *) malloc(n*sizeof(double));

//create a eyes local vector
  for(i=0; i<n; i++)
    pixpond[i] = 1.;

//communicate with the others processes to have the global reduce
  commScheme(&A, pixpond, 0);

// compute the inverse vector
  for(i=0; i<n; i++)
    pixpond[i] = 1./pixpond[i] ;


  return 0;
}
