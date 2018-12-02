#include "toeplitz.h"
extern int NFFT;
extern int FFTW_FLAG;
extern int PRINT_RANK;

//v1.1
//this is the updated algorithm for mpi_stbmm with the m_rowwise - fd@apc
//this parameter correspond to the number of column for the local rank rowwise order
//
//need to be done:
//- change some variables names
//- use 64bits integer
//nrow, idp, idv (this one maybe use local rank coordinates)
//idvm0, idvmn...to see 
//changer offset0 offsetn toSendLeft Right avec lambda-1
//delete v1_size
//offset -> OffsetLeft - Right




//=========================================================================

/// Performs the multiplication of a symmetric, Toeplitz block-diagonal matrix, T, by an arbitrary matrix, V, distributed over processes in the generalized column-wise way. 
/** @ingroup group12
    Each process performs the multiplication sequentially for each diagonal block and based on
    the sliding window algorithm. Prior to that MPI calls are used to exchange data between 
    neighboring process. Each of the diagonal blocks is a symmetric, band-diagonal Toeplitz
    matrix, which can be different for each block.
    The parameters are :
    \param V \b [input] distributed data matrix (with the convention V(i,j)=V[i+j*n]) ;
             \b [out] result of the product TV
    \param n number of rows for each Toeplitz block as stored in T
    \param m number of columns of the global data matrix V
    \param nrow number of rows of the global data matrix V
    \param T Toeplitz matrix composed of the non-zero entries of the first row of each Toeplitz
    block and concatenated together have to be arranged in the increasing order of n without
    repetitions and overlaps.
    \param nb_blocks_all number of all Toeplitz block on the diagonal of the full Toeplitz
    matrix
    \param nb_blocks_local number of Toeplitz blocks as stored in T 
    \param lambda half bandwith size for each Toeplitz block stroed in T
    \param idv global row index defining for each Toeplitz block as stored in the vector T 
    first element of the interval to which given Toeplitz matrix is to be applied.
    \param idp global index of the first element of the local part of V
    \param local_V_size a number of all elements in local V
    \param comm MPI communicator 
*/
int mpi_stbmm(double **V, int *n, int m, int nrow, int m_rowwise, double *T, int nb_blocks_local, int nb_blocks_all, int *lambda, int *idv, int idp, int local_V_size, MPI_Comm comm)
{

  //MPI parameters 
  int rank;   //process rank 
  int size;   //process number 

  MPI_Status status;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  PRINT_RANK=rank;

  FILE *file;
  file = stdout;

/*
  FILE *file;
  file = stdout;
  PRINT_RANK=rank ;

  char filename [1024];
  sprintf(filename,"output-mpi_stbmm_rowwise-%d.txt", rank);
  file = fopen(filename, "w");
*/

  //identification of the mpi neighbours process to communicate when there is a shared block
  int right = rank+1;
  int left  = rank-1;

  int i,j,k;  //some indexes 


  //Define the indices for each process
  int idv0, idvn;  //indice of the first and the last block of V for each processes

  int *nnew;
  nnew = (int*) calloc(nb_blocks_local, sizeof(int));
  int idpnew;
  int local_V_size_new;  
  int n_rowwise=local_V_size;


  int status_params = get_overlapping_blocks_params( nb_blocks_local, idv, n, local_V_size, nrow, idp, &idpnew, &local_V_size_new, nnew, &idv0, &idvn);

  if (VERBOSE==1)
    printf("status_params=%d\n", status_params);

  if( status_params == 0) {
  free(nnew);
    return(0); // no work to be done
  }

  if (lambda[idv0]==0 || lambda[idvn]==0)
    return print_error_message (2, __FILE__, __LINE__);


if (VERBOSE==1) {  //print on screen news parameters definition if VERBOSE
    fprintf(file, "new parameters caracteristics:\n");
    fprintf(file, "[%d] idp=%d ; idpnew=%d\n", rank, idp, idpnew);
    fprintf(file, "[%d] local_V_size=%d ; local_V_size_new=%d\n", rank, local_V_size, local_V_size_new);
    for(i=0;i<nb_blocks_local;i++)
      fprintf(file, "[%d] n[%d]=%d ; nnew[%d]=%d\n", rank, i, n[i], i, nnew[i]);
}


  int vShft=idpnew-idp;   //new first element of relevance in V

  //Define the column indices:
  //index of the first and the last column of V for the current process
  int idvm0 = idpnew/nrow;
  int idvmn = (idpnew+local_V_size_new-1)/nrow;
  //number of columns of V for the current process
  int ncol_rank = idvmn-idvm0+1;
  //number of blocks for the current process with possibly repetitions 
  int nb_blocks_rank;

  if(ncol_rank == 1) // Empty process not allowed 
    nb_blocks_rank = idvn - idv0 + 1;
  else
    nb_blocks_rank = (ncol_rank-2)*nb_blocks_local + (nb_blocks_local-idv0) + (idvn+1);  //in this case nb_blocks_local = nblocs_all

  if (VERBOSE==1) {  //print on screen 
  fprintf(file, "[%d] nb_blocks_rank=%d, nb_blocks_local=%d\n", rank, nb_blocks_rank, nb_blocks_local);
  }

  //Define the indices for the first and the last element in each blocks
  int idvp0 = idpnew%nrow-idv[idv0];  //index of the first element of the process in the first block
  int idvpn;  //reverse index of the last element of the process in the last block
              //It's the number of remaining elements needed to fully complete the last block
  idvpn = idv[idvn]+nnew[idvn]-1 - (idpnew+local_V_size_new-1)%nrow ;


  //Define the offsets for the first and last blocks of the process for V1
 // int v1_size = local_V_size_new;  //length of the new vector V1
  int offset0, offsetn;
  int distcorrmin_idv0 = lambda[idv0]-1;
  int distcorrmin_idvn = lambda[idvn]-1;

  //if(idvp0 != 0) 
    offset0 = min( idvp0, distcorrmin_idv0); 
  //if(idvpn != 0) 
    offsetn = min(idvpn, distcorrmin_idvn);  


  int toSendLeft=0;
  if(offset0!=0) {
    toSendLeft = min( idv[idv0]+nnew[idv0]-idpnew%nrow, distcorrmin_idv0); 
  }
  int toSendRight=0;
  if( offsetn != 0) {
    toSendRight = min( (idpnew+local_V_size_new)%nrow-idv[idvn], distcorrmin_idvn); 
  }

 int flag_optimlambda=1; //to allocate only the memory place needed

 int lambdaOut_offset;
 int lambdaIn_offset;
 double *LambdaOut;
 int lambdaOut_size, lambdaIn_size;

 if (flag_optimlambda==1) {
  LambdaOut=(double *) calloc((toSendLeft+toSendRight)*m_rowwise, sizeof(double));
  lambdaOut_offset = toSendLeft*m_rowwise;
  lambdaIn_offset = offset0*m_rowwise;
  lambdaOut_size = (toSendLeft+toSendRight)*m_rowwise ; 
  lambdaIn_size = (offset0+offsetn)*m_rowwise;
 }
 else {
  LambdaOut=(double *) calloc((lambda[idv0]+lambda[idvn])*m_rowwise, sizeof(double));
  lambdaOut_offset = lambda[idv0]*m_rowwise;
  lambdaIn_offset = lambda[idv0]*m_rowwise;
  lambdaOut_size = (lambda[idv0]+lambda[idvn])*m_rowwise;
  lambdaIn_size = (lambda[idv0]+lambda[idvn])*m_rowwise;
 }


  if(offset0!=0) {
    for (j=0;j<m_rowwise;j++)
    for (i=0;i<toSendLeft;i++)
      LambdaOut[i+j*toSendLeft]=(*V)[i+j*n_rowwise]; //good because toSendLeft=0 if it 
  }                                                   //doesnt start on a the first block. 
  if( offsetn != 0) {
    for (j=0;j<m_rowwise;j++)
    for (i=0;i<toSendRight;i++)
      LambdaOut[i+j*toSendRight+lambdaOut_offset]=(*V)[i+j*n_rowwise+local_V_size-toSendRight];
  }                                        //good too using same argument than for offset0!=0
                                           //if local_V_size!=local_V_size_new+vShft mean there is extra
                                           //terms a the end and so offsetn=0
                                           //idpnew+local_V_size_new = idp+local_V_size and vShft=idpnew-idp
                                           //so local_V_size=vShft+local_V_size_new
  if(rank==0 || offset0==0)
    left = MPI_PROC_NULL;
  if(rank==size-1 || offsetn==0)
    right = MPI_PROC_NULL;

  double *LambdaIn=(double *) calloc(lambdaIn_size, sizeof(double));

/*
//VERBOSE_DEV
  fprintf(file, "[%d] toSendLeft=%d ; toSendRight=%d\n", rank, toSendLeft, toSendRight);
  fprintf(file, "[%d] offset0=%d ; offsetn=%d\n", rank, offset0, offsetn);
  fprintf(file, "[%d] lambdaOut_size=%d ; lambdaIn_size=%d\n", rank, lambdaOut_size, lambdaIn_size);
  fprintf(file, "[%d] lambdaOut_offset=%d ; lambdaIn_offset=%d\n", rank, lambdaOut_offset, lambdaIn_offset);

  fprintf(file, "communications Out:\n");
    for (i=0;i<lambdaOut_size;i++)
      fprintf(file, "[%d] LambdaOut[%d]=%f\n", rank, i, LambdaOut[i]);

  fprintf(file, "communications In:\n");
    for (i=0;i<lambdaIn_size;i++)
      fprintf(file, "[%d] LambdaIn[%d]=%f\n", rank, i, LambdaIn[i]);
//END VERBOSE DEV
*/

 int flag_blockingcomm=0;  //to use blocking comm
 MPI_Request requestLeft_r, requestLeft_s;
 MPI_Request requestRight_r, requestRight_s;

  if (flag_blockingcomm==1) {
//send and receive data
  MPI_Sendrecv( LambdaOut, toSendLeft*m_rowwise, MPI_DOUBLE, left,  MPI_USER_TAG, (LambdaIn+lambdaIn_offset), offsetn*m_rowwise, MPI_DOUBLE, right, MPI_USER_TAG, comm, &status);
  MPI_Sendrecv( (LambdaOut+lambdaOut_offset), toSendRight*m_rowwise, MPI_DOUBLE, right,  MPI_USER_TAG, LambdaIn, offset0*m_rowwise, MPI_DOUBLE, left, MPI_USER_TAG, comm, &status);

  }
  else {
//to the Left
//! [communication Call example]
  MPI_Irecv((LambdaIn+lambdaIn_offset), offsetn*m_rowwise, MPI_DOUBLE, right, MPI_USER_TAG, comm, &requestLeft_r);
  MPI_Isend(LambdaOut, toSendLeft*m_rowwise, MPI_DOUBLE, left, MPI_USER_TAG, comm, &requestLeft_s);
//! [communication Call example]
//to the Right
  MPI_Irecv(LambdaIn, offset0*m_rowwise, MPI_DOUBLE, left, MPI_USER_TAG, comm, &requestRight_r);
  MPI_Isend((LambdaOut+lambdaOut_offset), toSendRight*m_rowwise, MPI_DOUBLE, right, MPI_USER_TAG, comm, &requestRight_s);

 }


//size of the first and the last block for the current process
  int v0rank_size, vnrank_size;
  if (nb_blocks_rank == 1) {  //only one block 
    v0rank_size = ((idpnew+local_V_size_new-1)%nrow +1) - idpnew%nrow + offset0 + offsetn;
    vnrank_size = 0; //just for convenience - no really need it
  }
  else { //more than one block
    v0rank_size = idv[idv0] + nnew[idv0] - idpnew%nrow + offset0;
    vnrank_size = ((idpnew+local_V_size_new-1)%nrow +1) - idv[idvn] + offsetn;
  }


if (flag_blockingcomm!=1) {
  //MPI_Wait for lambda comm
  MPI_Wait(&requestLeft_r, &status);
  MPI_Wait(&requestLeft_s, &status);
  MPI_Wait(&requestRight_r, &status);
  MPI_Wait(&requestRight_s, &status);

}


  free(LambdaOut);


//---------------------------------------
//initialization for the blocks loop

  int idv1=0;     //old index of *V1
  int idv2=0;     //index 


  int mid;  //local number of column for the current block
  //index of the first element of the process inside the first block
  int offset_id0;
  offset_id0 = idvp0;

//fftw variables
  fftw_complex *V_fft, *T_fft;
  double *V_rfft;
  fftw_plan plan_f, plan_b;
//init local block vector
  double *V1block;
  int lambdaShft;


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//loop on the blocks inside the process
  int nfft, blocksize;
  int iblock;  //index for the loop on the blocks
//  int loopindex;
  int id; //indice of the current block

  int vblock_size;
  int id0block;

  int jj;


  for(iblock=idv0;iblock<idv0+nb_blocks_rank;iblock++) {
    id = iblock%nb_blocks_local;  //index of current block


  if(nnew[id]>0) { //the block is ok


  for( lambdaShft=k=0; k<id; k++)
    lambdaShft += lambda[k];


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//first case : First block of the process
  if(iblock==idv0) {
  if(VERBOSE)
    fprintf(file, "[%d] First block...\n", rank);

  vblock_size=v0rank_size;
  id0block=(offset_id0-offset0);

  V1block = (double *) calloc(vblock_size*m_rowwise, sizeof(double));

  for (j=0;j<m_rowwise;j++) { 
#pragma omp parallel for                        
  for (i=0;i<offset0;i++)
    V1block[i+j*vblock_size] = LambdaIn[i+j*offset0];
  }
//note: check if copyblock could be used instead.


//if (nb_blocks_rank == 1) currentsize_middlepart=vblock_size-offset0-offsetn = local_V_size_new
//else currentsize_middlepart=vblock_size-offset0
  int currentsize_middlepart=min(vblock_size-offset0, local_V_size_new);

  for (j=0;j<m_rowwise;j++) {
#pragma omp parallel for   
  for (i=0;i<currentsize_middlepart;i++)
    V1block[offset0+i+j*vblock_size] = (*V)[i+vShft+j*n_rowwise];
  }

if (nb_blocks_rank == 1) {
  for (j=0;j<m_rowwise;j++) {
#pragma omp parallel for   
  for (i=0;i<offsetn;i++) {
    V1block[vblock_size-offsetn+i+j*vblock_size] = LambdaIn[i+lambdaIn_offset+j*offsetn];
  }}
}


  //init Toeplitz arrays
  tpltz_init(vblock_size, lambda[id], &nfft, &blocksize, &T_fft, (T+lambdaShft), &V_fft, &V_rfft, &plan_f, &plan_b);
  //Toeplitz computation
  if(VERBOSE)
    fprintf(file, "[%d] Before stmm call : nfft = %d, blocksize = %d\n", rank, nfft, blocksize);
  stmm(&V1block, vblock_size, m_rowwise, 0, m_rowwise*vblock_size, T_fft, lambda[id], V_fft, V_rfft, plan_f, plan_b, blocksize, nfft);

  tpltz_cleanup(&T_fft, &V_fft, &V_rfft, &plan_f, &plan_b);


  int currentsize=min(vblock_size-offset0, local_V_size_new);
  for (j=0;j<m_rowwise;j++) {
#pragma omp parallel for
  for (i=0;i<currentsize;i++)
    (*V)[vShft+i+j*n_rowwise] = V1block[offset0+i+j*vblock_size];
  }

  free(V1block);

  }//end (First case)


//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//Generic case : Generic block of the process
  else if(iblock!=idv0 && iblock!=idv0+nb_blocks_rank-1) {
  if(VERBOSE)
    fprintf(file, "[%d] generic block...\n");

  vblock_size=nnew[id];
  id0block=0;

  V1block = (double *) calloc(vblock_size*m_rowwise, sizeof(double));

  idv1 = idv[id]-idp%nrow - vShft + offset0 +nrow*( (iblock/nb_blocks_local) );  //no need
//  idv2 = idv[id]-idp%nrow + nrow*( (iblock/nb_blocks_local) );
  idv2 = idv[id]-(idpnew)%nrow+vShft + nrow*( (iblock/nb_blocks_local) );

  for (j=0;j<m_rowwise;j++) {
#pragma omp parallel for                        
  for (i=0;i<vblock_size;i++)
    V1block[i+j*vblock_size] = (*V)[i+idv2+j*n_rowwise];
//    V1block[i] = (*V)[i+idv1-offset0+vShft];
  }

  //init Toeplitz arrays
  tpltz_init(nnew[id], lambda[id], &nfft, &blocksize, &T_fft, (T+lambdaShft), &V_fft, &V_rfft, &plan_f, &plan_b);
  //Toeplitz computation
  if(VERBOSE)
    fprintf(file, "[%d] Before stmm call : nfft = %d, blocksize = %d\n", rank, nfft, blocksize);
  stmm(&V1block, vblock_size, m_rowwise, 0, m_rowwise*vblock_size, T_fft, lambda[id], V_fft, V_rfft, plan_f, plan_b, blocksize, nfft);


  tpltz_cleanup(&T_fft, &V_fft, &V_rfft, &plan_f, &plan_b);


  for (j=0;j<m_rowwise;j++) {
#pragma omp parallel for 
  for (i=0;i<vblock_size;i++) {
    (*V)[i+idv2+j*n_rowwise] = V1block[i+j*vblock_size];
  }}


  free(V1block);

  }  //end (Generic case)

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Last case : Last block of the process
  else if(iblock==idv0+nb_blocks_rank-1 && iblock!= idv0) {
  if(VERBOSE)  
    fprintf(file, "[%d] last block...\n");

  vblock_size=vnrank_size;
  id0block=0;

  V1block = (double *) calloc(vblock_size*m_rowwise, sizeof(double));

  idv1 = idv[id] - idp%nrow - vShft + offset0  + nrow*( (iblock/nb_blocks_local) );
  //int idv2_o = idv[id]+vShft-(idpnew)%nrow + nrow*( (iblock/nb_blocks_local) );
  idv2 = idv[id]-(idpnew)%nrow+vShft + nrow*( (iblock/nb_blocks_local) );


  for (j=0;j<m_rowwise;j++) {
#pragma omp parallel for   
  for (i=0;i<vblock_size-offsetn;i++)
    V1block[i+j*vblock_size] = (*V)[i+idv2+j*n_rowwise];
//    V1block[i] = (*V)[i+idv1-offset0+vShft];
  }

  for (j=0;j<m_rowwise;j++) {
#pragma omp parallel for   
  for (i=0;i<offsetn;i++)
    V1block[vblock_size-offsetn+i+j*vblock_size] = LambdaIn[i+lambdaIn_offset+j*offsetn];
  }


  //init Toeplitz arrays
  tpltz_init(vblock_size, lambda[id], &nfft, &blocksize, &T_fft, (T+lambdaShft), &V_fft, &V_rfft, &plan_f, &plan_b);
  //Toeplitz computation
  if(VERBOSE)
    printf("[%d] Before middle-level call : nfft = %d, blocksize = %d\n", rank, nfft, blocksize);
  stmm(&V1block, vblock_size, m_rowwise, 0, vblock_size*m_rowwise, T_fft, lambda[id], V_fft, V_rfft, plan_f, plan_b, blocksize, nfft);

  tpltz_cleanup(&T_fft, &V_fft, &V_rfft, &plan_f, &plan_b);

  for (j=0;j<m_rowwise;j++) {
#pragma omp parallel for                                                 
  for (i=0;i<vnrank_size-offsetn;i++) {
    (*V)[idv2+i+j*n_rowwise] = V1block[i+j*vblock_size];
  }}


  free(V1block);

  }//end of last block
  else { break; }//error  //we can put the generic case here instead of between first and last cases
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  }//end of if(nnew[id]>0)
  }//end of loop over the blocks


  free(LambdaIn);


  return 0;
}


//====================================================================

/// ..Copy a matrix block from an input matrix inside an output matrix. 
/** @ingroup group22
    Copies a matrix block of a size nblockrow x nblockcol from the input matrix Vin
    (size ninrow x nincol) starting with the element (inrow, incol) to the output
    matrix Vout (size notrow x noutcol) starting with the element (outrow, outcol) after
    multiplying by norm. If the output matrix is larger than the block the extra elements
    are either left as they were on the input or zeroed if zero_flag is set to 1. If the 
    block to be copied is larger than either the input or the output matrix an error occurs.
*/
// loop over the blocks to find which are useful for the local data

int get_overlapping_blocks_params(int nbloc, int *idv, int *n, int local_V_size, int nrow, int idp, int *idpnew, int *local_V_size_new, int *nnew, int *ifirstBlock, int *ilastBlock)
{

  int ib, nblockOK=0, firstrow, lastrow, nfullcol_data;
  int idptmp;


//check how many full columns input data have
  nfullcol_data = max(0, (local_V_size-(nrow-idp%nrow)%nrow-(idp+local_V_size)%nrow)/nrow );


  if( nfullcol_data > 0) {

  for( ib=0; ib<nbloc; ib++) {
    if( idv[ib] < nrow) {
      nnew[ib] = min( n[ib], nrow-idv[ib]);  //block used for the product
      nblockOK++;
    }
  }

  }
  else {  //no full column observed 

    firstrow = idp%nrow;
    lastrow = (idp+local_V_size-1)%nrow;

    if( firstrow < lastrow) {  //just one column partially observed   

    for( ib=0; ib<nbloc; ib++) {
    if( (idv[ib]+n[ib] > firstrow) && (idv[ib] < lastrow+1)) {
      nnew[ib] = min( n[ib], nrow-idv[ib]);  //block used for the product
      nblockOK++;
    }
    }

    }
    else {  //two columns partially observed   

      for( ib=0; ib<nbloc; ib++) {
        if( (idv[ib]+n[ib] > firstrow) && (idv[ib] < nrow)) {  //intersects first partial column
          nnew[ib] = min( n[ib], nrow-idv[ib]);  //block used for the product
          nblockOK++;
        }

        if( (idv[ib] < lastrow+1) && (idv[ib]+n[ib] > 0)) {  //intersects second partial column
          nnew[ib] = min( n[ib], nrow-idv[ib]);  //block used for the product
          nblockOK++;  //may overcount but we do not care
        }  //could use else insteed!
      }
     }
  }

  if (VERBOSE)
    printf("nblockOK=%d\n", nblockOK);


  if( nblockOK == 0) return(0);  //no blocks overlapping with the data 

  //find the first and last relevant blocks for the begining and end of the local data  V

 //first block
  idptmp = idp;

  for( *ifirstBlock = -1; *ifirstBlock == -1;     ) {
    for(ib=0;ib<nbloc;ib++) {
      if(nnew[ib] != 0 && idptmp%nrow < idv[ib]+nnew[ib]) break;
    }

    if (ib<nbloc && idv[ib] <= idptmp%nrow) {
      *ifirstBlock = ib;
      *idpnew = idptmp;
    }
    else if (ib<nbloc && idv[ib] > idptmp%nrow) {
      *ifirstBlock = ib;
      int extrabegining = idv[ib]-idp%nrow;
//      *idpnew = idp+extrabegining;//idv[ib];
      int idvfirstcolumn = idptmp/nrow;
      *idpnew = idv[ib]+idvfirstcolumn*nrow;
    }
    else { //ib=nb_blocs
      idptmp += (int) (nrow-idptmp%nrow);
//          idtmp = (int) ceil((1.0*idpnew)/(1.0*nrow))*nrow; // go to the first element of the next column
    }}


 //last block
  idptmp = idp+local_V_size-1;

  for( *ilastBlock = -1; *ilastBlock == -1; ) {
    for(ib=nbloc-1;ib>=0;ib--) {
      if(nnew[ib] != 0 && idv[ib] <= idptmp%nrow) break;
    }


    if (ib>=0 && idptmp%nrow < idv[ib]+nnew[ib]) {
      *ilastBlock = ib;
      *local_V_size_new = local_V_size-(*idpnew)+idp;
    }
    else if (ib>=0 && idv[ib]+nnew[ib] <= idptmp%nrow) {
      *ilastBlock = ib;
      int extraend = (local_V_size-1+idp)%nrow+1-(idv[ib]+nnew[ib]);
      //*local_V_size_new = (local_V_size+idp)%nrow-(idv[*ilastBlock]+nnew[*ilastBlock]);
     //idv[*ilastBlock]+nnew[*ilastBlock]-(*idpnew);
      *local_V_size_new = local_V_size-(*idpnew)+idp-extraend;

      int idvlastcolumn = idptmp/nrow;
      *local_V_size_new = idv[ib]+nnew[ib]+idvlastcolumn*nrow - (*idpnew);

    }
    else {
      idptmp = (int) idptmp - (idptmp%nrow)-1;
//        idtmp = (int) floor( (1.0*idpnew)/(1.0*nrow))*nrow-1; // go to the last element of the previous column
    }}


    return(1);
}




