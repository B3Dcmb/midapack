/** @file   mapmat.c
    @brief  Matrix routines implementation
    @note  Copyright (c) 2010-2012 APC CNRS Université Paris Diderot. This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, see http://www.gnu.org/licenses/lgpl.html

    @note For more information about ANR MIDAS'09 project see http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/index.html
    @note ACKNOWLEDGMENT: This work has been supported in part by the French National  Research Agency (ANR) through COSINUS program (project MIDAS no. ANR-09-COSI-009).
    @author Pierre Cargemel
    @date   November 2011*/
#ifdef W_MPI
 #include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mapmat.h"


/** Create a matrix specifying the number of local rows m,
    the number of non-zero elements per row nnz ,
    indices tab, values tab, ts_flags tab, flag for communication and a communicator comm.
    indices and values tabs must be allocated and contain at least m*nnz elements.
    It represents column indices of the nonzero elements. Respectively values tab
    represents the non-zero values.
    After call MatInit, all precomputation are done and the matrix structure is ready to use.
    That means you can start applying matrix operation.
    Another way to initialize a matrix structure is to apply step by step :
    - MatSetIndices
    - MatSetValues
    - MatLocalShape
    - MatComShape

    @warning do not modify indices tab until you will use the matrix structure.
    @warning [MPI COMM!] with Midapack sequential version, there is no communicator argument.
    @param A pointer to a Mat struct
    @param m number of local rows
    @param nnz number of non-zero per row
    @param indices input tab (modified)
    @param values input tab
    @param flag communication flag
    @param comm MPI communicator
    @ingroup matmap_group11
    @sa MatFree */
int MatInit(Mat *A, int m, int nnz, int *indices, double *values, int flag
#ifdef W_MPI
, MPI_Comm comm
#endif
){
  int err;
  MatSetIndices(A, m, nnz, indices);

  MatSetValues(A, m, nnz, values);


  err = MatLocalShape(A, 3);		// compute lindices (local columns) (method 3 = counting sort)

#ifdef W_MPI
  err = MatComShape(A, flag, comm);		// build communication scheme
#endif
  return err;
}


/** Set column indices of the nonzero elements.
    indices tab must be allocated and contains at least m*nnz elements.
    @param A pointer to a Mat struct
    @param m number of local rows
    @param nnz number of non-zero per row
    @param indices input tab
    @return void
    @ingroup matmap_group11*/
void MatSetIndices(Mat *A, int m, int nnz, int *indices){
  A->m    = m;				// set number of local rows
  A->nnz  = nnz;        		// set number of non-zero values per row
  A->indices = indices;			// point to indices
}


/** Set values of the nonzero elements.
    values tab must be allocated and contains at least m*nnz values.
    @param A pointer to a Mat struct
    @param m number of local rows
    @param nnz number of non-zero per row
    @param values input tab
    @return void
    @ingroup matmap_group11*/
void MatSetValues(Mat *A, int m, int nnz, double *values){
  int err;
  A->m    = m;				// set number of local rows
  A->nnz  = nnz;        		// set number of non-zero values per row
  A->values  = values;			// point to values
}


//===================Part added by Sebastien Cayrols to get amount of memory needed by communication algoritms
void CommInfo(Mat *A){
#if W_MPI
  int i=0,size,rank;
  double maxSizeR = 0.0;
  double maxSizeS = 0.0;
  double amountSizeR = 0.0;
  double amountSizeS = 0.0;
  double stepSum=0.0,stepAvg=0.0;
  //this value is based on data sent
  double *amountSizeByStep = NULL;
  double minStep=0.0, maxStep=0.0;
  double *s=NULL;
  double *r=NULL;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  s = malloc(4*sizeof(double));
  r = malloc(4*3*sizeof(double));
  amountSizeByStep = malloc(A->steps*sizeof(double));
  switch(A->flag){
    case NONE :
      break;
    case BUTTERFLY :
      for(i=0;i<A->steps;i++){
        amountSizeR +=A->nR[i];
        amountSizeS +=A->nS[i];
        if(A->nR[i]>maxSizeR)
          maxSizeR=A->nR[i];
        if(A->nS[i]>maxSizeS)
          maxSizeS=A->nS[i];
      }
    break;
    //==========================Modification added by Sebastien Cayrols : 01/09/2015 , Berkeley
    case BUTTERFLY_BLOCKING_1 :
      for(i=0;i<A->steps;i++){
        amountSizeR +=A->nR[i];
        amountSizeS +=A->nS[i];
        if(A->nR[i]>maxSizeR)
          maxSizeR=A->nR[i];
        if(A->nS[i]>maxSizeS)
          maxSizeS=A->nS[i];
      }
    break;
    case BUTTERFLY_BLOCKING_2 :
      for(i=0;i<A->steps;i++){
        amountSizeR +=A->nR[i];
        amountSizeS +=A->nS[i];
        if(A->nR[i]>maxSizeR)
          maxSizeR=A->nR[i];
        if(A->nS[i]>maxSizeS)
          maxSizeS=A->nS[i];
      }
    break;
    case NOEMPTYSTEPRING :
      for(i=0;i<A->steps;i++){
        amountSizeR +=A->nR[i];
        amountSizeS +=A->nS[i];
        if(A->nR[i]>maxSizeR)
          maxSizeR=A->nR[i];
        if(A->nS[i]>maxSizeS)
          maxSizeS=A->nS[i];
      }
    break;
    //==========================End modification
    case RING :
      for(i=0;i<A->steps;i++){
        amountSizeR +=A->nR[i];
        amountSizeS +=A->nS[i];
        if(A->nR[i]>maxSizeR)
          maxSizeR=A->nR[i];
        if(A->nS[i]>maxSizeS)
          maxSizeS=A->nS[i];
      }
    break;
    case NONBLOCKING :
      for(i=0;i<A->steps;i++){
        amountSizeR +=A->nR[i];
        amountSizeS +=A->nS[i];
        if(A->nR[i]>maxSizeR)
          maxSizeR=A->nR[i];
        if(A->nS[i]>maxSizeS)
          maxSizeS=A->nS[i];
      }
    break;
    case NOEMPTY :
      for(i=0;i<A->steps;i++){
        amountSizeR +=A->nR[i];
        amountSizeS +=A->nS[i];
        if(A->nR[i]>maxSizeR)
          maxSizeR=A->nR[i];
        if(A->nS[i]>maxSizeS)
          maxSizeS=A->nS[i];
      }
    break;
    case ALLTOALLV :          // added -- rs 2015/02/04
	  for(i=0;i<A->steps;i++){
	     amountSizeR +=A->nR[i];
	     amountSizeS +=A->nS[i];
	  }
    break;
    case ALLREDUCE :
      amountSizeR = A->com_count;
      amountSizeS = A->com_count;
      maxSizeR    = A->com_count;
      maxSizeS    = A->com_count;
    break;
  }

  if(A->flag != ALLREDUCE && A->flag != ALLTOALLV ){
    double *t=NULL;

    t=malloc(A->steps*sizeof(double));
    // Copy int array into double array
    for(i=0;i<A->steps;i++)
      t[i]=A->nS[i];

    MPI_Reduce(t,amountSizeByStep,A->steps,MPI_DOUBLE,MPI_SUM,0,comm);

    free(t);

    if(rank==0){
      stepSum=minStep=maxStep=amountSizeByStep[0];
      printf("\n[MEMORY]Step n°%4d, message size : %e",0,amountSizeByStep[0]);
      for(i=1;i<A->steps;i++){
        printf("\n[MEMORY]Step n°%4d, message size : %e",i,amountSizeByStep[i]);
        if(minStep>amountSizeByStep[i])
          minStep=amountSizeByStep[i];
        else if(maxStep<amountSizeByStep[i])
          maxStep=amountSizeByStep[i];
        stepSum+=amountSizeByStep[i];
      }
      stepAvg=stepSum/A->steps;
    }
  }
  s[0]=amountSizeR;
  s[1]=amountSizeS;
  s[2]=maxSizeR;
  s[3]=maxSizeS;
  MPI_Reduce(s,r,4,MPI_DOUBLE,MPI_SUM,0,comm);
  if(rank==0)
    for(i=0;i<4;i++)
      r[i]/=size;
  MPI_Reduce(s,&r[4],4,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(s,&r[8],4,MPI_DOUBLE,MPI_MAX,0,comm);
  if(rank==0){
    printf("\n[MEMORY]Step average             : %e\t[%e,%e]",stepAvg,minStep,maxStep);
    printf("\n[MEMORY]Amount of data received  : %e\t[%e,%e]", r[0],r[4],r[8]);
    printf("\n[MEMORY]Amount of data sent      : %e\t[%e,%e]", r[1],r[5],r[9]);
    printf("\n[MEMORY]Message size received    : %e\t[%e,%e]",r[2],r[6],r[10]);
    printf("\n[MEMORY]Message size sent        : %e\t[%e,%e]\n",r[3],r[7],r[11]);
  }
  free(s);
  free(r);
  free(amountSizeByStep );
#endif
}

void MatReset(Mat *A){
#if W_MPI
  switch(A->flag){
    case NONE :
      break;
    case BUTTERFLY :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case BUTTERFLY_BLOCKING_1 :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case BUTTERFLY_BLOCKING_2 :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case NOEMPTYSTEPRING :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case RING :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case NONBLOCKING :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case NOEMPTY :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case ALLTOALLV :    // added -- rs 2015/02/04
	  free(A->R);		//
	  free(A->nR);		//
	  free(A->S);		//
	  free(A->nS);
    break;
    case ALLREDUCE :
    break;
  }
#endif
}

//===================End

/** Free allocated tabs of a matrix structure including local indices tab and communication tabs.
    Do not free indices and values which user is responsible for.
    @param A pointer to a Mat struct
    @return void
    @sa MatInit MatLocalShape
    @ingroup matmap_group11 */
void MatFree(Mat *A){

  //get information about communication size
  CommInfo(A);

  free(A->lindices);
#if W_MPI
  switch(A->flag){
    case NONE :
      break;
    case BUTTERFLY :
      free(A->com_indices);	//
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    //==========================Modification added by Sebastien Cayrols : 01/09/2015 , Berkeley
    case BUTTERFLY_BLOCKING_1 :
      free(A->com_indices);	//
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case BUTTERFLY_BLOCKING_2 :
      free(A->com_indices);	//
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case NOEMPTYSTEPRING :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    //==========================End modification
    case RING :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case NONBLOCKING :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case NOEMPTY :
      free(A->R);		//
      free(A->nR);		//
      free(A->S);		//
      free(A->nS);
    break;
    case ALLTOALLV :    // Added: rs 2015/02/04
	  free(A->R);		//
	  free(A->nR);		//
	  free(A->S);		//
	  free(A->nS);
    break;
    case ALLREDUCE :
      free(A->com_indices);	//
//===================================Modification from Sebastien Cayrols : comment of these lines to avoid SEGSIGV
//      free(A->R);		//
//      free(A->nR);		//
//      free(A->S);		//
//      free(A->nS);
//===================================End modif
    break;
  }
#endif
}


/** Load matrix from a file.
    This is MatSave dual routine which loads data into matrix reading a specified file (or several specified files).
    Number of files should equal number of processor.
    File format should be ascii files (for examples look at files generated by MapMatSave routine).
    @warning Does not include gap samples flags functionality
    @todo Implement to read several file formats as basic ASCII, XML, HDF5...
    @param mat pointer to the Mat
    @param filename basename of a file, actually data are loaded from  several files denotes by "basename + processor number"
    @return error code
    @ingroup matmap_group11 */
int MatLoad(Mat *mat, char *filename){
  int err;
  int rank;
#if W_MPI
  MPI_Comm_rank(mat->comm, &rank);
#else
  rank=0;
#endif
  FILE *in;
  char fn[100];
  int i=0;
  sprintf(fn, "%s_%d.dat", filename, rank);
  printf("%s", fn);
  in=fopen(fn,"r");
  if(in==NULL){
     printf("cannot open file %s", fn);
     return 1;
  }
  while(feof(in)== 0 && i< (mat->m * mat->nnz)){
    if(mat->nnz==1){
      fscanf(in, "%d %lf", &(mat->indices[i]), &(mat->values[i]));
    }
    else if(mat->nnz==2){
      fscanf(in, "%d %lf %d %lf", &(mat->indices[i]), &(mat->values[i]), &(mat->indices[i+1]), &(mat->values[i+1]));
    }
    else{
      return 1;             //(nnz > 2) not implement
    }
    i+=mat->nnz;
  }
  if(i!= mat->m * mat->nnz ){
    printf("WARNNING data size doesn't fit\n");
  }
  fclose(in);
  return 0;
}



/** Write matrix into files
    This is the dual routine of MatLoad.
    It saves matrix data into files, and can be usefull to check that data are well stored.
    Obviously it performs IO, moreover with one file per processor.
    Therfore, just call this function in a development phase, with few data and few processors.
    @warning Does not include gap samples flags functionality.
    @param A pointer to the Mat
    @param filename file basename, for instance passing "toto" should produce the output files named "toto_$(rank)"
    @return error code
    @ingroup matmap_group11 */
int MatSave(Mat *mat, char *filename){
  FILE *out;
  char fn [100];
  int i,j;
  int rank;
#if W_MPI
  MPI_Comm_rank(mat->comm, &rank);
#else
  sprintf(fn,"%s_%d.dat", filename, rank);
#endif
  out=fopen(fn,"w");
  if(out==NULL){
     printf("cannot open file %s", fn);
     return 1;
  }
  for(i=0; i < (mat->nnz * mat->m); i+=mat->nnz){
    for(j=0; j< mat->nnz ; j++){
      fprintf(out,"%d ",mat->indices[i+j]);
      fprintf(out,"%f ", mat->values[i+j]);
    }
    fprintf(out, "\n");
  }
  fclose(out);
  return 0;
}



/** Compute a local indices into a dense vector, lindices,
    and reindices indices tab according the this local dense vector.
    For this three steps are performed :
    - sort and merge indices tab,
    - allocate lindices of size lcount and copy the sorted indices
    - reindex indices according the local indices

    @warning lindices is internally allocated ( to free it, use MatFree )
    @sa MatComShape MatFree MatSetIndices
    @ingroup matmap_group11 */
int MatLocalShape(Mat *A, int sflag){
  int *tmp_indices;

  tmp_indices = (int *) malloc((int64_t) (A->m) * A->nnz * sizeof(int));	//allocate a tmp copy of indices tab to sort
  memcpy(tmp_indices, A->indices, (int64_t) (A->m) * A->nnz * sizeof(int));	//copy

//  A->lcount = omp_psort(tmp_indices, A->m * A->nnz, sflag);		//sequential sort tmp_indices
  A->lcount = ssort(tmp_indices, A->m * A->nnz, sflag);		//sequential sort tmp_indices

  A->lindices = (int *) malloc( A->lcount * sizeof(int));
  memcpy(A->lindices, tmp_indices, A->lcount * sizeof(int));	//copy tmp_indices into lindices and free
  free(tmp_indices);

  sindex(A->lindices, A->lcount, A->indices, A->nnz * A->m);
  return 0;
}





#if W_MPI
/** Transform the matrix data structure, identifying columns shared by several processors
    @warning [MPI ONLY!] this function does not exist in Midapack sequential version
    @sa MatLocalShape MatInit TrMatVecProd
    @ingroup matmap_group11 */
int MatComShape(Mat *A, int flag, MPI_Comm comm){
  int size;
  int i, min, max, j;
  A->comm = comm;			// set communivcator
  A->flag=flag;
  MPI_Comm_size(A->comm, &size);
  if((A->flag==BUTTERFLY || A->flag==BUTTERFLY_BLOCKING_1 || A->flag==BUTTERFLY_BLOCKING_2) && is_pow_2(size)!=0)
    A->flag=RING;
  switch(A->flag){
    case BUTTERFLY :
      A->steps = log_2(size);
      A->S = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate sending maps tab
      A->R = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate receiving maps tab
      A->nS = (int* ) malloc(A->steps * sizeof(int));                   //allocate sending map sizes tab
      A->nR = (int* ) malloc(A->steps * sizeof(int));                   //allocate receiving map size tab
      butterfly_init(A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), A->R, A->nR, A->S, A->nS, &(A->com_indices), &(A->com_count), A->steps, A->comm);
      break;
    //==========================Modification added by Sebastien Cayrols : 01/09/2015 , Berkeley
    case BUTTERFLY_BLOCKING_1 :
      A->steps = log_2(size);
      A->S = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate sending maps tab
      A->R = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate receiving maps tab
      A->nS = (int* ) malloc(A->steps * sizeof(int));                   //allocate sending map sizes tab
      A->nR = (int* ) malloc(A->steps * sizeof(int));                   //allocate receiving map size tab
      butterfly_init(A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), A->R, A->nR, A->S, A->nS, &(A->com_indices), &(A->com_count), A->steps, A->comm);
      break;
    case BUTTERFLY_BLOCKING_2 :
      A->steps = log_2(size);
      A->S = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate sending maps tab
      A->R = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate receiving maps tab
      A->nS = (int* ) malloc(A->steps * sizeof(int));                   //allocate sending map sizes tab
      A->nR = (int* ) malloc(A->steps * sizeof(int));                   //allocate receiving map size tab
      butterfly_init(A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), A->R, A->nR, A->S, A->nS, &(A->com_indices), &(A->com_count), A->steps, A->comm);
      break;
    case NOEMPTYSTEPRING :
      A->steps = size;
      A->S = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate sending maps tab
      A->R = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate receiving maps tab
      A->nS = (int* ) malloc(A->steps * sizeof(int));                   //allocate sending map sizes tab
      A->nR = (int* ) malloc(A->steps * sizeof(int));                   //allocate receiving map size tab
      ring_init(A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), A->R, A->nR, A->S, A->nS, A->steps, A->comm);
      A->com_count = A->lcount-(A->nnz)*(A->trash_pix);
      A->com_indices = A->lindices+(A->nnz)*(A->trash_pix);
      break;
    //==========================End modification
    case RING :
      A->steps = size;
      A->S = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate sending maps tab
      A->R = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate receiving maps tab
      A->nS = (int* ) malloc(A->steps * sizeof(int));                   //allocate sending map sizes tab
      A->nR = (int* ) malloc(A->steps * sizeof(int));                   //allocate receiving map size tab
      ring_init(A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), A->R, A->nR, A->S, A->nS, A->steps, A->comm);
      A->com_count = A->lcount-(A->nnz)*(A->trash_pix);
      A->com_indices = A->lindices+(A->nnz)*(A->trash_pix);
      break;
    case NONBLOCKING :
      A->steps = size;
      A->S = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate sending maps tab
      A->R = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate receiving maps tab
      A->nS = (int* ) malloc(A->steps * sizeof(int));                   //allocate sending map sizes tab
      A->nR = (int* ) malloc(A->steps * sizeof(int));                   //allocate receiving map size tab
      ring_init(A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), A->R, A->nR, A->S, A->nS, A->steps, A->comm);
      A->com_count = A->lcount-(A->nnz)*(A->trash_pix);
      A->com_indices = A->lindices+(A->nnz)*(A->trash_pix);
      break;
    case NOEMPTY :
      A->steps = size;
      A->S = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate sending maps tab
      A->R = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate receiving maps tab
      A->nS = (int* ) malloc(A->steps * sizeof(int));                   //allocate sending map sizes tab
      A->nR = (int* ) malloc(A->steps * sizeof(int));                   //allocate receiving map size tab
      ring_init(A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), A->R, A->nR, A->S, A->nS, A->steps, A->comm);
      A->com_count = A->lcount-(A->nnz)*(A->trash_pix);
      A->com_indices = A->lindices+(A->nnz)*(A->trash_pix);
      break;
    case ALLTOALLV :
      A->steps = size;
      A->S = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate sending maps tab
      A->R = (int** ) malloc(A->steps * sizeof(int* ));                 //allocate receiving maps tab
      A->nS = (int* ) malloc(A->steps * sizeof(int));                   //allocate sending map sizes tab
      A->nR = (int* ) malloc(A->steps * sizeof(int));                   //allocate receiving map size tab
      ring_init(A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), A->R, A->nR, A->S, A->nS, A->steps, A->comm);
      A->com_count = A->lcount-(A->nnz)*(A->trash_pix);
      A->com_indices = A->lindices+(A->nnz)*(A->trash_pix);
      break;
    case ALLREDUCE :
      MPI_Allreduce(&A->lindices[A->lcount-1], &max, 1, MPI_INT, MPI_MAX, A->comm);	//maximum index
      MPI_Allreduce(&A->lindices[(A->nnz)*(A->trash_pix)], &min, 1, MPI_INT, MPI_MIN, A->comm);	//
      A->com_count=(max-min+1);
      A->com_indices = (int *) malloc((A->lcount-(A->nnz)*(A->trash_pix)) * sizeof(int)); //warning
      i=(A->nnz)*(A->trash_pix);
      j=0;
      while( j<A->com_count && i<A->lcount){ //same as subsetmap for a coutiguous set
        if(min+j < A->lindices[i]){
          j++;
        }
        else{
          A->com_indices[i-(A->nnz)*(A->trash_pix)]=j;
          i++;
          j++;
        }
      }
      break;
  }
 return 0;
}
#endif


/** Perform matrix-vector multiplication, \f$y \leftarrow A x\f$.
    @param A pointer to a Mat
    @param x input vector (overlapped)
    @param y output vector (distributed)
    @ingroup matmap_group11
    @ingroup matmap_group12a */
int MatVecProd(Mat *A, double *x, double* y, int pflag){
  int i, j, e;
  for(i=0; i<A->m; i++) 					//
      y[i] = 0.0;

  e=0;
  if(A->trash_pix){
    for(i=0; i<A->m*A->nnz; i+=A->nnz){
      if(A->indices[i]!=0){
        for(j=0; j<A->nnz; j++){					//
          y[e] += A->values[i+j] * x[A->indices[i+j]-(A->nnz)];
        }
      }
      e++;
    }
  }
  else{
    for(i=0; i<A->m*A->nnz; i+=A->nnz){
      for(j=0; j<A->nnz; j++){					//
        y[e] += A->values[i+j] * x[A->indices[i+j]];
      }
      e++;
    }
  }
  return 0;
};


#ifdef W_MPI
/** Perform transposed matrix-vector multiplication, \f$x \leftarrow A^t y\f$.
    This naive version does not require a precomputed communication structure.
    But communication volumes may be significant.
    Consequently in most of the cases is not optimized.
    @warning [MPI ONLY!] this function does not exist in Midapack sequential version
    @sa TrMatVecProd MatLocalShape
    @param mat pointer
    @param y local input vector (distributed)
    @param x local output vector (overlapped)
    @ingroup matmap_group11
    @ingroup matmap_group12b */
int TrMatVecProd_Naive(Mat *A, double *y, double* x, int pflag){
  int i, j, e, rank, size;
  int *rbuf, rbufcount;
  double *rbufvalues, *lvalues;
  int p, rp, sp, tag;
  MPI_Request s_request, r_request;
  MPI_Status status;

  MPI_Comm_rank(A->comm, &rank);				//get rank and size of the communicator
  MPI_Comm_size(A->comm, &size);				//
  lvalues = (double *) malloc( A->lcount *sizeof(double));	//allocate and set local values to 0.0
  for(i=0; i < A->lcount; i++)				//
    lvalues[i]=0.0;						//

  e=0;
  for(i=0; i<A->m; i++){					//local transform reduces
    for(j=0; j<A->nnz; j++){					//
      lvalues[A->indices[i*A->nnz+j]] += (A->values[i*A->nnz+j]) * y[i];
    }
  }

  memcpy(x, lvalues, (A->lcount)*sizeof(double)); 			//copy local values into the result*/
  MPI_Allreduce(&(A->lcount), &(rbufcount), 1, MPI_INT, MPI_MAX, A->comm);	//find the max communication buffer sizes, and allocate

  rbuf = (int *)    malloc(rbufcount * sizeof(int));
  rbufvalues  = (double *) malloc(rbufcount * sizeof(double));

  tag=0;
  for (p=1; p < size; p++){	//loop : collective global reduce in ring-like fashion
    rp = (size + rank - p)%size;
    sp = (rank + p)%size;
    MPI_Send(&(A->lcount), 1, MPI_INT, sp, 0, A->comm);				//exchange sizes
    MPI_Recv(&rbufcount, 1, MPI_INT, rp, 0, A->comm, &status);
    tag++;
    MPI_Irecv(rbuf, rbufcount, MPI_INT, rp, tag, A->comm, &r_request);		//exchange local indices
    MPI_Isend(A->lindices, A->lcount, MPI_INT, sp, tag, A->comm, &s_request);
    MPI_Wait(&r_request, &status);
    MPI_Wait(&s_request, &status);
    tag++;
    MPI_Irecv(rbufvalues, rbufcount, MPI_DOUBLE, rp, tag, A->comm, &r_request);	//exchange local values
    MPI_Isend(lvalues, A->lcount, MPI_DOUBLE, sp, tag, A->comm, &s_request);
    tag++;
    MPI_Wait(&r_request, &status);
    m2m_sum(rbufvalues, rbuf, rbufcount, x, A->lindices, A->lcount);		//sum in the result
    MPI_Wait(&s_request, &status);
  }
  free(lvalues);
  return 0;
}
#endif


/** Perform a transposed matrix-vector multiplication, \f$x \leftarrow A^t y\f$
    using a precomputed communication scheme. Before calling this routine,
    the communication structure should have been set, calling MatInit or MatComShape.
    The routine can be divided in two steps :
    - a local matrix vector multiplication
    - a collective-reduce. it consits in a sum reduce over all processes.

    The collective reduce is performed using algorithm previously defined : ring, butterfly ...
    @sa MatVecProd MatComShape TrMatVecProd_Naive MatInit
    @param A a pointer to a Mat
    @param y local input vector (distributed)
    @param x local output vector (overlapped)
    @ingroup matmap_group11 */
int TrMatVecProd(Mat *A, double *y, double* x, int pflag){
  double *sbuf, *rbuf;
  int i, j, k, e;
  int nSmax, nRmax;
  double *lvalues;
  if(A->trash_pix){
    for(i=0; i < A->lcount-A->nnz; i++)				//refresh vector
      x[i]=0.0;						//

      e=0;
      for(i=0; i< A->m*A->nnz; i+=A->nnz){
        if(A->indices[i]!=0){
          //local transform reduce
          for(j=0; j< A->nnz; j++){
            x[A->indices[i+j]-(A->nnz)] += A->values[i+j] * y[e];	//
          }
        }							//
        e++;
      }
  }
  else{
    for(i=0; i < A->lcount; i++)				//refresh vector
      x[i]=0.0;						//

      e=0;
      for(i=0; i< A->m*A->nnz; i+=A->nnz){//local transform reduce
        for(j=0; j< A->nnz; j++){
          x[A->indices[i+j]] += A->values[i+j] * y[e];	//
        }						//
        e++;
      }
  }


#ifdef W_MPI
  greedyreduce(A, x);					//global reduce
#endif
  return 0;
}






#ifdef W_MPI
/** @brief Print information about a matrix.

    @n Usefull function to check, debug or bench. It prints matrix array sizes.
    @sa MatSave
    @param A pointer to the Mat
    @ingroup matmap_group11*/
int MatInfo(Mat *mat, int verbose, char *filename){
  FILE *out;
  int *n;
  int *sr;
  int *s;
  int nnzline, sparsity, maxstep, maxsize, sumline, total;
  int i, j, k;
  char fn [100];
  int rank, size;
  int master=0;
  MPI_Comm_rank(mat->comm, &rank);
  MPI_Comm_size(mat->comm, &size);


  if(rank==master){			//master process saves data into filename_info.txt
    sprintf(fn,"%s_%s", filename, "info.txt");
    out=fopen(fn,"w");
    if(out==NULL){
      printf("cannot open file %s\n", fn);
      return 1;
    }
    printf("open file %s ...", fn);
    fprintf(out, "flag %d\n", mat->flag);	//print matirx main description : flag (communication scheme),
    fprintf(out, "rows %d\n ", mat->m);	//rows per process,
    fprintf(out, "nnz %d\n", mat->nnz);	//nnz (number of non zero per row).
    fprintf(out, "\n");	//separator
  }

  /*n = (int* ) calloc(mat->lcount,sizeof(int));		//allocate
  //printf("before gather %d\n", rank);
  MPI_Gather(&(mat->lcount), 1, MPI_INT, n, 1, MPI_INT, master, mat->comm);		//gather nbnonempty cols
  //printf("after gather %d\n", rank);

  if(rank==master){			//master process saves data into filename_info.txt
    fprintf(out, "cols :\n");	//nnz (number of non zero per row).
    for(i=0; i<size; i++)		//
      fprintf(out, "%d ", n[i]);	//non-empty columns per process.
    fprintf(out, "\n");			//
  }
  free(n); */				//free allocated tabs

  nnzline = 0;				//compute communication sparsity and maximum message size
  sumline = 0;
  for(i=0; i<mat->steps; i++){		//
    sumline+=mat->nS[i];
    if(mat->nS[i]==0){			//
      nnzline +=1;			//
    }					//
  }					//
  MPI_Reduce(&nnzline, &sparsity, 1, MPI_INT, MPI_SUM, 0, mat->comm);	//sparsity
  MPI_Reduce(&sumline, &total, 1, MPI_INT, MPI_SUM, 0, mat->comm);	//sparsity
  if(rank==master){         		//master process saves data into filename_info.txt
    fprintf(out, "sparsity %d\n", sparsity);	//
    fprintf(out, "total %d\n", total);	//
  }

  maxsize =0;
  for(i=0; i<mat->steps; i++){		//
    MPI_Reduce(&(mat->nS[i]), &maxstep, 1, MPI_INT, MPI_MAX, 0, mat->comm);	//maximum message size
    maxsize+=maxstep;			//
  }    				//
  if(rank==master){				//master process saves data into filename_info.txt
    fprintf(out, "maxsize %d\n ", maxsize);	//
    fprintf(out, "\n");				//separator
  }						//

 /* s = (int* ) calloc((mat->steps),sizeof(int));	//allocate steps
  MPI_Reduce(mat->nS, s, mat->steps, MPI_INT, MPI_SUM, 0, mat->comm);	//imaximum message size

  if(rank==master){			//master process saves data into filename_info.txt
    fprintf(out, "sumsteps :\n");	//nnz (number of non zero per row).
    for(i=0; i<mat->steps; i++)		//
      fprintf(out, "%d ", s[i]);	//non-empty columns per process.
    fprintf(out, "\n");			//
  }
  free(s);

  if(verbose==1){
    sr = (int* ) calloc((mat->steps)*size,sizeof(int));	//allocate send/receive matrix
    //printf("before gather %d\n", rank);
    MPI_Gather(mat->nS, mat->steps, MPI_INT, sr, mat->steps, MPI_INT, master, mat->comm);	//gather nbnonempty cols
    //printf("after gather %d\n", rank);

    if(rank==master){			//master process saves data into filename_info.txt
    fprintf(out, "send/receive matrix\n");	//separator
      for(i=0; i<size; i++){ 		//print collective description :
        if(mat->flag==BUTTERFLY){		//send-receive matrix
          for(j=0; j<size; j++){ 		//print send/receive matrix
            if(j>i){
              if(is_pow_2(j-i)==0)
                fprintf(out,"%d ", sr[i*(mat->steps)+log_2(j-i)]);
              else
                fprintf(out,"%d ", 0);
            }
            else if(i>j){
              if(is_pow_2(size+j-i)==0)
                fprintf(out,"%d ", sr[i*(mat->steps)+log_2(size+j-i)]);
              else
                fprintf(out,"%d ", 0);
            }
            else{
              fprintf(out,"%d ", 0);
            }
          }
          fprintf(out, "\n");
        }
        else{
          for(j=0; j<size; j++){ 		//print send/receive matrix
            if(j>i){
              fprintf(out,"%d ", sr[i*(mat->steps)+j-i]);
            }
            else if(i>j){
              fprintf(out,"%d ", sr[(i+1)*(mat->steps)-i+j]);
            }
            else{
              fprintf(out,"%d ", 0);
            }
          }
          fprintf(out, "\n");
        }
      }
    }
    free(sr);
  }*/

  if(rank==master){			//master process saves data into filename_info.txt
    fclose(out);
    printf("close %s\n", fn);
  }
  return 0;
}
#endif



#if W_MPI
int greedyreduce(Mat *A, double* x){
  int i, j, k;
  int nSmax, nRmax, nStot, nRtot;
  double *lvalues;
  lvalues = (double *) malloc((A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));	//allocate and set to 0.0 local values
  memcpy(lvalues, x, (A->lcount-(A->nnz)*(A->trash_pix)) *sizeof(double));		//copy local values into result values
  double *com_val;
  double *out_val;
  int ne=0;
  switch(A->flag){
    case BUTTERFLY :
      for(k=0; k< A->steps; k++)                                  //compute max communication buffer size
        if(A->nR[k] > nRmax)
          nRmax = A->nR[k];
      for(k=0; k< A->steps; k++)
        if(A->nS[k] > nSmax)
          nSmax = A->nS[k];
      com_val=(double *) malloc( A->com_count *sizeof(double));
      for(i=0; i < A->com_count; i++)
        com_val[i]=0.0;
      m2m(lvalues, A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), com_val, A->com_indices, A->com_count);
      butterfly_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
      m2m(com_val, A->com_indices, A->com_count, x, A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix));
      free(com_val);
      break;
    //==========================Modification added by Sebastien Cayrols : 01/09/2015 , Berkeley
    case BUTTERFLY_BLOCKING_1 :
      for(k=0; k< A->steps; k++)                                  //compute max communication buffer size
        if(A->nR[k] > nRmax)
          nRmax = A->nR[k];
      for(k=0; k< A->steps; k++)
        if(A->nS[k] > nSmax)
          nSmax = A->nS[k];
      com_val=(double *) malloc( A->com_count *sizeof(double));
      for(i=0; i < A->com_count; i++)
        com_val[i]=0.0;
      m2m(lvalues, A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), com_val, A->com_indices, A->com_count);
      butterfly_blocking_1instr_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
      m2m(com_val, A->com_indices, A->com_count, x, A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix));
      free(com_val);
      break;
    case BUTTERFLY_BLOCKING_2 :
      for(k=0; k< A->steps; k++)                                  //compute max communication buffer size
        if(A->nR[k] > nRmax)
          nRmax = A->nR[k];
      for(k=0; k< A->steps; k++)
        if(A->nS[k] > nSmax)
          nSmax = A->nS[k];
      com_val=(double *) malloc( A->com_count *sizeof(double));
      for(i=0; i < A->com_count; i++)
        com_val[i]=0.0;
      m2m(lvalues, A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix), com_val, A->com_indices, A->com_count);
      butterfly_blocking_1instr_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, com_val, A->steps, A->comm);
      m2m(com_val, A->com_indices, A->com_count, x, A->lindices+(A->nnz)*(A->trash_pix), A->lcount-(A->nnz)*(A->trash_pix));
      free(com_val);
      break;
    case NOEMPTYSTEPRING :
      for(k=1; k< A->steps; k++)				//compute max communication buffer size
        if(A->nR[k] > nRmax)
          nRmax = A->nR[k];
      nSmax = nRmax;
      ring_noempty_step_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, x, A->steps, A->comm);
      break;
    //==========================End modification
    case RING :
      for(k=1; k< A->steps; k++)				//compute max communication buffer size
        if(A->nR[k] > nRmax)
          nRmax = A->nR[k];
      nSmax = nRmax;
      ring_reduce(A->R, A->nR, nRmax, A->S, A->nS, nSmax, lvalues, x, A->steps, A->comm);
      break;
    case NONBLOCKING :
      ring_nonblocking_reduce(A->R, A->nR, A->S, A->nS, lvalues, x, A->steps, A->comm);
      break;
    case NOEMPTY :
      for(k=1; k< A->steps; k++)
        if(A->nR[k]!=0)
          ne++;
      ring_noempty_reduce(A->R, A->nR, ne, A->S, A->nS, ne, lvalues, x, A->steps, A->comm);
      break;
    case ALLREDUCE :
      com_val=(double *) malloc( A->com_count *sizeof(double));
      out_val=(double *) malloc( A->com_count *sizeof(double));
      for(i=0; i < A->com_count; i++){
        com_val[i]=0.0;
        out_val[i]=0.0;
      }
      s2m(com_val, lvalues, A->com_indices, A->lcount-(A->nnz)*(A->trash_pix));
      /*for(i=0; i < A->com_count; i++){
         printf("%lf ", com_val[i]);
      } */
      MPI_Allreduce(com_val, out_val, A->com_count, MPI_DOUBLE, MPI_SUM, A->comm);	//maximum index
      /*for(i=0; i < A->com_count; i++){
         printf("%lf ", out_val[i]);
      } */
      m2s(out_val, x, A->com_indices, A->lcount-(A->nnz)*(A->trash_pix));                                 //sum receive buffer into values
      free(com_val);
      free(out_val);
      break;
    case ALLTOALLV :
 	    nRtot=nStot=0;
	    for(k=0; k< A->steps; k++){				//compute buffer sizes
	       nRtot += A->nR[k];                // to receive
	       nStot += A->nS[k];                // to send
       }

      alltoallv_reduce(A->R, A->nR, nRtot, A->S, A->nS, nStot, lvalues, x, A->steps, A->comm);
      break;
  }
  free(lvalues);
  return 0;
}
#endif
