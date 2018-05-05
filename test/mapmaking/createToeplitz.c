// Create Toeplitz - fd 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "midapack.h"
#include <time.h>


char CHAR_RW='\0';  //global variable for write mode


int defineTpltz_avg( Tpltz *Nm1, int nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks_loc, int nb_blocks_tot, int idp, int local_V_size, Flag flag_stgy, MPI_Comm comm) 
{

// faire les allocs ici avec la structure Tpltz

  Nm1->nrow = nrow; //glob //recup du fichier params apres (en variables globales)
  Nm1->m_cw = m_cw; //glob
  Nm1->m_rw = m_rw; //glob
  Nm1->tpltzblocks = tpltzblocks; //toep
  Nm1->nb_blocks_loc = nb_blocks_loc; //toep
  Nm1->nb_blocks_tot = nb_blocks_tot;  //toep
  Nm1->idp = idp; //comput
  Nm1->local_V_size = local_V_size; //comput
  Nm1->flag_stgy = flag_stgy; //param
  Nm1->comm = comm; //param


  return 0;
};




int defineBlocks_avg(Block *tpltzblocks, double *T, int nb_blocks_loc, int n_block_avg, int lambda_block_avg, int64_t id0 )
{

int i;


  for ( i=0; i<nb_blocks_loc; i++)
    tpltzblocks[i].n = n_block_avg;

  for ( i=0; i<nb_blocks_loc; i++)
    tpltzblocks[i].lambda = lambda_block_avg;

  tpltzblocks[0].idv = (int64_t) (id0/n_block_avg) * n_block_avg ;
  for(i=1;i<nb_blocks_loc;i++)
    tpltzblocks[i].idv = (int64_t) tpltzblocks[i-1].idv + tpltzblocks[i-1].n;

  for( i=0; i<nb_blocks_loc; i++) {
    tpltzblocks[i].T_block = (T);
  }

  for(i=0; i<nb_blocks_loc; i++)  
    printf("tpltzblocks[%d].idv=%ld\n", i, tpltzblocks[i].idv);

  return 0;
}



//Create T from input files
int createTfromfiles(double *T, int Tsize)
{

//  oo



  return 0;
}



int createT(double *T, int Tsize)
{

  int i;
  //srand (time (NULL));  //init seed
  srand (Tsize); 

  if (1==0) {
  //input matrix definition of T
    for(i=0;i<Tsize;i++)
      T[i]= 1.0 + rand()/((double) RAND_MAX);


  }
  else if (1==1) {

  //input matrix definition of T
    for(i=0;i<Tsize;i++) {
      if (i == 0) {
        T[i]=10.;}
      else if (i == 1) {
        T[i]=2.;}
      else if (i == 2) {
        T[i]=3.;}
      else {
        T[i]=rand()/((double) RAND_MAX);
     }}
 }
  else if (1==0) {

  //input matrix definition of T
    for(i=0;i<Tsize;i++) {
      if (i == 0) {
        T[i]=2.;}
      else if (i == 1) {
        T[i]=-1.;}
      else if (i == 2) {
        T[i]=0.;}
      else {
        T[i]=0.;//rand()/((double) RAND_MAX);
     }}
 }
  else {

//#include "invtt_params.h"

  double *invtt;

  T = invtt;
//  createinvtt(invtt);

  } //end if


  return 0;
}

