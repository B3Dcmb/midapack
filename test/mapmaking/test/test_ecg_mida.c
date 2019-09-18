/**
 * \file    test_ecg_prealps_op.c
 * \author  Olivier Tissot
 * \date    2016/06/23
 * \brief   Example of usage of E(nlarged) C(onjugate) G(radient)
 */

/******************************************************************************/
/*                                  INCLUDE                                   */
/******************************************************************************/
/* STD */
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
/* MPI */
#include <mpi.h>
/* MKL */
#include <mkl.h>

/* CPaLAMeM */
//#include <cpalamem_macro.h>
//#include <cpalamem_instrumentation.h>

/* midapack */
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "midapack.h"

/* preAlps */
#include "overlap_ecg.h"

/* Command line parser */
#include <ctype.h>
#include <getopt.h>

/******************************************************************************/
/* midapack cluster and problem options */

double Opmmpreconditioner(Mat A, Mat BJ,  double *X, double *Y, int ncol);
double Opmmmatrix(Mat A, Tpltz Nm1, double *X, double *Y, int ncol);


extern const char *WORKDIR="/global/cscratch1/sd/elbouha/data_TOAST/test4_clean/";
double FKNEE=1.0;//0.25;

/******************************************************************************/
/*                                    CODE                                    */
/******************************************************************************/
/* Private function to print the help message */
void _print_help() {
  printf("DESCRIPTION\n");
  printf("\tSolves Ax = b using a Parallel Enlarged Conjugate Gradient."
          " A must be symmetric positive definite.\n");
  printf("USAGE\n");
  printf("\tmpirun -n nb_proc"
         " ./test_ecg_prealps_op"
         " -e/--enlarging-factor int"
         " [-h/--help]"
         " [-i/--iteration-maximum int]"
         " -o/--ortho-alg int"
         " -r/--search-dir-red int"
         " [-t/--tolerance double]\n");
  printf("OPTIONS\n");
  printf("\t-e/--enlarging-factor : enlarging factor"
                                  " (cannot exceed nprocs)\n");
  printf("\t-h/--help             : print this help message\n");
  printf("\t-i/--iteration-maximum: maximum of iteration count"
                                  " (default is 1000)\n");
  printf("\t-o/--ortho-alg        : orthogonalization scheme"
                                  " (0: odir, 1: omin)\n");
  printf("\t-r/--search-dir-red   : adaptive reduction of the search"
                                  " directions (0: no, 1: yes)\n");
  printf("\t-t/--tolerance        : tolerance of the method"
                                  " (default is 1e-5)\n");
}

/* midapack private function */
void usage(){
    printf("usage...\n");
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  CPLM_SetEnv();

  /*================ Command line parser ================*/
  int c;
  static struct option long_options[] = {
    {"enlarging-factor" , required_argument, NULL, 'e'},
    {"help"             , no_argument      , NULL, 'h'},
    {"iteration-maximum", optional_argument, NULL, 'i'},
    {"ortho-alg"        , required_argument, NULL, 'o'},
    {"search-dir-red"   , required_argument, NULL, 'r'},
    {"tolerance"        , optional_argument, NULL, 't'},
    {NULL               , 0                , NULL, 0}
  };

  int opterr = 0;
  int option_index = 0;

  // Set global parameters for ECG
  double tol = 1e-6;
  int maxIter = 1000;
  int enlFac = 1, ortho_alg = 0, bs_red = 0;

  while ((c = getopt_long(argc, argv, "e:hi:o:r:t:", long_options, &option_index)) != -1)
    switch (c) {
      case 'e':
        enlFac = atoi(optarg);
        break;
      case 'h':
        _print_help();
        MPI_Abort(MPI_COMM_WORLD, opterr);
      case 'i':
        if (optarg != NULL)
          maxIter = atoi(optarg);
        break;
      case 'o':
        ortho_alg = atoi(optarg);
        break;
      case 'r':
        bs_red = atoi(optarg);
        break;
      case 't':
        if (optarg != NULL)
          tol = atof(optarg);
        break;
      case '?':
        if (optopt == 'e'
            || optopt == 'i'
            || optopt == 'o'
            || optopt == 'r'
            || optopt == 't')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        _print_help();
        MPI_Abort(MPI_COMM_WORLD, opterr);
      default:
        MPI_Abort(MPI_COMM_WORLD, opterr);
    }

// CPLM_OPEN_TIMER
  /*================ Initialize ================*/
  int rank, size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //printf("rank = %d, size = %d \n",rank,size);

  // Force sequential execution on each MPI process
  // OT: I tested and it still works with OpenMP activated
  MKL_Set_Num_Threads(1);

  /*======== midapack initialization ==========*/
  int64_t	M;       //Global number of rows
  int 	N, Nnz;  //of columns, of non-zeros values per column for the pointing matrix A
  int		m, n;  //local number of rows, of columns for the pointing matrix A
  int64_t	gif;			//global indice for the first local line
  int		i, j, k;

  int 		*indices, *id0pix, *ll;
  double 	*values;
  int 		pointing_commflag ;	//option for the communication scheme for the pointing matrix
  double	*b, *Ag, *Ad, *pol_ang, *wghts; 	 	//temporal domain vectors
  double	*x, *g, *d, *Ax_b, *tmp;	//pixel domain vectors
  double        alpha, beta, gamma, resold, resnew;
  double 	localreduce;
  double	st, t;		 	//timer, start time
  int 		output, timer, info;

//communication scheme for the pointing matrix  (to move in .h)
  pointing_commflag=6; //2==BUTTERFLY - 1==RING

//global data caracteristics
  int Nb_t_Intervals = 128;           //total number of stationnary intervals
  int t_Interval_length = 1749900;          //length for each stationnary interval
  int t_Interval_length_true = 17499;
  int LambdaBlock = pow(2,0);//pow(2,14)+1;  //lambda length for each stationnary interval

  double fknee = FKNEE; //0.25;
  Nnz = 3;

  //Number of loop we need to read the all t_Interval_length
    int t_Interval_loop = t_Interval_length/t_Interval_length_true ;
    //printf("[rank %d] t_Interval_loop=%d\n", rank, t_Interval_loop );

  //total length of the time domaine signal
    M = (int64_t) Nb_t_Intervals*t_Interval_length ;
    N = 0;

    //printf("[rank %d] M=%ld\n", rank, M);

  //compute distribution indexes over the processes
    partition(&gif, &m, M, rank, size);

    double Nb_t_Intervals_loc_dble = Nb_t_Intervals/size;
    int Nb_t_Intervals_loc = ceil( Nb_t_Intervals*1.0/size );
    int nb_proc_shared_one_interval = max(1, size/Nb_t_Intervals ); //same as ceil( (size*1.0)/Nb_t_Intervals );
    //printf("[rank %d] nb_proc_shared_one_interval=%d\n", rank, nb_proc_shared_one_interval );

    int t_Interval_length_loc = t_Interval_length/nb_proc_shared_one_interval; //just to check
  //should be equal to min(m ,t_Interval_length)

    //printf("[rank %d] size=%d \t m=%d \t Nb_t_Intervals=%d \t t_Interval_length=%d\n", rank, size, m, Nb_t_Intervals, t_Interval_length );
    //printf("[rank %d] Nb_t_Intervals_loc=%d \t t_Interval_length_loc=%d\n", rank, Nb_t_Intervals_loc , t_Interval_length_loc);

  //input data memory allocation
  indices  = (int *) malloc(Nnz*m * sizeof(int));     //for pointing matrix indices
  b   = (double *) malloc(m * sizeof(double));    //full raw data vector for the signal
  wghts = (double *) malloc(Nnz*m * sizeof(double)); //for poiting matrix weights


//Read data from files:
//note: work only if the number of processes is a multiple of the number of stationary intervals

//Definition for the input data
  int part_id;      // stationnaly period id number
  int *point_data;  // scann strategy input data for the pointing matrix
  double *signal;   // signal input data
  double *weights; // weights of the pointing matrix

  int number_in_interval;

  int flag_bigdata=1;
  st=MPI_Wtime();

  if (flag_bigdata==1 && nb_proc_shared_one_interval>1) {

    if (rank==0) {
      printf("#######  ENTER BiG DATA MODE   ################\n");
    }
      part_id = rank/nb_proc_shared_one_interval; //floor(rank/nb_proc_shared_one_interval);
      number_in_interval = rank%nb_proc_shared_one_interval;
      printf("[rank %d] part_id=%d\n", rank, part_id );  //interval id number
      printf("[rank %d] number_in_interval=%d\n", rank, number_in_interval );

      int t_Interval_loop_loc = ceil( t_Interval_loop*1.0/nb_proc_shared_one_interval);
      printf("[rank %d] t_Interval_loop_loc=%d\n", rank, t_Interval_loop_loc );
      printf("[rank %d] m=%d \t t_Interval_length_true*t_Interval_loop_loc=%d\n", rank, m, t_Interval_length_true*t_Interval_loop_loc );

      point_data  = (int *) malloc(Nnz*t_Interval_length_true*t_Interval_loop_loc * sizeof(int));
      signal      = (double *) malloc(t_Interval_length_true*t_Interval_loop_loc  * sizeof(double));
      weights = (double *) malloc(Nnz*t_Interval_length_true*t_Interval_loop_loc * sizeof(double));

      int jump=0;
      if(t_Interval_loop_loc==1){
        jump = t_Interval_length_true*floor(number_in_interval/t_Interval_loop);
      }
      else
        jump = t_Interval_length_true*t_Interval_loop_loc*number_in_interval;
      for (i=0; i < t_Interval_loop_loc; ++i) {
        ioReadTOAST_data(jump, i, t_Interval_length_true, part_id, point_data+t_Interval_length_true*Nnz*i, signal+t_Interval_length_true*i, weights+t_Interval_length_true*Nnz*i);
      }
      //just keep the relevant part of the stationary interval for the local process
      int nb_proc_shared_one_subinterval = max(1, size/(Nb_t_Intervals*t_Interval_loop) );
      //same as ceil( (size*1.0)/Nb_t_Intervals );
      int number_in_subinterval = rank%nb_proc_shared_one_subinterval;

      int t_Interval_length_subinterval_loc = t_Interval_length_true/nb_proc_shared_one_subinterval;
      //note: we must have this to be exactly an integer.

      printf("[rank %d] nb_proc_shared_one_subinterval=%d\n", rank, nb_proc_shared_one_subinterval );
      printf("[rank %d] number_in_subinterval=%d\n", rank, number_in_subinterval );
      fflush(stdout);


      for (i=0; i<(Nnz*m); i++){
        indices[i]=point_data[i+Nnz*number_in_subinterval*t_Interval_length_subinterval_loc];
        wghts[i] = weights[i+ Nnz*number_in_subinterval*t_Interval_length_subinterval_loc];
      }
      for(i=0; i<m; i++){
        b[i] = signal[i+number_in_subinterval*t_Interval_length_subinterval_loc];
      }

      free(point_data);
      free(signal);
      free(weights);


  }
  else {

    //for the case we share the stationary intervals in severals processes with no big data flag
    //NB: this case is not up-to-date with the latest data format and should never happen
    if (nb_proc_shared_one_interval>1) {
      part_id = rank/nb_proc_shared_one_interval; //floor(rank/nb_proc_shared_one_interval);
      number_in_interval = rank%nb_proc_shared_one_interval;
      printf("[rank %d] part_id=%d\n", rank, part_id );  //interval id number
      printf("[rank %d] number_in_interval=%d\n", rank, number_in_interval );

      point_data = (int *) malloc(Nnz*t_Interval_length * sizeof(int));
      signal     = (double *) malloc(t_Interval_length * sizeof(double));

      //read the all stationary interval
      for (i=0; i < t_Interval_loop; ++i) {
        ioReadfile(t_Interval_length, part_id, point_data+t_Interval_length_true*Nnz*i, signal+t_Interval_length_true*i);
      }
      //just keep the relevant part of the stationary interval for the local process
      for (i=0; i<(Nnz*m); i++)
      indices[i]=point_data[i+number_in_interval*t_Interval_length_loc];

      for(i=0; i<(m); i++)
      b[i] = signal[i+number_in_interval*t_Interval_length_loc];

      free(point_data);
      free(signal);

    }
    else { //for the case we dont need to share
      //Read the relevants raw inputs data from files distributed by stationary period
      //note: Work only for no sharing stationnary interval. 1 proc for 1 or more stationary intervals

      for (k=0; k < Nb_t_Intervals_loc; ++k) {
        point_data = indices + t_Interval_length*Nnz*k;
        signal = b + t_Interval_length*k;
        weights = wghts + t_Interval_length*Nnz*k;
        part_id = Nb_t_Intervals_loc*rank + k;
        int jump=0;
        for (i=0; i < t_Interval_loop; ++i) {
          ioReadTOAST_data(jump, i, t_Interval_length_true, part_id, point_data+t_Interval_length_true*Nnz*i, signal+t_Interval_length_true*i, weights+t_Interval_length_true*Nnz*i);
        }
      }//end of the loop over the intervals
    }//end if
  }//End if

  /*======== Construct the operator ========*/

  Mat	A;			        //pointing matrix structure
//Pointing matrix init
  A.trash_pix =0;
  MatInit( &A, m, Nnz, indices, wghts, pointing_commflag, MPI_COMM_WORLD);

  t=MPI_Wtime();

  if (rank==0) {
    printf("[rank %d] Reading input data and initializing pointing matrix time=%lf \n", rank, t-st);
  }
  fflush(stdout);
  // printf("A.lcount = %d\n", A.lcount);
  st=MPI_Wtime();
//Build pixel-to-time domain mapping
  id0pix = (int *) malloc(A.lcount/(A.nnz) * sizeof(int)); //index of the first time sample pointing to each pixel
  ll = (int *) malloc(m * sizeof(int)); //linked list of time samples indexes

  //initialize the mapping arrays to -1
  for(i=0; i<m; i++){
    ll[i] = -1;
  }
  for(j=0; j<A.lcount/(A.nnz); j++){
    id0pix[j] = -1;
  }
  //build the linked list chain of time samples corresponding to each pixel
  // for(i=0; i<m; i++){
  //   if(id0pix[A.indices[i*A.nnz]/(A.nnz)] == -1)
  //     id0pix[A.indices[i*A.nnz]/(A.nnz)] = i;
  //   else
  //     set_link(ll, id0pix[A.indices[i*A.nnz]/(A.nnz)], i);
  // }
  for(i=0;i<m;i++){
    if(id0pix[A.indices[i*A.nnz]/(A.nnz)] == -1)
      id0pix[A.indices[i*A.nnz]/(A.nnz)] = i;
    else{
      ll[i] = id0pix[A.indices[i*A.nnz]/(A.nnz)];
      id0pix[A.indices[i*A.nnz]/(A.nnz)] = i;
    }
  }

  A.id0pix = id0pix;
  A.ll = ll;
  t=MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==0) {
    printf("[rank %d] Total pixel-to-time domain mapping time=%lf \n", rank, t-st);
  }
  fflush(stdout);

  m=A.m;					//number of local time samples
  n=A.lcount;					//number of local pixels


  // PCG starting vector input definition for the pixel domain map (MatInit gives A.lcount)
  int *lhits;
  double *cond;

  x   = (double *) malloc(n*sizeof(double));
  cond = (double *) malloc((int)(n/3)*sizeof(double));
  lhits = (int *) malloc((int)(n/3) * sizeof(int));
  for(j=0; j<n; j++){
    x[j] = 0.;
    if(j%3 == 0){
      lhits[(int)(j/3)] = 0;
      cond[(int)(j/3)] = 0.;
    }
  }

//Create piecewise Toeplitz matrix

//specifics parameters:
  int nb_blocks_tot = Nb_t_Intervals;
  int n_block_avg = M/nb_blocks_tot;  //should be equal to t_Intervals_length in this example
                                      //because we dont have flotting blocks
  int lambda_block_avg = LambdaBlock;

//flags for Toeplitz product strategy
  Flag flag_stgy;
  flag_stgy_init_auto(&flag_stgy);

//to print something on screen
  flag_stgy.flag_verbose=1;

//to define fixed bs:
//  flag_stgy.flag_bs = 1;
//  flag_stgy.fixed_bs = pow(2,17);

//define Toeplitz blocks list and structure for Nm1
  Block *tpltzblocks;
  Tpltz Nm1;

//dependants parameters:
  int64_t nrow = M;
  int mcol = 1;

  int64_t id0 = gif;
  int local_V_size = m;

  int Tsize = lambda_block_avg;
  double *T;  //toeplitz data storage
  T  = (double *) calloc(Tsize ,sizeof(double));

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank==0){
    ioReadTpltzfile( Tsize, T);
    printf("Tsize = %d",Tsize);
    printf("\n correlation = [%f,...,%f]",T[0],T[Tsize-1]);
  }
  MPI_Bcast(T, Tsize,  MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);


  int nb_blocks_loc;
  nb_blocks_loc = ceil( local_V_size*1.0/n_block_avg );

  double nb_blocks_loc_part =  (local_V_size*1.0)/(n_block_avg) ;

// check special cases to have exact number of local blocks
  if ((id0/n_block_avg + nb_blocks_loc) * n_block_avg < (id0+local_V_size))
    nb_blocks_loc=nb_blocks_loc+1;

  if (rank==0 | rank==1) {
    printf("M=%ld, m=%d \n", M, m);
    printf("gif = %ld \n", gif);
    fflush(stdout);
  }

  int nb_proc_shared_a_block = ceil( size*1.0/nb_blocks_tot );
  int nb_comm = (nb_proc_shared_a_block)-1 ;

//Block definition
  tpltzblocks = (Block *) malloc(nb_blocks_loc * sizeof(Block));
  defineBlocks_avg(tpltzblocks, T, nb_blocks_loc, n_block_avg, lambda_block_avg, id0 );
  defineTpltz_avg( &Nm1, nrow, 1, mcol, tpltzblocks, nb_blocks_loc, nb_blocks_tot, id0, local_V_size, flag_stgy, MPI_COMM_WORLD);

//print Toeplitz parameters for information
  if (rank==0 | rank==1) {
    printf("[rank %d] size=%d, nrow=%ld, local_V_size=%ld, id0=%ld \n", rank, size, nrow, local_V_size, id0);
    printf("[rank %d] nb_blocks_tot=%d, nb_blocks_loc=%d, n_block_avg=%d, lambda_block_avg=%d \n", rank, nb_blocks_tot, nb_blocks_loc, n_block_avg, lambda_block_avg);
    printf("[rank %d] nb_proc_shared_a_block=%d, nb_comm=%d \n", rank, nb_proc_shared_a_block, nb_comm);
    fflush(stdout);
  }

  /*======== Construct the preconditioner ========*/

  Mat BJ;
  precondblockjacobilike(&A, Nm1, &BJ, b, cond, lhits);

  n=A.lcount-(A.nnz)*(A.trash_pix);
// Reallocate memory for well-conditioned map
  tmp = realloc(x, n * sizeof(double));
  if(tmp !=NULL){
    x = tmp;
  }

  N += n;
  MPI_Allreduce(MPI_IN_PLACE,&N,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (rank == 0) {
    printf("GLOBAL N = %ld \n", N);
  }

  double *pixpond;
  pixpond = (double *) malloc(n*sizeof(double));

  get_pixshare_pond( &A, pixpond);
  /*============= Construct the rhs =============*/
  double *rhs;
  rhs = (double *) malloc(n*sizeof(double));
  get_rhs(A, Nm1, b, x, rhs);

  double normres_init = 0.0;
  for(i = 0; i < n; i++) {
    normres_init += rhs[i] * rhs[i] * pixpond[i];
  }
  //printf("rank = %d, ||r_0||^2_loc = %.6f \n",rank, normres_init);

  MPI_Allreduce(MPI_IN_PLACE,&normres_init,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  normres_init = sqrt(normres_init);
  if (rank == 0) {
    printf("GLOBAL ||r_0|| = %.6f \n", normres_init);
    fflush(stdout);
  }

  /*================ ECG solve ================*/
  preAlps_ECG_t ecg;
  // Set parameters
  ecg.comm = MPI_COMM_WORLD;
  ecg.globPbSize = N;
  ecg.locPbSize = n;
  ecg.maxIter = maxIter;
  ecg.enlFac = enlFac;
  ecg.tol = tol;
  ecg.ortho_alg = (ortho_alg == 0 ? ORTHODIR : ORTHOMIN);
  ecg.bs_red = (bs_red == 0 ? NO_BS_RED : ADAPT_BS);
  int rci_request = 0;
  int stop = 0;
  double* sol = NULL;       // solution = one vector
  sol = (double*) malloc(n*sizeof(double));
  double* rel_res = NULL;   // relative residual
  rel_res = (double*) malloc((maxIter)*sizeof(double));
  // timings
  double time_ECG_total = 0;
  double time_AV = 0;
  double time_invMV = 0;

// CPLM_TIC(step1,"ECGSolve")

  time_ECG_total = MPI_Wtime();
  // Allocate memory and initialize variables
  preAlps_oECGInitialize(&ecg, rhs, &rci_request, A.lindices+(A.nnz)*(A.trash_pix));
  ecg.normb = normres_init;
  rel_res[0] = 1.0;

  // Finish initialization
  time_invMV += Opmmpreconditioner(A, BJ, ecg.R_p, ecg.P_p, ecg.bs);
  // preconditioner ecg.R -> ecg.P
  time_AV += Opmmmatrix(A, Nm1, ecg.P_p, ecg.AP_p, ecg.bs);
  // block operator ecg.P -> ecg.AP

  // Main loop
  while (stop != 1) {
    preAlps_oECGIterate(&ecg,&rci_request, pixpond);
    if (rci_request == 0) {
      time_AV += Opmmmatrix(A, Nm1, ecg.P_p, ecg.AP_p, ecg.bs);
      // block operator ecg.P -> ecg.AP
    }
    else if (rci_request == 1) {

      preAlps_oECGStoppingCriterion(&ecg,&stop, pixpond);

      rel_res[ecg.iter] = ecg.res/ecg.normb;
      if (stop == 1) break;
      if (ecg.ortho_alg == ORTHOMIN) {
        time_invMV += Opmmpreconditioner(A, BJ, ecg.R_p, ecg.Z_p, ecg.enlFac);
        // preconditioner ecg.R -> ecg.Z
      }
      else if (ecg.ortho_alg == ORTHODIR) {
        time_invMV += Opmmpreconditioner(A, BJ, ecg.AP_p, ecg.Z_p, ecg.bs);
        // preconditioner ecg.AP -> ecg.Z
      }
    }
  }
  time_ECG_total = MPI_Wtime() - time_ECG_total;

  if (rank == 0) {
    preAlps_oECGPrint(&ecg,5);

    printf("*** TIMING, total ECG time = %e s\n", time_ECG_total);
    printf("*** TIMING,     A x V time = %e s\n", time_AV);
    printf("*** TIMING,  invM x V time = %e s\n", time_invMV);

    char filename[40];
    sprintf(filename, "rel_res_enlfac_%d_o=%d.txt", ecg.enlFac, ortho_alg);
    FILE* fp = fopen(filename,"w");
    for(i = 0; i <= ecg.iter;i++){
       fprintf (fp, "%.15e\n",rel_res[i]);
    }
    fclose(fp);
  }
  free(rel_res);

  // Retrieve solution and free memory
  preAlps_oECGFinalize(&ecg,sol);
// CPLM_TAC(step1)

  /*================ Finalize ================*/

  // Free arrays
  if (rhs != NULL) free(rhs);
  if (sol != NULL) free(sol);
// CPLM_CLOSE_TIMER

  // CPLM_printTimer(NULL);
  MPI_Finalize();
  return 0;
}
/******************************************************************************/

int partition(int64_t *gif, int *m, int64_t M, int rank, int size){
  int64_t r, k;
  k = M / size;
  r = M - k*size;
  if( rank < r){
    *gif = (k+1) * rank;
    *m = k+1;
  }
  else{
    *gif = r*(k+1) + k*(rank-r);
    *m = k;
    }
  return 0;
}

int get_rhs(Mat A, Tpltz Nm1, double *b, double *x, double *rhs)
{ // rhs = A^T*Nm1* (b - A*x_0 )
  int i; //some indexes
  int m, n;
  double *_g;   // time domain vector

  m = A.m;    //number of local time samples
  n = A.lcount-(A.nnz)*(A.trash_pix);   //number of local pixels

  _g = (double *) malloc(m*sizeof(double));

  MatVecProd(&A, x, _g, 0);
  for(i=0; i<m; i++)	//
    _g[i] = b[i] - _g[i];		//
  stbmmProd(Nm1,_g);
  TrMatVecProd(&A, _g, rhs, 0);

  return 0;
}

double Opmmmatrix(Mat A, Tpltz Nm1, double *X, double *Y, int ncol)
{ // Y = A^T*Nm1*A * X
  double timing = MPI_Wtime();

  int i, j; //some indexes
  int m, n;
  double *_g;   // time domain vector
  double *x, *g; // map domain vector

  m = A.m;    //number of local time samples
  n = A.lcount-(A.nnz)*(A.trash_pix);   //number of local pixels

  _g = (double *) malloc(m * sizeof(double));
  g = (double *) malloc(n * sizeof(double));

  for(i=0; i<ncol; i++){
    // get column vector x
    x = X + i*n;
    MatVecProd(&A, x, _g, 0);
    stbmmProd(Nm1,_g);
    TrMatVecProd(&A, _g, g, 0);
    for(j=0; j<n; j++){
      Y[i*n+j] = g[j];
    }
  }
  free(_g);
  free(g);
  return MPI_Wtime() - timing;
}


double Opmmpreconditioner(Mat A, Mat BJ,  double *X, double *Y, int ncol)
{ // Y = M^-1 X
  double timing = MPI_Wtime();

  int i, j; //some indexes
  int n;
  double *x, *Cg; // map domain vector

  n = A.lcount-(A.nnz)*(A.trash_pix);   //number of local pixels

  Cg = (double *) malloc(n * sizeof(double));

  for(i=0; i<ncol; i++){
    // get column vector x
    x = X + i*n;
    MatVecProd(&BJ, x, Cg, 0);
    for(j=0; j<n; j++){
      Y[i*n+j] = Cg[j];
    }
  }
  free(Cg);
  return MPI_Wtime() - timing;
}
