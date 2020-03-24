/** @file mappraiser.h
    @brief <b> Declaration of the backbone routines of the map-making code.</b>
    @author Hamza El Bouhargani
    @date May 2019 */

/* Define the banded block Toeplitz matrix */
int defineTpltz_avg( Tpltz *Nm1, int64_t nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks_loc, int nb_blocks_tot, int64_t idp, int local_V_size, Flag flag_stgy, MPI_Comm comm);

int defineBlocks_avg(Block *tpltzblocks, double *T, int nb_blocks_loc, void *local_blocks_sizes, int lambda_block_avg, int64_t id0 );

/* IO routines */
int ioWritebinfile( int mapsize, int mappart_id, int *lstid, double *map, double *cond, int *lhits);

int ioReadbinfile( int mapsize, int mappart_id, int *lstid, double *map, double *cond, int *lhits);

//Import HEALPix routines with little tweaks to pass datatypes in arguments
static void util_fail_ (const char *file, int line, const char *func,
  const char *msg);

#if defined (__GNUC__)
#define UTIL_FUNC_NAME__ __func__
#else
#define UTIL_FUNC_NAME__ "unknown"
#endif
#define UTIL_ASSERT(cond,msg) \
  if(!(cond)) util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)
#define UTIL_FAIL(msg) \
  util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)

static void setCoordSysHP(char coordsys,char *coordsys9);

static void printerror (int status);

void write_map (void *signal, int type, long nside, const char *filename,
  char nest, const char *coordsys);

/* Preconditioner routines */

// Block-Jacobi preconditioner
int precondblockjacobilike(Mat *A, Tpltz Nm1, Mat *BJ, double *b, double *cond, int *lhits);

// Point Jacobi preconditioner
int precondjacobilike(Mat A, Tpltz Nm1, int *lhits, double *cond, double *vpixDiag);

// Local product A^T * diagNM1 * A
int getlocalW(Mat *A, Tpltz Nm1, double *vpixBlock, int *lhits);

int getlocDiagN(Mat *A, Tpltz Nm1, double *vpixDiag);

int DiagAtA(Mat *A, double *diag, int pflag);

// Communication routine for building the pixel blocks of the preconditioner
int commScheme(Mat *A, double *vpixDiag, int pflag);

/* PCG routines */

// Pixel share ponderation to deal with overlapping pixels between multiple MPI procs
int get_pixshare_pond(Mat *A, double *pixpond);

//PCG routine
int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz Nm1, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int K);

int PCG_GLS_templates(char *outpath, char *ref, Mat *A, Tpltz Nm1, TemplateClass *X, double *B, int **sweeptstamps, int npoly, int *nsweeps, int **az_binned, int n_sss_bins, int nces, int nb_blocks_loc, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int K, double sampling_freq);

//ECG routine
int ECG_GLS(char *outpath, char *ref, Mat *A, Tpltz Nm1, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int maxIter, int enlFac, int ortho_alg, int bs_red);
