/** @file mappraiser.h
    @brief <b> Declaration of the backbone routines of the map-making code.</b>
    @author Hamza El Bouhargani
    @date May 2019 */

#ifndef MAPPRAISER_H
#define MAPPRAISER_H

#include "midapack.h"

/* Define the banded block Toeplitz matrix */

int defineTpltz_avg(Tpltz *Nm1, int64_t nrow, int m_cw, int m_rw, Block *tpltzblocks, int nb_blocks_loc, int nb_blocks_tot, int64_t idp, int local_V_size, Flag flag_stgy, MPI_Comm comm);

int defineBlocks_avg(Block *tpltzblocks, double *T, int nb_blocks_loc, void *local_blocks_sizes, int lambda_block_avg, int64_t id0);

/* IO routines */

int ioWritebinfile(int mapsize, int mappart_id, int *lstid, double *map, double *cond, int *lhits);

int ioReadbinfile(int mapsize, int mappart_id, int *lstid, double *map, double *cond, int *lhits);

// Import HEALPix routines with little tweaks to pass datatypes in arguments
static void util_fail_(const char *file, int line, const char *func,
                       const char *msg);

#if defined(__GNUC__)
#define UTIL_FUNC_NAME__ __func__
#else
#define UTIL_FUNC_NAME__ "unknown"
#endif
#define UTIL_ASSERT(cond, msg) \
    if (!(cond))               \
    util_fail_(__FILE__, __LINE__, UTIL_FUNC_NAME__, msg)
#define UTIL_FAIL(msg) \
    util_fail_(__FILE__, __LINE__, UTIL_FUNC_NAME__, msg)

static void setCoordSysHP(char coordsys, char *coordsys9);

static void printerror(int status);

void write_map(void *signal, int type, long nside, const char *filename,
               char nest, const char *coordsys);

/* Preconditioner routines */

struct Precond
{
    int precond; // 0 = BJ, 1 = 2lvl a priori, 2 = 2lvl a posteriori
    int n;
    int Zn;
    Mat BJ_inv;
    Mat BJ;
    double *pixpond;

    /* 2 lvl only (NULL otherwise) */
    double **Z;
    double **AZ;
    double *Em1; // size Zn*Zn
    double *Qg;  // size n
    double *AQg; // size n
    double *Qtx; // size Zn
    double *w;   // size Zn
};
typedef struct Precond Precond;

// Block-Jacobi preconditioner
int precondblockjacobilike(Mat *A, Tpltz *Nm1, Mat *BJ_inv, Mat *BJ, double *b, double *cond, int *lhits);

// Preconditioner constructor
void build_precond(struct Precond **out_p, double **out_pixpond, int *out_n, Mat *A, Tpltz *Nm1, double **in_out_x, double *b, const double *noise, double *cond, int *lhits, double tol, int Zn, int precond);

// Product of the preconditioner with a map vector
void apply_precond(struct Precond *p, const Mat *A, Tpltz *Nm1, double *g, double *Cg);

// Free memory of the preconditioner
void free_precond(struct Precond **in_out_p);

/* PCG routines */

// Pixel share ponderation to deal with overlapping pixels between multiple MPI procs
void get_pixshare_pond(Mat *A, double *pixpond);

// PCG routine
int PCG_GLS_true(char *outpath, char *ref, Mat *A, Tpltz Nm1, PCG_var *PCG_variable, double *b, double *noise, double *cond, int *lhits, double tol, int K, int precond, int Z_2lvl, S2HAT_parameters *S2HAT_params);

// ECG routine
#ifdef WITH_ECG
int ECG_GLS(char *outpath, char *ref, Mat *A, Tpltz *Nm1, double *x, double *b, double *noise, double *cond, int *lhits, double tol, int maxIter, int enlFac, int ortho_alg, int bs_red);
#endif

#endif /* MAPPRAISER_H */
