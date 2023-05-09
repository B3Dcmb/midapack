/**
 * \file    ecg.h
 * \author  Olivier Tissot
 * \date    2016/06/24
 * \brief   Enlarged Preconditioned C(onjugate) G(radient) solver
 *
 * \details Implements Orthomin, Orthodir as well as their dynamic
 *          counterparts (BF-Omin and D-Odir).
 */

/******************************************************************************/
/*                                  INCLUDE                                   */
/******************************************************************************/
#ifndef OVERLAP_ECG_H
#define OVERLAP_ECG_H

#include <operator.h>
#include <cplm_v0_timing.h>

/**
 * \enum preAlps_ECG_Ortho_Alg_t
 * \brief A-orthonormalization algorithm
 * \author Olivier Tissot
 * \date 2016/06/24
 */
typedef enum {
  ORTHOMIN,
  ORTHODIR,
  ORTHODIR_FUSED
} preAlps_ECG_Ortho_Alg_t;
/**
 * \enum preAlps_ECG_Block_Size_Red_t
 * \brief Block size reduction
 * \author Olivier Tissot
 * \date 2016/06/24
 */
typedef enum {
  ADAPT_BS,
  NO_BS_RED
} preAlps_ECG_Block_Size_Red_t;

/**
* \struct preAlps_ECG_t
* \brief Enlarged Conjugate Gradient solver
* \author Olivier Tissot
* \date 2016/06/24
*/
typedef struct {
  /* Input variable */
  double* b;                    /**< Right hand side */

  /* Internal symbolic variables */
  CPLM_Mat_Dense_t* X;     /**< Approximated solution */
  CPLM_Mat_Dense_t* R;     /**< Residual */
  CPLM_Mat_Dense_t* V;     /**< Descent directions ([P,P_prev] or P) */
  CPLM_Mat_Dense_t* AV;    /**< A*V */
  CPLM_Mat_Dense_t* Z;     /**< Preconditioned residual (Omin) or AP (Odir) */
  CPLM_Mat_Dense_t* alpha; /**< Descent step */
  CPLM_Mat_Dense_t* beta;  /**< Step to construt search directions */

  /** User interface variables */
  CPLM_Mat_Dense_t* P;      /**< Search directions */
  CPLM_Mat_Dense_t* AP;     /**< A*P */
  double* R_p;              /**< Residual */
  double* P_p;              /**< Search directions */
  double* AP_p;             /**< A*P_p */
  double* Z_p;              /**< Preconditioned residual (Omin) or AP (Odir) */

  /** Working arrays */
  double*           work;
  int*              iwork;

  /* Single value variables */
  double            normb;     /**< norm_2(b) */
  double            res;       /**< norm_2 of the residual */
  int               iter;      /**< Iteration */
  int               bs;        /**< Block size */
  int               kbs;       /**< Krylov basis size */

  /* Options and parameters */
  int                          globPbSize; /**< Size of the global problem */
  int                          locPbSize;  /**< Size of the local problem */
  int                          maxIter;    /**< Maximum number of iterations */
  int                          enlFac;     /**< Enlarging factor */
  double                       tol;        /**< Tolerance */
  preAlps_ECG_Ortho_Alg_t      ortho_alg;  /**< A-orthonormalization algorithm */
  preAlps_ECG_Block_Size_Red_t bs_red;     /**< Block size reduction */
  MPI_Comm             comm;               /**< MPI communicator */

  /* Timings */
  double tot_t ;  /* Total */
  double comm_t;  /* Communication */
  double trsm_t;  /* trsm  */
  double gemm_t;  /* gemm  */
  double potrf_t; /* potrf */
  double pstrf_t; /* pstrf */
  double lapmt_t; /* lapmt */
  double gesvd_t; /* gesvd */
  double geqrf_t; /* geqrf */
  double ormqr_t; /* ormqr */
  double copy_t;  /* copy */

} preAlps_ECG_t;
/******************************************************************************/

/******************************************************************************/
/*                                    CODE                                    */
/******************************************************************************/

/**
 * \brief Create the solver and allocate memory
 *
 * \param[in, out] ecg solver structure
 * \param[in]      rhs the local part of the right-hand side
 * \param[out]     rci_request the initialized RCI flag
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  preAlps_oECGInitialize(preAlps_ECG_t* ecg, double* rhs, int* rci_request, int* glob_indices);

/**
 * \brief Performs different steps in ECG iteration according to the value of
 * rci_request
 *
 * \param[in, out] ecg solver structure
 * \param[in, out] rci_request the RCI flag
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  preAlps_oECGIterate(preAlps_ECG_t* ecg, int* rci_request, double* w);

/**
 * \brief Check for the residual norm and return a boolean that is true if the
 * normalized residual is lower than the specified tolerance
 *
 * \param[in, out] ecg solver structure
 * \param[in, out] rci_request the RCI flag
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  preAlps_oECGStoppingCriterion(preAlps_ECG_t* ecg, int* stop, double* w);

/**
 * \brief Releases the internal memory and returns the solution
 *
 * \param[in, out] ecg solver structure
 * \param[out]     solution the local part of the solution
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  preAlps_oECGFinalize(preAlps_ECG_t* ecg, double* solution);

/**
 * \brief Print informations on the solver
 *
 * \param[in] ecg solver structure
 * \param[in] verbosity the level of information printed: if verbosity <= 1
 * then the iteration count, the residual and the block size are printed, if
 * verbosity > 1 then the full detail of the memory footprint of the solver is
 * also printed
 */
void preAlps_oECGPrint(preAlps_ECG_t* ecg, int verbosity);

/* "Private" functions */

/**
 * \brief Private function
 * \details Allocate memory for the solver
 *
 * \param[in, out] ecg solver structure
 */
int  _preAlps_oECGMalloc(preAlps_ECG_t* ecg);

/**
 * \brief Private function
 * \detail Initialize the solver assuming that the memory has been allocated
 *
 * \param[in, out] ecg solver structure
 * \param[in]      rhs the local part of the right-hand side
 * \param[out]     rci_request the initialized RCI flag
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  _preAlps_oECGReset(preAlps_ECG_t* ecg, double* rhs, int* rci_request, int* glob_indices);

/**
 * \brief Private function
 * \detail Returns the solution
 *
 * \param[in, out] ecg solver structure
 * \param[out]     solution the local part of the solution
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  _preAlps_oECGWrapUp(preAlps_ECG_t* ecg, double* solution);

/**
 * \brief Private function
 * \detail Release memory of the solver
 *
 * \param[in, out] ecg solver structure
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
void _preAlps_oECGFree(preAlps_ECG_t* ecg);

/**
 * \brief Private function
 * \detail Enlarge the vector x
 *
 * \param[in] x the local part of the vector to enlarge
 * \param[out] XSplit the local part of the enlarged vector
 * \param[in] colIndex the index of the column of XSplit where the local vector * is put
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  _preAlps_oECGSplit(double* x, CPLM_Mat_Dense_t* XSplit, int* glob_indices);

/**
 * \brief Private function
 * \detail Orthomin iteration
 *
 * \param[in, out] ecg solver structure
 * \param[in, out] rci_request the RCI flag
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  _preAlps_oECGIterateOmin(preAlps_ECG_t* ecg, int* rci_request, double* w);

/**
 * \brief Private function
 * \detail Orthodir iteration
 *
 * \param[in, out] ecg solver structure
 * \param[in, out] rci_request the RCI flag
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  _preAlps_oECGIterateOdir(preAlps_ECG_t* ecg, int* rci_request, double* w);

/**
 * \brief Private function
 * \detail Orthodir fused iteration (1 MPI_Allreduce per iteration)
 *
 * \param[in, out] ecg solver structure
 * \param[in, out] rci_request the RCI flag
 * \return 0 if the execution succeeded
 * \return 1 if the execution failed
 */
int  _preAlps_oECGIterateOdirFused(preAlps_ECG_t* ecg, int* rci_request, double* w);
/******************************************************************************/

int _preAlps_oECGweightedMatMult(CPLM_Mat_Dense_t* X, CPLM_Mat_Dense_t* Y, double* w, CPLM_Mat_Dense_t* C);

#endif

