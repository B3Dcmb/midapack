/** @file mapmat.h
    @brief <b> Declarations of the matrix type and his associated routines.</b>
    @n these routines are developed to handle sparse matrices.
    Typically, in the CMB Data Analysis context, it is especially developed
   handle pointing or unpointing matrices. Thus, the unpointing matrix @a A can
   be defined as a MIDAS_Mat. Operating with the pointing matrices can be done
   without redefining a new matrix.
    @author Pierre Cargemel
    @date November 2011 */

/* Update by Hamza El Bouhargani
    @date February 2019 */

#ifndef MAPMAT_H
#define MAPMAT_H

#include "alm.h"
#include "als.h"
#include "bitop.h"
#include "butterfly.h"
#include "cindex.h"
#include "csort.h"
#include "ring.h"
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#define NONE        0
#define RING        1
#define BUTTERFLY   2
#define NONBLOCKING 3
#define NOEMPTY     4
#define ALLTOALLV   5
#define ALLREDUCE   6
//================Modification introduced by Sebastien Cayrols : 01/09/2015 ;
// Berkeley
#define BUTTERFLY_BLOCKING_1 7
#define BUTTERFLY_BLOCKING_2 8
#define NOEMPTYSTEPRING      9
//================End modification
#define SEQ 0
#define OMP 1
#define GPU 2

/** @brief Matrix structure
    @n A* = (A0* | A1* | ... | Ap-1* )
    */
typedef struct {
    int  flag;      // flag for communication scheme (NONE, RING, BUTTERFLY ...)
    int  m;         // number local rows
    int  nnz;       // number non-zero per rows
    int  trash_pix; // flag for presence of trash pixel
    int *indices;   // column indices tab; size = m * nnz; can be a global or
                  // local numbering
    double *values; // non-zero values tab; size = m * nnz
    int    *id_last_pix; // index of the last time sample pointing to each pixel
                      // (no nnz repeat factor)
    int *ll;             // linked list of time samples indexes linked by pixels
    //--------local shaping---------------
    int  lcount;
    int *lindices; // local indices tab (monotony with global numbering); size =
                   // lcount
#ifdef W_MPI
    MPI_Comm comm; // MPI communicator
    //--------com shaping-----------------
    int  *com_indices, com_count; // communicated indices tab, and size
    int   steps;                  // number of steps in the communication scheme
    int  *nS, *nR; // number of indices (to send and to receive); size = steps
    int **R, **S;  // sending or receiving indices tab
#endif
} Mat;

int MatInit(Mat *A, int m, int nnz, int *indices, double *values, int flag
#ifdef W_MPI
            ,
            MPI_Comm comm
#endif
);

void MatSetIndices(Mat *A, int m, int nnz, int *indices);

void MatSetValues(Mat *A, int m, int nnz, double *values);

void MatFree(Mat *A);

int MatLocalShape(Mat *A, int sflag);

#if W_MPI

int MatComShape(Mat *A, int flag, MPI_Comm comm);

#endif

int MatVecProd(Mat *A, double *x, double *y, int pflag);

int TrMatVecProd(Mat *A, double *y, double *x, int pflag);

#if W_MPI

int TrMatVecProd_Naive(Mat *A, double *y, double *x, int pflag);

#endif

int MatLoad(Mat *A, char *filename);

int MatSave(Mat *A, char *filename);

#if W_MPI

int MatInfo(Mat *A, int master, char *filename);

#endif

int greedyreduce(Mat *A, double *x);

#ifdef __cplusplus
}
#endif

#endif // MAPMAT_H

// Doxygen definitions

/** @defgroup matmap Pointing module
 *  Pointing operations module
 */

/** @defgroup matmap_group1 user interface (API)
 *  These are routines "officially" accessible by a user.
 *  @ingroup matmap
 */

/** @defgroup matmap_group11 utility routines
 *	These are auxiliary utility routines. They may be used in both
 *  memory-shared/sequential and dsitributed contexts, though their syntax
 *	may change depending on that.
 *
 *  @ingroup matmap_group1
 */

/** @defgroup matmap_group12 main routines
 *	These are core routines of the pointing library.
 *  They come in two flavors:
 *  - shared-memory: multithreaded (openMP/sequential) routines
 *  - distributed-memory (MPI) routines
 *
 *  @ingroup matmap_group1
 */

/** @defgroup matmap_group12a multithreaded/sequential routines
 *  These are sequential and/or shared-memory routines.
 *  @ingroup matmap_group12
 */

/** @defgroup matmap_group12b distributed memory (MPI) routines
 *  These are distributed-memory routines.
 *
 *  <b>Note</b> that if the MPI flag is not set
 *  on the compilation stage, most of the routines listed here will run
 *  sequentially unless explicitly stated to the contrary in their description.
 *
 *  @ingroup matmap_group12
 */

/** @defgroup matmap_group2 internal routines
 *  These are auxiliary, internal routines, not intended to be used by no-expert
 * user. They are divided in two groups:
 *  - low level computation routines
 *  - internal routines
 *  @ingroup matmap
 */

/** @defgroup matmap_group21 low-level computation routines
 *  These are low-level computation routines not really expected to be used by
 * users, e.g., the syntax may evolve unexpectedly, but may be be found on
 * occasions useful ...
 *  @ingroup matmap_group2
 */

/** @defgroup matmap_group22 lower internal routines
 *  These are low level internal routines. These are generally not to be used by
 * external users.
 *  @ingroup matmap_group2
 */
