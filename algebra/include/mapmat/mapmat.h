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

#ifdef W_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <bits/stdint-uintn.h>
#include <stdbool.h>

/** @brief Matrix structure
    @n A* = (A0* | A1* | ... | Ap-1* )
    */
typedef struct mat_t {
    int m;             // number of local rows
    int nnz;           // number of non-zero values per row
    int trash_pix;     // number of extra pixels
    bool ignore_extra; // if true, do not include extra pixels in estimated map
    int *indices;      // column indices tab; size = m; can be a global or
                       // local numbering
    double *values;    // non-zero values tab; size = m * nnz
    uint8_t *flags;    // array indicating flagged samples (0=valid, 1=flagged)
    int *id_last_pix;  // index of the last sample pointing to each pixel
    int *ll;           // linked list of time samples indexes linked by pixels
    //--------local shaping---------------
    int lcount;    // number of local pixels (including extra)
    int *lindices; // local indices tab (monotony with global numbering)
#ifdef W_MPI
    MPI_Comm comm;               // MPI communicator
    int *com_indices, com_count; // communicated indices tab, and size
#endif
} Mat;

void MatInit(Mat *A, int m, int nnz, int *indices, double *values,
             uint8_t *flags
#ifdef W_MPI
             ,
             MPI_Comm comm
#endif
);
void MatSetIndices(Mat *A, int m, int *indices);
void MatSetValues(Mat *A, int m, int nnz, double *values);
void MatFree(Mat *A);
void MatLocalShape(Mat *A, int sflag);

// communication between processes
#if W_MPI
void MatComShape(Mat *A, MPI_Comm comm);
void greedyreduce(Mat *A, double *x);
#endif

// matrix-vector operations
void MatVecProd(Mat *A, const double *x, double *y);
void TrMatVecProd(Mat *A, const double *y, double *x);

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
