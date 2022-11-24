/** @file mapmatc.h
    @author Pierre Cargemel
    @date October 2012 */

#ifndef MAPMATC_H
#define MAPMATC_H

#include "butterfly.h"
#include "alm.h"
#include "als.h"
#include "bitop.h"
#include "cindex.h"
#include "csort.h"
#include "ring.h"



#define NONE 0
#define RING 1
#define BUTTERFLY 2
#define NONBLOCKING 3
#define NOEMPTY 4
#define SEQ 0
#define OMP 1
#define GPU 2



/** @brief Matrix structure
    @n A* = (A0* | A1* | ... | Ap-1* )
    */
typedef struct {
  int		flag;			// flag for communication scheme (NONE, RING, BUTTERFLY ...)
  int		r;			// number of local coarse space
  int		*m;			// table containing number of rows in each caorse space
  int		*nnz;                   // number non-zero per row in each corse space
  int		*disp;                  // displacement
  int		**indices;		// column rows indices tab;
  double	**values;		// non-zero values tab;
  //--------local shaping---------------
  int		lcount;
  int		*lindices;		// local indices tab (monotony with global numbering); size = lcount
#ifdef W_MPI
  MPI_Comm	comm;                   // MPI communicator
  //--------com shaping-----------------
  int		*com_indices, com_count;// communicated indices tab, and size
  int		steps;			// number of steps in the communication scheme
  int		*nS, *nR;		// number of indices (to send and to receive); size = steps
  int		**R, **S;		// sending or receiving indices tab
#endif
}CMat;


int CMatInit(CMat *A, int r, int *m, int *nnz, int **indices, double **values, int flag
#ifdef W_MPI
  ,MPI_Comm comm
#endif
);

int CMatFree(CMat *A);

#ifdef W_MPI
int CMatComShape(CMat *A, int flag);
#endif

int CMatVecProd(CMat *A, double *x, double *y, int pflag);

int CTrMatVecProd(CMat *A, double *y, double* x, int pflag);

#endif /* MAPMATC_H */