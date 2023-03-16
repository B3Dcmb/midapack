/** @file cindex.h 
    @author Pierre Cargemel
    @date May 2011 */

#ifndef MAPMAT_CINDEX_H
#define MAPMAT_CINDEX_H

int sindex(int *array1, int nb1, int *array2, int nb2);

#if OPENMP
int omp_pindex(int *array1, int nb1, int *array2, int nb2);
#endif

#endif //MAPMAT_CINDEX_H
