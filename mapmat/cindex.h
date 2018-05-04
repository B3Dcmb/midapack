/** @file cindex.h 
    @author Pierre Cargemel
    @date May 2011 */



int sindex(int *array1, int nb1, int *array2, int nb2);

#if OPENMP
int omp_pindex(int *array1, int nb1, int *array2, int nb2);
#endif
