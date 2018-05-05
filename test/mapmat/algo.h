/** @file algo.h 
    @brief <b> Declaration of subroutines handling sets of indices or sets of values. </b>
    @n Working with sparse matrices involves to uses specific data storage. 
    Thus it requires to handle several sets of indices (rows, columns, local, global).
    That's why, this file gathers couples subroutines which handle those sets of indices, as sorting methods, dichotmomic search, intersectioni, union, and so on ...
    These subroutines are called inside the upper routines declared in mapmat.h or partitioning.h. 
    @author Pierre Cargemel
    @date November 2011 */


/** @brief basic sort

    Sorts array
    @param count number of elements to sort
    @param array pointer to element array*/
void sort(int count, int *array);

/** @brief count

    Count number of distinct elements assuming array is in ascending order
    @param count number of elements to sort
    @param array pointer to element array*/
int count(int count, int *array);

/** @brief

    Count number of identical element between 2 arrays assuming arrays are strictly ascending ordered
    @param count1 number of elements in first array (local array)
    @param array1 pointer to element first array
    @param count2 number of elements in second array (parameter array)
    @param array2 pointer to element second array
    @param array3 indices array of the intersection (according local arrayindices)*/
void intersection(int count1, int *array1, int count2, int *array2, int count3, int *array3);

/** @brief

    Count number of identical element between 2 arrays assuming arrays are strictly ascending ordered
    @param count1 number of elements in first array
    @param array1 pointer to element first array
    @param count2 number of elements in second array
    @param array2 pointer to element second array*/
int shared(int count1, int *array1, int count2, int *array2);


/** @brief parses 2 indices/values arrays and sum values into the first one whenever indices are equal. 

    @param count1 number of elements in first array
    @param indices_array1 pointer to element first indices array
    @param values_array1 pointer to element first values array
    @param count2 number of elements in second array
    @param indices_array2 pointer to element second indices array
    @param values_array2 pointer to element second values array*/
void reduce(int count1, int *indices_array1, double *values_array1, int count2, int* indices_array2, double* values_array2);


/** @brief

    Copy an array with possible element overlaps, into a compressed one (no overlaps). 
    Number of distinct elements nbout has tob known. 
    Obviously nbin is less than nbout*/
void condense(int nbin, int *arrayin, int nbout, int *arrayout);


/** @brief partition

    Basic patrtitionning tool, usually called from MatMapSetSize.
    Spread Euclideanly N on \param us numbers well ballanced.
    Note that, providing a communicator is not necessary.
    As a result, this function can also be called in a sequantial or OpenMP context.
    @todo move this function out of algo.h*/
int partition(int *fi, int *wi, int N, int me, int us);
