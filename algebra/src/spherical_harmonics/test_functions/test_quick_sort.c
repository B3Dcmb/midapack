#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
// #include <mkl.h>
// #include "fitsio.h"
// #include <mpi.h>
#include <unistd.h>

void quick_sort(int *indices, int left, int right);
void quick_sort_with_indices(int *indices, int *index_of_indices, int left, int right);
int projection(int *values_in, int number_values, int *projector_in2out, int *values_out);

int main()
{
    struct timespec start, end;
    uint64_t delta_us;

    int size_array = 100000;

    int *worst_case = (int *)malloc(size_array*sizeof(int));
    int *best_case = (int *)malloc(size_array*sizeof(int));

    int i;
    for(i=0; i<size_array; i++){
        worst_case[i] = size_array - i;
        best_case[i] = i;
    }

    int *worst_case_copy = (int *)malloc(size_array*sizeof(int));
    int *best_case_copy = (int *)malloc(size_array*sizeof(int));

    int *worst_case_copy_2 = (int *)malloc(size_array*sizeof(int));
    int *best_case_copy_2 = (int *)malloc(size_array*sizeof(int));

    int *projector_best_case = (int *)malloc(size_array*sizeof(int));
    int *projector_worst_case = (int *)malloc(size_array*sizeof(int));

    memcpy(worst_case_copy, worst_case, size_array);
    memcpy(best_case_copy, best_case, size_array);
    memcpy(worst_case_copy_2, worst_case, size_array);
    memcpy(best_case_copy_2, best_case, size_array);
    memcpy(projector_best_case, best_case, size_array);
    memcpy(projector_worst_case, best_case, size_array);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    quick_sort(best_case, 0, size_array-1);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Time quick_sort best_case : %ld ms \n", delta_us);
    
    
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    quick_sort(worst_case, 0, size_array-1);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Time quick_sort worst_case : %ld ms \n", delta_us);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    quick_sort_with_indices(best_case_copy, projector_best_case, 0, size_array-1);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Time quick_sort_with_indices best_case : %ld ms \n", delta_us);
    
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    quick_sort_with_indices(worst_case_copy, projector_worst_case, 0, size_array-1);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Time quick_sort_with_indices worst_case : %ld ms \n", delta_us);

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    projection(best_case_copy_2, size_array, projector_best_case, best_case);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Time projection worst_case : %ld ms \n", delta_us);
    
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    projection(worst_case_copy_2, size_array, projector_worst_case, worst_case);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    delta_us = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Time projection worst_case : %ld ms \n", delta_us);

    free(best_case);
    free(best_case_copy);
    free(best_case_copy_2);
    free(projector_best_case);
    
    free(worst_case);
    free(worst_case_copy);
    free(worst_case_copy_2);
    free(projector_worst_case);

    return 0;
}

void quick_sort(int *indices, int left, int right){
  int pivot;
  int tmp, key;
  int i,j;
  if (left<right){
    key=indices[left];
    i=left+1;
    j=right;
    while(i<=j){
      while((i<=right) && (indices[i]<=key))  i++;
      while ((j>left) && (indices[j]>key))   j--;
      if(i<j){
        tmp=indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
        i++;
        j--;
      }
    }
    tmp=indices[left];
    indices[left] = indices[j];
    indices[j] = tmp;
    pivot = j;
    quick_sort(indices, left, pivot-1);
    quick_sort(indices, pivot+1, right);
  }
}

void quick_sort_with_indices(int *indices, int *index_of_indices, int left, int right){
  int pivot;
  int tmp, key;
  int tmp_indices;
  int i,j;
  if (left<right){
    key=indices[left];
    i=left+1;
    j=right;
    while(i<=j){
      while((i<=right) && (indices[i]<=key))  i++;
      while ((j>left) && (indices[j]>key))   j--;
      if(i<j){
        tmp=indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
        tmp_indices = index_of_indices[i];
        index_of_indices[i] = index_of_indices[j];
        index_of_indices[j] = tmp_indices;
        i++;
        j--;
      }
    }
    tmp=indices[left];
    indices[left] = indices[j];
    indices[j] = tmp;
    tmp_indices = index_of_indices[left];
    index_of_indices[left] = index_of_indices[j];
    index_of_indices[j] = tmp_indices;
    pivot = j;
    quick_sort_with_indices(indices, index_of_indices, left, pivot-1);
    quick_sort_with_indices(indices, index_of_indices, pivot+1, right);
  }
}

int projection(int *values_in, int number_values, int *projector_in2out, int *values_out)
{
  int i;
  for (i=0; i<number_values; i++)
    values_out[i] = values_in[projector_in2out[i]];
  
  return 0;
}
