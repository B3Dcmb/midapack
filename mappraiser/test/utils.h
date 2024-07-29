//
// Created by sbiquard on 5/17/23.
//

#ifndef MAPPRAISER_TEST_UTILS_H
#define MAPPRAISER_TEST_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

// Create a directory if it does not exist
int create_directory(const char *directoryName, mode_t mode);

// Helper functions for arrays

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

void fillArrayFromFile(const char *filename, void *array, size_t size,
                       size_t elementSize);

// Function pointer type for printing elements
typedef void (*PrintFunc)(const void *);
__attribute__((unused)) void printArray(const void *array, size_t size,
                                        size_t elementSize,
                                        PrintFunc printElement);
__attribute__((unused)) void printInt(const void *element);
__attribute__((unused)) void printDouble(const void *element);
__attribute__((unused)) void printFloat(const void *element);

#ifdef __cplusplus
}
#endif

#endif // MAPPRAISER_TEST_UTILS_H
