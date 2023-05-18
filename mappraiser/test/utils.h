//
// Created by sbiquard on 5/17/23.
//

#ifndef MAPPRAISER_TEST_UTILS_H
#define MAPPRAISER_TEST_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

// Helper functions for directory manipulation

int directoryExists(const char *directoryName);
int createDirectory(const char *directoryName);

// Helper functions for arrays

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

void fillArrayFromFile(const char *filename, void *array, size_t size, size_t elementSize);

typedef void (*PrintFunc)(const void *); // Function pointer type for printing elements
__attribute__((unused)) void printArray(const void *array, size_t size, size_t elementSize, PrintFunc printElement);
__attribute__((unused)) void printInt(const void *element);
__attribute__((unused)) void printDouble(const void *element);
__attribute__((unused)) void printFloat(const void *element);

#endif // MAPPRAISER_TEST_UTILS_H
