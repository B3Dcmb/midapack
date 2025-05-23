//
// Created by sbiquard on 5/17/23.
//

#include "utils.h"

int create_directory(const char *path, mode_t mode) {
    struct stat st;

    // Check if the directory already exists
    if (stat(path, &st) == 0) {
        if (S_ISDIR(st.st_mode)) {
            // Directory already exists
            return 0;
        } else {
            // File exists with the same name, handle the error
            fprintf(stderr, "Error: A file with the same name exists: %s\n",
                    path);
            return -1;
        }
    }

    // Directory doesn't exist, create it
    if (mkdir(path, mode) == -1) {
        perror("mkdir");
        return -1;
    }

    return 0;
}

void fillArrayFromFile(const char *filename, void *array, size_t size,
                       size_t elementSize) {
    FILE *file = fopen(filename, "rb"); // Open the file in binary mode

    if (file == NULL) {
        fprintf(stderr, "Failed to open the file %s\n", filename);
        return;
    }

    // Read binary data from the file into the array
    size_t elementsRead = fread(array, elementSize, size, file);

    if (elementsRead != size) {
        printf("Failed to read the entire array from the file %s\n", filename);
        printf("Elements read: %ld / %ld\n", elementsRead, size);
    }

    fclose(file); // Close the file
}

__attribute__((unused)) void printArray(const void *array, size_t size,
                                        size_t elementSize,
                                        PrintFunc printElement) {
    const char *base =
        (const char *)array; // Treat the array as an array of characters

    for (size_t i = 0; i < size; ++i) {
        const void *element =
            base +
            i * elementSize; // Calculate the address of the current element
        printElement(
            element); // Call the callback function to print the element
    }
}

// Example callback function for printing integers
__attribute__((unused)) void printInt(const void *element) {
    const int *intValue =
        (const int *)element; // Treat the element as an integer
    printf("%d ", *intValue);
}

__attribute__((unused)) void printDouble(const void *element) {
    const double *doubleValue =
        (const double *)element; // Treat the element as a double
    printf("%lf ", *doubleValue);
}

__attribute__((unused)) void printFloat(const void *element) {
    const float *floatValue =
        (const float *)element; // Treat the element as a float
    printf("%f ", *floatValue);
}
