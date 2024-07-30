#ifndef MEMUTILS_H
#define MEMUTILS_H

#include <stdio.h>
#include <stdlib.h>

/* safe memory allocation macros, inspired by
 * https://stackoverflow.com/questions/16297142/safe-malloc-realloc-wrapping-the-call-into-a-macro
 */

static void *safe_malloc(size_t n, const char *file, unsigned long line) {
    void *p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "[%s:%lu]Out of memory(%lu bytes) at %s\n", file, line,
                n, file);
        exit(EXIT_FAILURE);
    }
    return p;
}
#define SAFEMALLOC(n) safe_malloc(n, __FILE__, __LINE__)

static void *safe_calloc(size_t nmemb, size_t size, const char *file,
                         unsigned long line) {
    void *p = calloc(nmemb, size);
    if (p == NULL) {
        fprintf(stderr,
                "[%s:%lu] Out of memory for %lu elements of size %lu at %s\n",
                file, line, nmemb, size, file);
        exit(EXIT_FAILURE);
    }
    return p;
}
#define SAFECALLOC(nmemb, size) safe_calloc(nmemb, size, __FILE__, __LINE__)

static void *safe_realloc(void *ptr, size_t new_size, const char *file,
                          unsigned long line) {
    void *new_ptr = realloc(ptr, new_size);
    if (new_ptr == NULL) {
        fprintf(stderr,
                "[%s:%lu] Out of memory when reallocating %lu bytes at %s\n",
                file, line, new_size, file);
        free(ptr); // Free the original memory to avoid leaks
        exit(EXIT_FAILURE);
    }
    return new_ptr;
}
#define SAFEREALLOC(ptr, new_size)                                             \
    safe_realloc(ptr, new_size, __FILE__, __LINE__)

// free a pointer and set if to NULL to avoid double free
#define FREE(ptr) ((void)(free(ptr), (ptr) = NULL))

#endif // MEMUTILS_H
