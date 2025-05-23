#ifndef MAPPRAISER_IOFILES_H
#define MAPPRAISER_IOFILES_H

#include <stdlib.h>

void write_map(void *signal, int type, long nside, const char *filename,
               char nest, const char *coordsys);

int saveArrayToFile(const char *filename, const void *array, size_t arraySize,
                    size_t elementSize);

#endif // MAPPRAISER_IOFILES_H
