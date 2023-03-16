#ifndef MAPPRAISER_IOFILES_H
#define MAPPRAISER_IOFILES_H

int ioWritebinfile(int mapsize, int mappart_id, int *lstid, double *map, double *cond, int *lhits);

int ioReadbinfile(int mapsize, int mappart_id, int *lstid, double *map, double *cond, int *lhits);

// Import HEALPix routines with little tweaks to pass datatypes in arguments
static void util_fail_(const char *file, int line, const char *func,
                       const char *msg);

#if defined(__GNUC__)
#define UTIL_FUNC_NAME__ __func__
#else
#define UTIL_FUNC_NAME__ "unknown"
#endif
#define UTIL_ASSERT(cond, msg) \
    if (!(cond))               \
    util_fail_(__FILE__, __LINE__, UTIL_FUNC_NAME__, msg)
#define UTIL_FAIL(msg) \
    util_fail_(__FILE__, __LINE__, UTIL_FUNC_NAME__, msg)

static void setCoordSysHP(char coordsys, char *coordsys9);

static void printerror(int status);

void write_map(void *signal, int type, long nside, const char *filename,
               char nest, const char *coordsys);

#endif //MAPPRAISER_IOFILES_H
