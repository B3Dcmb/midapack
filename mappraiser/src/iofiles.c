/**
 * @file iofiles.c
 * @brief IO routines for writing products of the map maker
 * @author Frederic Dauvergne
 * @date November 2012
 * @last_update May 2019 by Hamza El Bouhargani
 */

#include <fitsio.h>
#include <stdio.h>
#include <string.h>

#include <mappraiser/iofiles.h>

static void util_fail_(const char *file, int line, const char *func,
                       const char *msg) {
    fprintf(stderr, "%s, %i (%s):\n%s\n", file, line, func, msg);
    exit(1);
}

#if defined(__GNUC__)
#define UTIL_FUNC_NAME__ __func__
#else
#define UTIL_FUNC_NAME__ "unknown"
#endif
#define UTIL_ASSERT(cond, msg)                                                 \
    if (!(cond))                                                               \
    util_fail_(__FILE__, __LINE__, UTIL_FUNC_NAME__, msg)
#define UTIL_FAIL(msg) util_fail_(__FILE__, __LINE__, UTIL_FUNC_NAME__, msg)

static void printerror(int status);
static void setCoordSysHP(char coordsys, char *coordsys9);

// Copied from HEALPix write_healpix_map routine with little tweaks to pass
// datatypes in arguments
void write_map(void *signal, int type, long nside, const char *filename,
               char nest, const char *coordsys) {
    fitsfile *fptr; /* pointer to the FITS file, defined in fitsio.h */
    int status = 0, hdutype;

    long naxes[] = {0, 0};

    char order[9]; /* HEALPix ordering */
    char *ttype[] = {"SIGNAL"};
    char *tform[] = {"1D"};
    char *tunit[] = {" "};
    char coordsys9[9];

    /* create new FITS file */
    fits_create_file(&fptr, filename, &status);
    fits_create_img(fptr, SHORT_IMG, 0, naxes, &status);
    fits_write_date(fptr, &status);
    fits_movabs_hdu(fptr, 1, &hdutype, &status);
    fits_create_tbl(fptr, BINARY_TBL, 12L * nside * nside, 1, ttype, tform,
                    tunit, "BINTABLE", &status);
    fits_write_key(fptr, TSTRING, "PIXTYPE", "HEALPIX", "HEALPIX Pixelisation",
                   &status);

    strcpy(order, nest ? "NESTED  " : "RING    ");
    fits_write_key(fptr, TSTRING, "ORDERING", order,
                   "Pixel ordering scheme, either RING or NESTED", &status);
    fits_write_key(fptr, TLONG, "NSIDE", &nside,
                   "Resolution parameter for HEALPIX", &status);

    UTIL_ASSERT(strlen(coordsys) >= 1, "bad ccordsys value");
    setCoordSysHP(coordsys[0], coordsys9);
    fits_write_key(fptr, TSTRING, "COORDSYS", coordsys9,
                   "Pixelisation coordinate system", &status);

    fits_write_comment(fptr,
                       "G = Galactic, E = ecliptic, C = celestial = equatorial",
                       &status);

    fits_write_col(fptr, type, 1, 1, 1, 12 * nside * nside, signal, &status);
    fits_close_file(fptr, &status);
    printerror(status);
}

static void printerror(int status) {
    if (status == 0)
        return;

    fits_report_error(stderr, status);
    UTIL_FAIL("FITS error");
}

static void setCoordSysHP(char coordsys, char *coordsys9) {
    strcpy(coordsys9, "C       ");
    if (coordsys == 'G')
        strcpy(coordsys9, "G       ");
    else if (coordsys == 'E')
        strcpy(coordsys9, "E       ");
    else if ((coordsys != 'C') && (coordsys != 'Q'))
        fprintf(stderr,
                "%s (%d): System Cordinates are not correct"
                "(Galactic,Ecliptic,Celestial=Equatorial). "
                " Celestial system was set.\n",
                __FILE__, __LINE__);
}

int saveArrayToFile(const char *filename, const void *array, size_t arraySize,
                    size_t elementSize) {
    // Open the file in binary write mode
    FILE *file = fopen(filename, "wb");

    if (file == NULL) {
        perror("Error opening the file");
        return 1;
    }

    size_t elementsWritten = fwrite(array, elementSize, arraySize, file);

    if (elementsWritten != arraySize) {
        perror("Error writing to the file");
        fclose(file);
        return 1;
    }

    fclose(file);

    return 0;
}
