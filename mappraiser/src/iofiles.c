/**
 * @file   iofiles.c
 * @author Frederic Dauvergne
 * @date   November 2012
 * @last_update May 2019 by Hamza El Bouhargani
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fitsio.h>

#include "mappraiser/iofiles.h"

char *WORKDIR;
int IOVERBOSE = 0;

// Copied from HEALPix write_healpix_map routine with little tweaks to pass datatypes in arguments
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
                       "G = Galactic, E = ecliptic, C = celestial = equatorial", &status);

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
        fprintf(stderr, "%s (%d): System Cordinates are not correct"
                        "(Galactic,Ecliptic,Celestial=Equatorial). "
                        " Celestial system was set.\n",
                __FILE__, __LINE__);
}

static void util_fail_(const char *file, int line, const char *func,
                       const char *msg) {
    fprintf(stderr, "%s, %i (%s):\n%s\n", file, line, func, msg);
    exit(1);
}

int ioReadfile(int block_size, int part_id, unsigned int *point_data, double *signal) {
    int i;

    char p_vectorFile[256];
    char *p_vectorFileNameBasis = "point_data_0_";

    char s_vectorFile[256];
    //        char *s_vectorFileNameBasis = "signal_";
    char *s_vectorFileNameBasis = "pure_signal_";

    FILE *fp;

    sprintf(p_vectorFile, "%s%s%01d.dat", WORKDIR, p_vectorFileNameBasis, part_id);
    sprintf(s_vectorFile, "%s%s%01d.dat", WORKDIR, s_vectorFileNameBasis, part_id);

    if (IOVERBOSE > 1) {
        printf(" Pointing file name: %s\n", p_vectorFile);
        printf("   Signal file name: %s\n", s_vectorFile);
    }

    fp = fopen(p_vectorFile, "rb");
    fread(point_data, sizeof(unsigned int), block_size, fp);
    fclose(fp);

    fp = fopen(s_vectorFile, "rb");
    fread(signal, sizeof(double), block_size, fp);
    fclose(fp);

    return 0;
}

int ioReadfile_pol(int jump, int loop, int block_size, int part_id, unsigned int *point_data, double *signal,
                   double *pol_ang) {
    int i;

    char p_vectorFile[256];
    char *p_vectorFileNameBasis = "point_data_0_";

    char s_vectorFile[256];
    //        char *s_vectorFileNameBasis = "signal_";
    char *s_vectorFileNameBasis = "pure_signal_";

    char pol_vectorFile[256];
    char *pol_vectorFileNameBasis = "pol_angle_";

    FILE *fp;

    sprintf(p_vectorFile, "%s%s%01d.dat", WORKDIR, p_vectorFileNameBasis, part_id);
    sprintf(s_vectorFile, "%s%s%01d.dat", WORKDIR, s_vectorFileNameBasis, part_id);
    sprintf(pol_vectorFile, "%s%s%01d.dat", WORKDIR, pol_vectorFileNameBasis, part_id);

    if (IOVERBOSE > 1) {
        printf(" Pointing file name: %s\n", p_vectorFile);
        printf("   Signal file name: %s\n", s_vectorFile);
        printf(" Polarisation angles file name: %s\n", pol_vectorFile);
    }

    fp = fopen(p_vectorFile, "rb");
    fseek(fp, 3 * (jump + loop * block_size) * sizeof(unsigned int), SEEK_SET);
    fread(point_data, sizeof(unsigned int), 3 * block_size, fp);
    fclose(fp);

    fp = fopen(s_vectorFile, "rb");
    fseek(fp, (jump + loop * block_size) * sizeof(double), SEEK_SET);
    fread(signal, sizeof(double), block_size, fp);
    fclose(fp);

    fp = fopen(pol_vectorFile, "rb");
    fseek(fp, (jump + loop * block_size) * sizeof(double), SEEK_SET);
    fread(pol_ang, sizeof(double), block_size, fp);
    fclose(fp);

    return 0;
}

int ioReadTOAST_data(int jump, int loop, int block_size, int part_id, unsigned int *point_data, double *signal,
                     double *wghts) {
    int i;

    char p_vectorFile[256];
    char *p_vectorFileNameBasis = "pixels_";

    char s_vectorFile[256];
    //        char *s_vectorFileNameBasis = "signal_";
    char *s_vectorFileNameBasis = "pure_signal_";

    char wght_vectorFile[256];
    char *wght_vectorFileNameBasis = "weights_";

    FILE *fp;

    sprintf(p_vectorFile, "%s%s%01d.dat", WORKDIR, p_vectorFileNameBasis, part_id);
    sprintf(s_vectorFile, "%s%s%01d.dat", WORKDIR, s_vectorFileNameBasis, part_id);
    sprintf(wght_vectorFile, "%s%s%01d.dat", WORKDIR, wght_vectorFileNameBasis, part_id);

    if (IOVERBOSE > 1) {
        printf(" Pixels file name: %s\n", p_vectorFile);
        printf("   Signal file name: %s\n", s_vectorFile);
        printf(" Pointing weights file name: %s\n", wght_vectorFile);
    }

    fp = fopen(p_vectorFile, "rb");
    fseek(fp, 3 * (jump + loop * block_size) * sizeof(unsigned int), SEEK_SET);
    fread(point_data, sizeof(unsigned int), 3 * block_size, fp);
    fclose(fp);

    fp = fopen(s_vectorFile, "rb");
    fseek(fp, (jump + loop * block_size) * sizeof(double), SEEK_SET);
    fread(signal, sizeof(double), block_size, fp);
    fclose(fp);

    fp = fopen(wght_vectorFile, "rb");
    fseek(fp, 3 * (jump + loop * block_size) * sizeof(double), SEEK_SET);
    fread(wghts, sizeof(double), 3 * block_size, fp);
    fclose(fp);

    return 0;
}

int ioReadfilePure(int block_size, int part_id, unsigned int *point_data, double *signal) {
    int i;

    char p_vectorFile[256];
    char *p_vectorFileNameBasis = "point_data_0_";

    char s_vectorFile[256];
    //        char *s_vectorFileNameBasis = "signal_";
    char *s_vectorFileNameBasis = "pure_signal_";

    FILE *fp;

    sprintf(p_vectorFile, "%s%s%01d.dat", WORKDIR, p_vectorFileNameBasis, part_id);
    sprintf(s_vectorFile, "%s%s%01d.dat", WORKDIR, s_vectorFileNameBasis, part_id);

    printf(" Pointing file name: %s\n", p_vectorFile);
    printf("   Signal file name: %s\n", s_vectorFile);

    fp = fopen(p_vectorFile, "rb");
    fread(point_data, sizeof(unsigned int), block_size, fp);
    fclose(fp);

    fp = fopen(s_vectorFile, "rb");
    fread(signal, sizeof(double), block_size, fp);
    fclose(fp);

    return 0;
}

int ioReadrandom(int block_size, int part_id, unsigned int *point_data, double *signal) {
    int i;

    // Random generator:
    srand(part_id); // initialize the random generator
    for (i = 0; i < block_size; i++)
        signal[i] = 1.0 + (10 * ((double) rand()) / RAND_MAX - 1);

    for (i = 0; i < block_size; i++)
        point_data[i] = i;

    return 0;
}

int ioReadTpltzfile(int lambda, double *Tblock) {

    int i;

    char N_vectorFile[256];
    char *N_vectorFileNameBasis = "inv_tt_x3";

    FILE *fp;

    sprintf(N_vectorFile, "%s%s.bin", WORKDIR, N_vectorFileNameBasis);

    printf(" Block Toeplitz values file name: %s\n", N_vectorFile);

    fp = fopen(N_vectorFile, "rb");
    fread(Tblock, sizeof(double), lambda, fp);

    fclose(fp);

    return 0;
}

int ioReadTpltzrandom(int lambda, double *Tblock) {
    //        int i;
    // double lambdareduce=10;
    //
    //   srand (lambda); //init seed
    //
    // //input matrix definition of T
    //   for(i=1;i<lambdareduce;i++) {
    //     Tblock[i]= -1.0/((double) i);
    //   }
    //   for(i=lambdareduce;i<lambda;i++) {
    //     Tblock[i]= rand()/((double) 100*RAND_MAX);
    // }
    // for(i=1;i<lambda;i++)
    //   Tblock[i]= rand()/((double)RAND_MAX);

    Tblock[0] = 1;

    return 0;
}

int ioWritebinfile(int mapsize, int mappart_id, int *lstid, double *map, double *cond, int *lhits) {
    int i;

    char lstid_vectorFile[256];
    char lhits_vectorFile[256];
    char cond_vectorFile[256];
    char x_vectorFile[256];
    char *lstid_vectorFileNameBasis = "/global/cscratch1/sd/elbouha/data_TOAST/output/mapout_lstid_";
    char *lhits_vectorFileNameBasis = "/global/cscratch1/sd/elbouha/data_TOAST/output/mapout_lhits_";
    char *cond_vectorFileNameBasis = "/global/cscratch1/sd/elbouha/data_TOAST/output/mapout_lcond_";
    char *x_vectorFileNameBasis = "/global/cscratch1/sd/elbouha/data_TOAST/output/mapout_";

    FILE *fp;

    sprintf(lstid_vectorFile, "%s%01d.dat", lstid_vectorFileNameBasis, mappart_id);
    sprintf(lhits_vectorFile, "%s%01d.dat", lhits_vectorFileNameBasis, mappart_id);
    sprintf(cond_vectorFile, "%s%01d.dat", cond_vectorFileNameBasis, mappart_id);
    sprintf(x_vectorFile, "%s%01d.dat", x_vectorFileNameBasis, mappart_id);

    //        printf(" Map file name: %s\n", lstid_vectorFile);
    //        printf(" Map file name: %s\n", x_vectorFile);

    fp = fopen(lstid_vectorFile, "wb");
    fwrite(lstid, sizeof(int), mapsize, fp);
    fclose(fp);

    fp = fopen(lhits_vectorFile, "wb");
    fwrite(lhits, sizeof(int), (int) (mapsize / 3), fp);
    fclose(fp);

    fp = fopen(cond_vectorFile, "wb");
    fwrite(cond, sizeof(double), (int) (mapsize / 3), fp);
    fclose(fp);

    fp = fopen(x_vectorFile, "wb");
    fwrite(map, sizeof(double), mapsize, fp);
    fclose(fp);

    return 0;
}

int ioReadbinfile(int mapsize, int mappart_id, int *lstid, double *map, double *cond, int *lhits) {

    int i;

    char lstid_vectorFile[256];
    char lhits_vectorFile[256];
    char cond_vectorFile[256];
    char x_vectorFile[256];
    char *lstid_vectorFileNameBasis = "mapout_lstid_";
    char *lhits_vectorFileNameBasis = "mapout_lhits_";
    char *cond_vectorFileNameBasis = "mapout_lcond_";
    char *x_vectorFileNameBasis = "mapout_";

    FILE *fp;

    sprintf(lstid_vectorFile, "%s%01d.dat", lstid_vectorFileNameBasis, mappart_id);
    sprintf(lhits_vectorFile, "%s%01d.dat", lhits_vectorFileNameBasis, mappart_id);
    sprintf(cond_vectorFile, "%s%01d.dat", cond_vectorFileNameBasis, mappart_id);
    sprintf(x_vectorFile, "%s%01d.dat", x_vectorFileNameBasis, mappart_id);

    printf(" Map id file name: %s\n", lstid_vectorFile);
    printf(" Hits map file name: %s\n", lhits_vectorFile);
    printf(" Condition map file name: %s\n", cond_vectorFile);
    printf(" Map file name: %s\n", x_vectorFile);

    fp = fopen(lstid_vectorFile, "rb");
    fread(lstid, sizeof(int), mapsize, fp);
    fclose(fp);

    fp = fopen(lhits_vectorFile, "rb");
    fread(lhits, sizeof(int), (int) (mapsize / 3), fp);
    fclose(fp);

    fp = fopen(cond_vectorFile, "rb");
    fread(cond, sizeof(double), (int) (mapsize / 3), fp);
    fclose(fp);

    fp = fopen(x_vectorFile, "rb");
    fread(map, sizeof(double), mapsize, fp);
    fclose(fp);

    return 0;
}
