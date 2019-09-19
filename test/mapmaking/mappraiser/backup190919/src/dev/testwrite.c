#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "fitsio.h"

static void util_fail_ (const char *file, int line, const char *func,
  const char *msg);

#if defined (__GNUC__)
#define UTIL_FUNC_NAME__ __func__
#else
#define UTIL_FUNC_NAME__ "unknown"
#endif
#define UTIL_ASSERT(cond,msg) \
  if(!(cond)) util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)
#define UTIL_FAIL(msg) \
  util_fail_(__FILE__,__LINE__,UTIL_FUNC_NAME__,msg)

static void setCoordSysHP(char coordsys,char *coordsys9);

static void printerror (int status);


void write_map (void *signal, int type, long nside, const char *filename,
  char nest, const char *coordsys);

int main(int argc, char *argv[])
{
  int i;
  double *map;
  map = (double *) malloc(12*sizeof(double));

  float *map2;
  map2 = (float *) malloc(12*sizeof(float));

  int* hits;
  hits = (int *) malloc(12*sizeof(int));

  for(i=0;i<12;i++){
    map[i] = i;
    map2[i] = i;
    hits[i] = 2*i;
  }

  char *mapname = "map.fits";
  char *map2name = "map2.fits";
  char *hitsname = "hits.fits";
  char nest = 1;
  char *cordsys = "C";
  int type1 = TDOUBLE;
  int type2 = TFLOAT;
  int type3 = TINT;
  write_map(map, type1, 1, mapname, nest, cordsys);
  write_map(map2, type2, 1, map2name, nest, cordsys);
  write_map(hits, type3, 1, hitsname, nest, cordsys);


  return 0;
}

void write_map (void *signal, int type, long nside, const char *filename,
  char nest, const char *coordsys)
  {
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status=0, hdutype;

  long naxes[] = {0,0};

  char order[9];                 /* HEALPix ordering */
  char *ttype[] = { "SIGNAL" };
  char *tform[] = { "1E" };
  char *tunit[] = { " " };
  char coordsys9[9];

  /* create new FITS file */
  fits_create_file(&fptr, filename, &status);
  fits_create_img(fptr, SHORT_IMG, 0, naxes, &status);
  fits_write_date(fptr, &status);
  fits_movabs_hdu(fptr, 1, &hdutype, &status);
  fits_create_tbl( fptr, BINARY_TBL, 12L*nside*nside, 1, ttype, tform,
                        tunit, "BINTABLE", &status);
  fits_write_key(fptr, TSTRING, "PIXTYPE", "HEALPIX", "HEALPIX Pixelisation",
    &status);

  strcpy(order, nest ? "NESTED  " : "RING    ");
  fits_write_key(fptr, TSTRING, "ORDERING", order,
    "Pixel ordering scheme, either RING or NESTED", &status);
  fits_write_key(fptr, TLONG, "NSIDE", &nside,
    "Resolution parameter for HEALPIX", &status);

  UTIL_ASSERT(strlen(coordsys)>=1,"bad ccordsys value");
  setCoordSysHP(coordsys[0],coordsys9);
  fits_write_key(fptr, TSTRING, "COORDSYS", coordsys9,
    "Pixelisation coordinate system", &status);

  fits_write_comment(fptr,
    "G = Galactic, E = ecliptic, C = celestial = equatorial", &status);

  fits_write_col(fptr, type, 1, 1, 1, 12*nside*nside, signal,&status);
  fits_close_file(fptr, &status);
  printerror(status);
  }

static void printerror (int status)
{
  if (status==0) return;

  fits_report_error(stderr, status);
  UTIL_FAIL("FITS error");
}

static void setCoordSysHP(char coordsys,char *coordsys9)
{
  strcpy(coordsys9,"C       ");
  if (coordsys=='G')
    strcpy (coordsys9,"G       ");
  else if (coordsys=='E')
    strcpy (coordsys9,"E       ");
  else if ((coordsys!='C')&&(coordsys!='Q'))
    fprintf(stderr, "%s (%d): System Cordinates are not correct"
                    "(Galactic,Ecliptic,Celestial=Equatorial). "
                    " Celestial system was set.\n", __FILE__, __LINE__);
}

static void util_fail_ (const char *file, int line, const char *func,
  const char *msg)
  {
  fprintf(stderr,"%s, %i (%s):\n%s\n",file,line,func,msg);
  exit(1);
  }
