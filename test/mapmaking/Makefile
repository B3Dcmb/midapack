DIR = /global/homes/e/elbouha/midapack
TOEPLITZ = $(DIR)/src/toeplitz/
GCC = gcc
# MPICC = /global/homes/e/elbouha/.conda/envs/hamzaenv/bin/mpicc
MPICC = cc
CC = $(MPICC) -qopenmp -g $(OPT)
INC = -I/opt/cray/pe/fftw/3.3.6.2/haswell/include
LIB = -L$(DIR)/lib -lmidapack
# FFTW_LIB = -L/opt/cray/pe/fftw/3.3.6.2/haswell/lib -lfftw3 -lfftw3_threads
#FFTW     = -D FFTW_MULTITHREADING
FFTW_LIB = /opt/cray/pe/fftw/3.3.6.2/haswell/lib/libfftw3.a /opt/cray/pe/fftw/3.3.6.2/haswell/lib/libfftw3_threads.a

SRC = ../../src/toeplitz/

MAPMAT = $(DIR)/src/mapmat/
F90 = ftn
OPT =  -D W_MPI -D W_OPENMP -O3

#supplementary libraries
# PTSCOTCH_INC = -I/home/cargemel/scotch_5.1.12/include/
HEALPIX_INC = -I$(HOME)/Healpix_3.31/include
CFITSIO_INC = -I/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/include
# CFITSIO_INC = -I/global/common/sw/cray/cnl6/ivybridge/cfitsio/3.410/intel/17.0.2.174/hgce4lk/include

# PTSCOTCH_LIB = -L/home/cargemel/scotch_5.1.12/lib/ -lptscotch -lptscotcherr
HEALPIX_LIB = -L$(HOME)/Healpix_3.31/lib -lchealpix -lhealpix
CFITSIO_LIB = -L/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/lib -lcfitsio
# CFITSIO_LIB = -L/global/common/sw/cray/cnl6/ivybridge/cfitsio/3.410/intel/17.0.2.174/hgce4lk/lib -lcfitsio


INTEL_DIR = /opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/
INTEL_LIB = $(INTEL_DIR)libsvml.a $(INTEL_DIR)libirc.a $(INTEL_DIR)libimf.a $(INTEL_DIR)libifcore.a $(INTEL_DIR)libiomp5.a

##############   t e s t s   ############

all:
	make clean
	make example_pcg
	make exemple_iofilepix
	make exemple_iofilepix_diff
	make exemple_iofilepix_diff2
	make mk_unpointing

example_pcg: $(TOEPLITZ) $(MAPMAT)
	@echo "start compiling pcg example ........"
	make createToeplitz.o \
	precond.o \
	pcg_true.o \
	iofiles.o \
	test_pcg \
	test_pcg_pol \
	toast_pipeline \
	test_pcg_rand \

exemple_iofilepix:
	make map2gif.o \
	mk_iofilespix \
	mk_iofilespix_pol

exemple_iofilepix_diff:
	make map2gif.o \
	mk_iofilespix_diff

exemple_iofilepix_diff2:
	make map2gif.o \
	mk_iofilespix_diff2

exemple_unpointing:
	make map2gif.o \
	mk_unpointing

#

createToeplitz.o: $(SRC)
	$(CC) -c createToeplitz.c $(INC) $(LIB) $(FFTW_LIB)

precond.o: $(SRC)
	$(CC) -c precond.c $(INC) $(LIB)


pcg_true.o: $(SRC)
	$(CC) -c pcg_true.c $(INC) $(LIB)

pcg_true_alt.o: $(SRC)
	$(CC) -c pcg_true_alt.c $(INC) $(LIB)

pcg_like.o: $(SRC)
	$(CC) -c pcg_like.c $(INC) $(LIB)

iofiles.o: $(SRC)
	$(CC) -c iofiles.c $(INC) $(LIB)


map2gif.o : map2gif.f90
	$(F90) -c map2gif.f90 -o map2gif.o $(HEALPIX_INC)

mk_iofilespix : mk_iofilespix.c map2gif.o iofiles.o
	$(CC) mk_iofilespix.c -o mk_iofilespix $(HEALPIX_INC) $(INC) $(LIB) map2gif.o iofiles.o -lm -lintlc -Wl,--start-group $(INTEL_LIB) $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB) -Wl,--end-group
mk_iofilespix_pol : mk_iofilespix_pol.c iofiles.o
	$(CC) mk_iofilespix_pol.c -o mk_iofilespix_pol $(INC) $(LIB) iofiles.o -lm
mk_iofilespix_diff : mk_iofilespix_diff.c map2gif.o iofiles.o
	$(CC) mk_iofilespix_diff.c -o mk_iofilespix_diff $(HEALPIX_INC) $(INC) $(LIB) map2gif.o iofiles.o -lm -lintlc -Wl,--start-group $(INTEL_LIB) $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB) -Wl,--end-group

mk_iofilespix_diff2 : mk_iofilespix_diff2.c map2gif.o iofiles.o
	$(CC) mk_iofilespix_diff2.c -o mk_iofilespix_diff2 $(HEALPIX_INC) $(INC) $(LIB) map2gif.o iofiles.o -lm -lintlc -Wl,--start-group $(INTEL_LIB) $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB) -Wl,--end-group

mk_unpointing : mk_unpointing.c map2gif.o iofiles.o
	$(CC) mk_unpointing.c -o mk_unpointing $(HEALPIX_INC) $(INC) $(LIB) map2gif.o iofiles.o -lm -lintlc -Wl,--start-group $(INTEL_LIB) $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB) -Wl,--end-group
#sky_scan : sky_scan.c $(MAPMAT) algo.o map2gif.o
#	$(CC) sky_scan.c -o sky_scan $(INC) $(LIB) algo.o map2gif.o -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.4/ -lgfortran $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB)

test_pcg: createToeplitz.o precond.o test_pcg.c iofiles.o pcg_true.o
	$(CC) test_pcg.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -lm -o pcg

test_pcg_pol: createToeplitz.o precond.o test_pcg_polar.c iofiles.o pcg_true.o
		$(CC) test_pcg_polar.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -lm -o pcg_pol

toast_pipeline: createToeplitz.o precond.o toast_pipeline.c iofiles.o pcg_true.o
		$(CC) toast_pipeline.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -lm -o toast_pcg

test_pcg_rand: createToeplitz.o precond.o test_pcg_rand.c iofiles.o pcg_true.o
	$(CC) test_pcg_rand.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -lm -o pcg_rand

test_com: test_com.c iofiles.o
			$(CC) test_com.c iofiles.o $(INC) $(LIB) $(FFTW_LIB) -lm -o test_com

##############   u t i l s   #############

clean:
	rm -f *.o ;
