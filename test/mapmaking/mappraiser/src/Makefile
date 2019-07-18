DIR = /global/homes/e/elbouha/midapack
TOEPLITZ = $(DIR)/src/toeplitz/
MPICC = cc
CC = $(MPICC) -qopenmp -g -Wall -fPIC $(OPT)
INC = -I$(DIR)/include
LIB = -L$(DIR)/lib -lmidapack
FFTW_LIB = -L/opt/cray/pe/fftw/3.3.6.3/haswell/lib -lfftw3 -lfftw3_threads
FFTW_INC = -I/opt/cray/pe/fftw/3.3.6.3/haswell/include
# FFTW_LIB = -L/opt/cray/pe/fftw/3.3.6.3/haswell/lib/libfftw3.a /opt/cray/pe/fftw/3.3.6.3/haswell/lib/libfftw3_threads.a

SRC = $(DIR)/src/toeplitz/

MAPMAT = $(DIR)/src/mapmat/
F90 = ftn
OPT =  -D W_MPI -D W_OPENMP -O3

#supplementary libraries
HEALPIX_INC = -I$(HOME)/Healpix_3.50/include
CFITSIO_INC = -I/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/include

HEALPIX_LIB = -L$(HOME)/Healpix_3.50/lib -lchealpix -lhealpix
CFITSIO_LIB = -L/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/lib -lcfitsio

##############   t e s t s   ############

all:
	make mappraiser

mappraiser: $(TOEPLITZ) $(MAPMAT)
	@echo "start compiling MAPPRAISER ........"
	make createToeplitz.o \
	precond.o \
	pcg_true.o \
	iofiles.o \
	mappraiser.o \

createToeplitz.o: $(SRC)
	$(CC) -c createToeplitz.c $(INC) $(LIB) $(FFTW_INC) $(FFTW_LIB)

precond.o: $(MAPMAT)
	$(CC) -c precond.c $(INC) $(LIB) $(FFTW_INC) $(FFTW_LIB)

pcg_true.o: $(SRC) $(MAPMAT)
	$(CC) -c pcg_true.c $(INC) $(LIB) $(FFTW_INC) $(FFTW_LIB)

iofiles.o:
	$(CC) -c iofiles.c $(INC) $(LIB) $(CFITSIO_INC) $(CFITSIO_LIB) $(FFTW_INC) $(FFTW_LIB)

mappraiser.o: $(SRC) $(MAPMAT)
	$(CC) -c mappraiser.c $(INC) $(LIB) $(CFITSIO_INC) $(FFTW_INC) $(FFTW_LIB)

# s4cmb_pipeline: createToeplitz.o precond.o s4cmb_pipeline.c iofiles.o pcg_true.o
# 	$(CC) s4cmb_pipeline.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -lm -o s4cmb_pcg -Wl,@/global/homes/e/elbouha/midapack/test/mapmaking/allinea-profiler.ld

# lib: $(LIB) ./src
# 	cc -shared createToeplitz.o precond.o pcg_true.o iofiles.o mappraiser.o $(LIB) -o ./lib/libmappraiser.so
##############   u t i l s   #############

clean:
	rm -f *.o ;
