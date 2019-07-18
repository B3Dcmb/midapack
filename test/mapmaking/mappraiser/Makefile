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
	make map
	make library

map: ./src/
	make -C ./src/

library:
	cc -shared ./src/createToeplitz.o ./src/precond.o ./src/pcg_true.o ./src/iofiles.o ./src/mappraiser.o $(LIB) $(CFITSIO_LIB) -o ./lib/libmappraiser.so
##############   u t i l s   #############

clean:
	make clean -C ./src/
	rm -f lib/*midapack.so
