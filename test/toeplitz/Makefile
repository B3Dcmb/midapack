DIR = /global/homes/a/ajamal/CMB_TOAST/midapack
TOEPLITZ = $(DIR)/src/toeplitz/
GCC = gcc
# MPICC = /global/homes/e/elbouha/.conda/envs/hamzaenv/bin/mpicc
MPICC = cc
CC = $(MPICC) -qopenmp $(OPT)
OPT =  -D W_MPI -D W_OPENMP -D MPI -D OPENMP -O3
INC = -I$(DIR)/include -I/opt/cray/pe/fftw/3.3.8.3/haswell/include
LIB = -L$(DIR)/lib -lmidapack
FFTW_LIB = -L/opt/cray/pe/fftw/3.3.8.3/haswell/lib -lfftw3 -lfftw3_threads
FFTW     = -D FFTW_MULTITHREADING



SRC = ../../src/toeplitz/


##############   t e s t s   ############

all:
	make example

example:
	@echo "start compiling toeplitz examples ........"
	make mpi_stmm
	make mpi_stbmm
	make mpi_gstbmm


#test:
#	@echo "start compiling toeplitz tests ........"
#	make test_mpi_stmm

mpi_stmm: $(SRC)
	$(CC) testexple_mpi_stmm.c $(INC) $(LIB) $(FFTW_LIB) -lm -o testexple_mpi_stmm

mpi_stbmm: $(SRC)
	$(CC) testexple_mpi_stbmm.c $(INC) $(LIB) $(FFTW_LIB) -lm -o testexple_mpi_stbmm

mpi_gstbmm: $(SRC)
	$(CC) testexple_mpi_gstbmm.c $(INC) $(LIB) $(FFTW_LIB) -lm -o testexple_mpi_gstbmm



##############   u t i l s   #############

clean:
	rm -f *.o ;
