GCC = gcc
# MPICC = /global/homes/e/elbouha/.conda/envs/hamzaenv/bin/mpicc
MPICC = cc
OPT =  -D W_MPI -D W_OPENMP -O3
CC = $(MPICC) -qopenmp -g $(OPT)
FFTW_INC = -I/opt/cray/pe/fftw/3.3.6.3/haswell/include
##############   s r c   ############

all:
	make toeplitz.o \
	toeplitz_seq.o \
	toeplitz_block.o \
	toeplitz_gappy.o \
	toeplitz_nofft.o \
	toeplitz_rshp.o \
	toeplitz_params.o \
	toeplitz_utils.o \
	toeplitz_wizard.o


toeplitz.o: toeplitz.c toeplitz_nofft.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz.c

toeplitz_seq.o: toeplitz_seq.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz_seq.c

toeplitz_block.o: toeplitz_block.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz_block.c

toeplitz_gappy.o: toeplitz_gappy.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz_gappy.c

toeplitz_nofft.o: toeplitz_nofft.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz_nofft.c

toeplitz_rshp.o: toeplitz_rshp.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz_rshp.c

toeplitz_params.o: toeplitz_params.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz_params.c

toeplitz_utils.o: toeplitz_utils.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz_utils.c

toeplitz_wizard.o: toeplitz_wizard.c
	$(CC) $(OPTC) $(FFTW_INC) -c toeplitz_wizard.c

########   u t i l s   ########

clean:
	rm -f *.o ;
