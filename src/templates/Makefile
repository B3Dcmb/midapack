MIDAPACK_ROOT=$(PREFIX)/midapack
MIDAPACK_OBJ=$(MIDAPACK_ROOT)/obj
TEMPLATES_OBJ=$(MIDAPACK_OBJ)/templates

GCC = gcc
# MPICC = /global/homes/e/elbouha/.conda/envs/hamzaenv/bin/mpicc
MPICC = cc
OPT =  -D W_MPI -D W_OPENMP -O3
CC = $(MPICC) -qopenmp -g -Wall -fPIC $(OPT)
# FFTW_INC = -I/opt/cray/pe/fftw/3.3.6.3/haswell/include
# linking MKL library
MKL_DIR=/opt/intel/mkl
MKL_INC              = $(MKL_DIR)/include
MKL_LIB              = $(MKL_DIR)/lib/intel64
LIBMKL                                                  = -L$(MKL_LIB) -lmkl_rt

all : $(TEMPLATES_OBJ)/templates.o

$(TEMPLATES_OBJ)/templates.o : templates.c
	$(CC) -c $< -I$(MKL_INC) $(LIBMKL) -o $@

clean :
	rm $(TEMPLATES_OBJ)/*.o
	rm -r $(TEMPLATES_OBJ)
