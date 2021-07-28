MIDAPACK_ROOT=$(PREFIX)/midapack_bis
MIDAPACK_OBJ=$(MIDAPACK_ROOT)/obj
MAPMAT_OBJ=$(MIDAPACK_OBJ)/mapmat

GCC = gcc
# MPICC = /global/homes/e/elbouha/.conda/envs/hamzaenv/bin/mpicc
MPICC = cc
OPT =  -D W_MPI -D W_OPENMP -O3
CC = $(MPICC) -qopenmp -g -Wall -Werror -fPIC $(OPT)
# FFTW_INC = -I/opt/cray/pe/fftw/3.3.6.3/haswell/include

all : $(MAPMAT_OBJ)/mapmat.o $(MAPMAT_OBJ)/mapmatc.o

$(MAPMAT_OBJ)/bitop.o : bitop.c
	$(CC) -c $< -o $@

$(MAPMAT_OBJ)/als.o : als.c
	$(CC) -c $< -o $@

$(MAPMAT_OBJ)/alm.o : alm.c $(MAPMAT_OBJ)/als.o
	$(CC) -c $< -o $@

$(MAPMAT_OBJ)/cindex.o : cindex.c
	$(CC) -c $< -o $@

$(MAPMAT_OBJ)/csort.o : csort.c $(MAPMAT_OBJ)/als.o
	$(CC) -c $< -o $@

$(MAPMAT_OBJ)/ring.o : ring.c $(MAPMAT_OBJ)/alm.o
	$(CC) -c $< -o $@

$(MAPMAT_OBJ)/butterfly.o : butterfly.c $(MAPMAT_OBJ)/alm.o
	$(CC) -c $< -o $@

$(MAPMAT_OBJ)/mapmat.o : mapmat.c $(MAPMAT_OBJ)/cindex.o $(MAPMAT_OBJ)/csort.o $(MAPMAT_OBJ)/bitop.o $(MAPMAT_OBJ)/ring.o $(MAPMAT_OBJ)/butterfly.o
	$(CC) -c $< -o $@

$(MAPMAT_OBJ)/mapmatc.o : mapmatc.c $(MAPMAT_OBJ)/cindex.o $(MAPMAT_OBJ)/csort.o $(MAPMAT_OBJ)/bitop.o $(MAPMAT_OBJ)/ring.o $(MAPMAT_OBJ)/butterfly.o
	$(CC) -c $< -o $@

clean :
	rm $(MAPMAT_OBJ)/*.o
	rm -r $(MAPMAT_OBJ)
