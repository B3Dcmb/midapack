GCC = gcc
# MPICC = /global/homes/e/elbouha/.conda/envs/hamzaenv/bin/mpicc
MPICC = cc
OPT =  -D W_MPI -D W_OPENMP -O3
CC = $(MPICC) -qopenmp -g -Wall -Werror -fPIC $(OPT)
FFTW_INC = -I/opt/cray/pe/fftw/3.3.6.3/haswell/include


all : mapmat.o mapmatc.o

bitop.o : bitop.c
	$(CC) -c bitop.c

als.o : als.c
	$(CC) -c als.c

alm.o : alm.c als.o
	$(CC) -c alm.c

cindex.o : cindex.c
	$(CC) -c cindex.c

csort.o : csort.c als.o
	$(CC)  -c csort.c

ring.o : ring.c alm.o
	$(CC) -c ring.c

butterfly.o : butterfly.c alm.o
	$(CC) -c butterfly.c

mapmat.o : mapmat.c cindex.o csort.o bitop.o ring.o butterfly.o
	$(CC) -c mapmat.c

mapmatc.o : mapmatc.c cindex.o csort.o bitop.o ring.o butterfly.o
	$(CC) -c mapmatc.c

clean :
	rm *.o
