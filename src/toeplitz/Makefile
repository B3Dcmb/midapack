MIDAPACK_ROOT=$(PREFIX)/midapack_2lvl
MIDAPACK_OBJ=$(MIDAPACK_ROOT)/obj
TOEPLITZ_OBJ=$(MIDAPACK_OBJ)/toeplitz

GCC = gcc
MPICC = cc
OPT =  -D W_MPI -D W_OPENMP -O3
CC = $(MPICC) -qopenmp -g -Wall -Werror -fPIC $(OPT)

##############   s r c   ############
all:
	make $(TOEPLITZ_OBJ)/toeplitz.o \
	$(TOEPLITZ_OBJ)/toeplitz_seq.o \
	$(TOEPLITZ_OBJ)/toeplitz_block.o \
	$(TOEPLITZ_OBJ)/toeplitz_gappy.o \
	$(TOEPLITZ_OBJ)/toeplitz_nofft.o \
	$(TOEPLITZ_OBJ)/toeplitz_rshp.o \
	$(TOEPLITZ_OBJ)/toeplitz_params.o \
	$(TOEPLITZ_OBJ)/toeplitz_utils.o \
	$(TOEPLITZ_OBJ)/toeplitz_wizard.o

$(TOEPLITZ_OBJ)/toeplitz.o: toeplitz.c toeplitz_nofft.c
	$(CC) -c $< -o $@

$(TOEPLITZ_OBJ)/toeplitz_seq.o: toeplitz_seq.c
	$(CC) -c $< -o $@

$(TOEPLITZ_OBJ)/toeplitz_block.o: toeplitz_block.c
	$(CC) -c $< -o $@

$(TOEPLITZ_OBJ)/toeplitz_gappy.o: toeplitz_gappy.c
	$(CC) -c $< -o $@

$(TOEPLITZ_OBJ)/toeplitz_nofft.o: toeplitz_nofft.c
	$(CC) -c $< -o $@

$(TOEPLITZ_OBJ)/toeplitz_rshp.o: toeplitz_rshp.c
	$(CC) -c $< -o $@

$(TOEPLITZ_OBJ)/toeplitz_params.o: toeplitz_params.c
	$(CC) -c $< -o $@

$(TOEPLITZ_OBJ)/toeplitz_utils.o: toeplitz_utils.c
	$(CC) -c $< -o $@

$(TOEPLITZ_OBJ)/toeplitz_wizard.o: toeplitz_wizard.c
	$(CC) -c $< -o $@

########   u t i l s   ########

clean:
	rm -f $(TOEPLITZ_OBJ)/*.o
	rm -r $(TOEPLITZ_OBJ)
