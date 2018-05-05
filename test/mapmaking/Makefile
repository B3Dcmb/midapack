DIR = /home/dauvergn/midas-svn-OK/midas/trunk
TOEPLITZ = $(DIR)/src/toeplitz/
GCC = gcc
MPICC = /usr/local/openmpi/bin/mpicc
CC = $(MPICC) -fopenmp $(OPT)
INC = -I$(DIR)/include -I/usr/include
LIB = -L$(DIR)/lib -lmidapack
FFTW_LIB = -lfftw3 -lfftw3_threads
#FFTW     = -D FFTW_MULTITHREADING

SRC = ../../src/toeplitz/

MAPMAT = $(DIR)/src/mapmat/
F90 = gfortran
OPT =  -D W_MPI -D W_OPENMP -O3

#supplementary libraries
PTSCOTCH_INC = -I/home/cargemel/scotch_5.1.12/include/
HEALPIX_INC = -I/home/cargemel/Healpix_2.20a/include
CFITSIO_INC = -I/home/cargemel/cfitsio/include

#HEALPIX_INC = -I/home/dauvergn/midas-svn/midas/trunk/test/mapmaking/healpix/include

PTSCOTCH_LIB = -L/home/cargemel/scotch_5.1.12/lib/ -lptscotch -lptscotcherr
HEALPIX_LIB = -L/home/cargemel/Healpix_2.20a/lib -lchealpix -lhealpix
CFITSIO_LIB = -L/home/cargemel/cfitsio/lib -lcfitsio

#HEALPIX_LIB = -L/home/dauvergn/midas-svn/midas/trunk/test/mapmaking/healpix/lib -lchealpix -lhealpix


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
	test_pcg_rand

exemple_iofilepix:
	make map2gif.o \
	mk_iofilespix

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
	$(CC) mk_iofilespix.c -o mk_iofilespix $(INC) $(LIB) map2gif.o iofiles.o -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.4/ -lgfortran $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB) 

mk_iofilespix_diff : mk_iofilespix_diff.c map2gif.o iofiles.o
	$(CC) mk_iofilespix_diff.c -o mk_iofilespix_diff $(INC) $(LIB) map2gif.o iofiles.o -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.4/ -lgfortran $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB)

mk_iofilespix_diff2 : mk_iofilespix_diff2.c map2gif.o iofiles.o
	$(CC) mk_iofilespix_diff2.c -o mk_iofilespix_diff2 $(INC) $(LIB) map2gif.o iofiles.o -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.4/ -lgfortran $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB)

mk_unpointing : mk_unpointing.c map2gif.o iofiles.o
	$(CC) mk_unpointing.c -o mk_unpointing $(INC) $(LIB) map2gif.o iofiles.o -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.4/ -lgfortran $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB)
#sky_scan : sky_scan.c $(MAPMAT) algo.o map2gif.o
#	$(CC) sky_scan.c -o sky_scan $(INC) $(LIB) algo.o map2gif.o -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.4/ -lgfortran $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB) 

test_pcg: createToeplitz.o precond.o test_pcg.c iofiles.o pcg_true.o
	$(CC) test_pcg.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -o pcg


test_pcg_rand: createToeplitz.o precond.o test_pcg_rand.c iofiles.o pcg_true.o
	$(CC) test_pcg_rand.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -o pcg_rand


##############   u t i l s   #############

clean:
	rm -f *.o ;


