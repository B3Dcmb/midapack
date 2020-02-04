DIR = /global/homes/a/ajamal/CMB_TOAST/midapack
MAPMAT = $(DIR)/src/mapmat/
F90 = ftn
GCC = gcc
# MPICC = /global/homes/e/elbouha/.conda/envs/hamzaenv/bin/mpicc
MPICC = cc
OPT =  -D W_MPI -D W_OPENMP -D MPI -D OPENMP -O3
CC = $(MPICC) -qopenmp $(OPT)
INC = -I$(DIR)/include -I/opt/cray/pe/fftw/3.3.8.3/haswell/include
LIB = -L$(DIR)/lib -lmidapack
#supplementary libraries
# PTSCOTCH_INC = -I$(HOME)/scotch_5.1.12/include/
HEALPIX_INC = -I$(HOME)/Healpix_3.50/include
CFITSIO_INC = -I/usr/common/software/cfitsio/3.47/include

# PTSCOTCH_LIB = -L$(HOME)/scotch_5.1.12/lib/ -lptscotch -lptscotcherr
HEALPIX_LIB = -L$(HOME)/Healpix_3.50/lib -lchealpix -lhealpix
CFITSIO_LIB = -L/usr/common/software/cfitsio/3.47/lib -lcfitsio
# HDF5_LIB = -L/usr/lib/ -lhdf5
# -lm -L/usr/lib64/gcc/x86_64-suse-linux/4.8/
INTEL_DIR = /opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64/
INTEL_LIB = $(INTEL_DIR)libsvml.a $(INTEL_DIR)libirc.a $(INTEL_DIR)libimf.a $(INTEL_DIR)libifcore.a $(INTEL_DIR)libiomp5.a
##############################################
all:
	make test
	make example

algo.o: algo.c
	$(CC) -c algo.c

io_hdf5.o: io_hdf5.c
	$(CC) -c io_hdf5.o

cg.o: cg.c
	$(CC) $(INC) -c cg.c


##############################################
test: $(MAPMAT)
	@echo "start compiling tests ........"
	make test_als
	make test_alm
	make test_cindex
	make test_csort
	make test_ring
	make test_butterfly


test_als: $(MAPMAT)
	$(CC) test_als_0.c $(INC) $(LIB) -o als_0
	$(CC) test_als_1.c $(INC) $(LIB) -o als_1
	$(CC) test_als_2.c $(INC) $(LIB) -o als_2
	$(CC) test_als_3.c $(INC) $(LIB) -o als_3

test_alm: $(MAPMAT)
	$(CC) test_alm_0.c $(INC) $(LIB) -o alm_0

test_cindex: $(MAPMAT)
	$(CC) test_cindex_0.c $(INC) $(LIB) -o cindex_0

test_csort: $(MAPMAT) test_csort_0.c test_csort_1.c
	$(CC) test_csort_0.c $(INC) $(LIB) -o csort_0
	$(CC) test_csort_1.c $(INC) $(LIB) -o csort_1

test_ring: $(MAPMAT)
	$(CC) test_ring_0.c $(INC) $(LIB) -o ring_0

test_butterfly: $(MAPMAT)
	$(CC) test_butterfly_0.c $(INC) $(LIB) -o butterfly_0
	$(CC) test_butterfly_1.c $(INC) $(LIB) -o butterfly_1
	$(CC) test_butterfly_2.c $(INC) $(LIB) -o butterfly_2

##############################################
example: $(MAPMAT)
	make sky_scan_loop
	make sky_scan
	# make test_init_mapmat
	# make test_cg_mapmat
	# make test_init_mapmatc
	# make test_cg_mapmatc

test_init_mapmat: $(MAPMAT) fem1D.c algo.o multiI.c test_spI_mapmat.c
	$(CC) fem1D.c $(INC) $(LIB) algo.o -o fem1D
	$(CC) multiI.c $(INC) $(LIB) algo.o -o multiI
	$(CC) -o test_spI_mapmat test_spI_mapmat.c $(INC) $(LIB) algo.o -o spI_mat

test_cg_mapmat: $(MAPMAT) test_rand_mapmat.c cg.o
	$(CC) test_rand_mat_0.c cg.o $(INC) $(LIB) -o rand_mat_0
	$(CC) test_diags_mat_0.c cg.o $(INC) $(LIB)  -o diags_mat_0

test_init_mapmatc: $(MAPMAT) algo.o
	@echo 'nothing to do'

test_cg_mapmatc: $(MAPMAT) test_rand_mapmat.c cg.o
	$(CC) test_rand_matc.c cg.o $(INC) $(LIB) -o rand_matc
	$(CC) test_diags_matc.c cg.o $(INC) $(LIB)  -o diags_matc

test_spI_mapmat: $(MAPMAT)mapmat.o test_spI_mapmat.c

sky_scan : sky_scan.c $(MAPMAT) algo.o map2gif.o
	$(CC) sky_scan.c -o sky_scan $(HEALPIX_INC) $(INC) $(LIB) algo.o map2gif.o -lintlc -Wl,--start-group $(INTEL_LIB) $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB) -Wl,--end-group
sky_scan_loop : sky_scan_loop.c $(MAPMAT) algo.o map2gif.o
	$(CC) sky_scan_loop.c -o sky_scan_loop $(HEALPIX_INC) $(INC) $(LIB) algo.o map2gif.o -lintlc -Wl,--start-group $(INTEL_LIB) $(HEALPIX_LIB) -lhpxgif -lgif $(CFITSIO_LIB) -Wl,--end-group

map2gif.o : map2gif.f90
	$(F90) -c map2gif.f90 -o map2gif.o $(HEALPIX_INC)




clean :
	rm -f *.o ;
