DIR = /global/homes/e/elbouha/midapack
TOEPLITZ = $(DIR)/src/toeplitz/
MPICC = cc
CC = $(MPICC) -qopenmp -g -Wall -fPIC $(OPT)
INC = -I$(DIR)/include
LIB = -L$(DIR)/lib -lmidapack

CPLM_CORE     = utils/cplm_core
CPLMDIR      = utils/cplm_light
CPLM_V0_DIR  = utils/cplm_v0
preAlpsDIR = /global/homes/e/elbouha/preAlps

LIBpreAlps = -L/global/homes/e/elbouha/preAlps/lib/ -lpreAlps_release -lcplm_light -lcplm_v0 -lcplm_core
INCpreAlps = -I$(preAlpsDIR) -I../ -I../include -I$(preAlpsDIR)/$(CPLM_CORE) -I$(preAlpsDIR)/$(CPLMDIR) -I$(preAlpsDIR)/$(CPLM_V0_DIR) -I$(preAlpsDIR)/utils -I$(preAlpsDIR)/src/preconditioners -I$(preAlpsDIR)/src/solvers

METIS_DIR=/global/common/cori/software/metis/5.1.0
METIS_INC            = $(METIS_DIR)/include
METIS_LIB            = $(METIS_DIR)/lib
INCpreAlps += $(addprefix -I,$(METIS_INC))
LIBpreAlps += -L$(METIS_LIB) -lmetis

MKL_DIR=/opt/intel/mkl
MKL_INC              = $(MKL_DIR)/include
MKL_LIB              = $(MKL_DIR)/lib/intel64
INCpreAlps       += $(addprefix -I,$(MKL_INC))
LIBpreAlps            += -L$(MKL_LIB) -lmkl_rt

# FFTW_LIB = -L/opt/cray/pe/fftw/3.3.6.3/haswell/lib -lfftw3 -lfftw3_threads
# FFTW_INC = -I/opt/cray/pe/fftw/3.3.6.3/haswell/include
# FFTW_LIB = -L/opt/cray/pe/fftw/3.3.6.3/haswell/lib/libfftw3.a /opt/cray/pe/fftw/3.3.6.3/haswell/lib/libfftw3_threads.a

SRC = $(DIR)/src/toeplitz/

MAPMAT = $(DIR)/src/mapmat/
F90 = ftn
OPT =  -D W_MPI -D W_OPENMP -O3

#supplementary libraries
HEALPIX_INC = -I$(HOME)/Healpix_3.50/include
CFITSIO_INC = -I/usr/common/software/cfitsio/3.47/include

HEALPIX_LIB = -L$(HOME)/Healpix_3.50/lib -lchealpix -lhealpix
CFITSIO_LIB = -L/usr/common/software/cfitsio/3.47/lib -lcfitsio

##############   t e s t s   ############

all:
	make mappraiser

mappraiser: $(TOEPLITZ) $(MAPMAT)
	@echo "start compiling MAPPRAISER ........"
	make createToeplitz.o \
	precond.o \
	pcg_true.o \
	ecg.o \
	iofiles.o \
	mappraiser.o \

createToeplitz.o: $(SRC)
	$(CC) -c createToeplitz.c $(INC) $(LIB)

precond.o: $(MAPMAT)
	$(CC) -c precond.c $(INC) $(LIB)

pcg_true.o: $(SRC) $(MAPMAT)
	$(CC) -c pcg_true.c $(INC) $(LIB)

ecg.o:
	$(CC) -c ecg.c $(INC) $(LIB) $(INCpreAlps) $(LIBpreAlps)

iofiles.o:
	$(CC) -c iofiles.c $(INC) $(LIB) $(CFITSIO_INC) $(CFITSIO_LIB)

mappraiser.o: $(SRC) $(MAPMAT)
	$(CC) -v -c mappraiser.c $(INC) $(LIB) $(CFITSIO_INC) $(INCpreAlps) $(LIBpreAlps)

# s4cmb_pipeline: createToeplitz.o precond.o s4cmb_pipeline.c iofiles.o pcg_true.o
# 	$(CC) s4cmb_pipeline.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -lm -o s4cmb_pcg -Wl,@/global/homes/e/elbouha/midapack/test/mapmaking/allinea-profiler.ld

##############   u t i l s   #############

clean:
	rm -f *.o ;
