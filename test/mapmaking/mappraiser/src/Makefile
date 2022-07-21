#DIR = $(HOME)/midapack
DIR = $(MIDADIR)
MIDAPACK_ROOT=$(PREFIX)/midapack
MIDAPACK_LIB=$(MIDAPACK_ROOT)/lib
MPICC = cc
OPT =  -D W_MPI -D W_OPENMP -O3
CC = $(MPICC) -qopenmp -g -Wall -fPIC $(OPT)
INC = -I$(DIR)/include
LIB = -L$(MIDAPACK_LIB) -lmidapack
MAPPRAISER_ROOT=$(PREFIX)/mappraiser
MAPPRAISER_OBJ=$(MAPPRAISER_ROOT)/obj
MAPPRAISER_LIB=$(MAPPRAISER_ROOT)/lib

# Linking preAlps library for ECG solver utilities
CPLM_CORE     = utils/cplm_core
CPLMDIR      = utils/cplm_light
CPLM_V0_DIR  = utils/cplm_v0
preAlpsDIR = $(HOME)/preAlps

LIBpreAlps = -L$(preAlpsDIR)/lib/ -lpreAlps_release -lcplm_light -lcplm_v0 -lcplm_core
INCpreAlps = -I$(preAlpsDIR)/$(CPLM_CORE) -I$(preAlpsDIR)/$(CPLMDIR) -I$(preAlpsDIR)/$(CPLM_V0_DIR) -I$(preAlpsDIR)/utils -I$(preAlpsDIR)/src/preconditioners -I$(preAlpsDIR)/src/solvers

# Linking metis library needed by preAlps
METIS_DIR=/global/common/cori/software/metis/5.1.0
METIS_INC            = $(METIS_DIR)/include
METIS_LIB            = $(METIS_DIR)/lib
INCpreAlps += $(addprefix -I,$(METIS_INC))
LIBpreAlps += -L$(METIS_LIB) -lmetis

# linking MKL library
MKL_DIR=/opt/intel/mkl
MKL_INC              = $(MKL_DIR)/include
MKL_LIB              = $(MKL_DIR)/lib/intel64
INCpreAlps       += $(addprefix -I,$(MKL_INC))
LIBpreAlps            += -L$(MKL_LIB) -lmkl_rt

#supplementary libraries
HEALPIX_INC = -I$(HOME)/Healpix_3.50/include
CFITSIO_INC = -I/usr/common/software/cfitsio/3.47/include

HEALPIX_LIB = -L$(HOME)/Healpix_3.50/lib -lchealpix -lhealpix
CFITSIO_LIB = -L/usr/common/software/cfitsio/3.47/lib -lcfitsio

##############   t e s t s   ############

all:
	make mappraiser

mappraiser:
	@echo "start compiling MAPPRAISER ........"
	make $(MAPPRAISER_OBJ)/createToeplitz.o \
	$(MAPPRAISER_OBJ)/precond.o \
	$(MAPPRAISER_OBJ)/pcg_true.o \
	$(MAPPRAISER_OBJ)/ecg.o \
	$(MAPPRAISER_OBJ)/iofiles.o \
	$(MAPPRAISER_OBJ)/mappraiser.o \

$(MAPPRAISER_OBJ)/createToeplitz.o: createToeplitz.c
	$(CC) -c $< $(INC) $(LIB) -o $@

$(MAPPRAISER_OBJ)/precond.o: precond.c
	$(CC) -c $< $(INC) $(LIB) -o $@

$(MAPPRAISER_OBJ)/pcg_true.o: pcg_true.c
	$(CC) -c $< $(INC) $(LIB) $(CFITSIO_INC) $(CFITSIO_LIB) -o $@

$(MAPPRAISER_OBJ)/ecg.o: ecg.c
	$(CC) -c $< $(INC) $(LIB) $(INCpreAlps) $(LIBpreAlps) -o $@

$(MAPPRAISER_OBJ)/iofiles.o: iofiles.c
	$(CC) -c $< $(INC) $(LIB) $(CFITSIO_INC) $(CFITSIO_LIB) -o $@

$(MAPPRAISER_OBJ)/mappraiser.o: mappraiser.c
	$(CC) -c $< $(INC) $(LIB) $(CFITSIO_INC) $(INCpreAlps) $(LIBpreAlps) -o $@

# s4cmb_pipeline: createToeplitz.o precond.o s4cmb_pipeline.c iofiles.o pcg_true.o
# 	$(CC) s4cmb_pipeline.c createToeplitz.o precond.o pcg_true.o iofiles.o $(INC) $(LIB) $(FFTW_LIB) -lm -o s4cmb_pcg -Wl,@/global/homes/t/tcimic/midapack/test/mapmaking/allinea-profiler.ld

##############   u t i l s   #############

clean:
	rm -f $(MAPPRAISER_OBJ)/*.o ;
