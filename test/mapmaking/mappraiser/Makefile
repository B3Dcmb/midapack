DIR = /global/homes/t/tcimic/midapack
TOEPLITZ = $(DIR)/src/toeplitz/
MPICC = cc
CC = $(MPICC) -qopenmp -g -Wall -fPIC $(OPT)
INC = -I$(DIR)/include
LIB = -L$(DIR)/lib -lmidapack

CPLM_CORE     = utils/cplm_core
CPLMDIR      = utils/cplm_light
CPLM_V0_DIR  = utils/cplm_v0
preAlpsDIR = /global/homes/t/tcimic/preAlps

LIBpreAlps = -L/global/homes/t/tcimic/preAlps/lib/ -lpreAlps_release -lcplm_light -lcplm_v0 -lcplm_core
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
	make map
	make library

map: ./src/
	make -C ./src/

library:
	cc -shared ./src/createToeplitz.o ./src/precond.o ./src/pcg_true.o ./src/ecg.o ./src/iofiles.o ./src/mappraiser.o $(LIB) $(CFITSIO_LIB) $(LIBpreAlps) -o ./lib/libmappraiser.so

##############   u t i l s   #############

clean:
	make clean -C ./src/
	rm -f lib/*midapack.so
