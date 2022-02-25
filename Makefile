PROJECT = Midapack_
VERSION = 1.1b
DIR = $(HOME)/midapack
DIRTAR = export_tar
LIBNAME = libmidapack
MIDAPACK_ROOT=$(PREFIX)/midapack
MIDAPACK_OBJ=$(MIDAPACK_ROOT)/obj
MAPMAT_OBJ=$(MIDAPACK_OBJ)/mapmat
TOEPLITZ_OBJ=$(MIDAPACK_OBJ)/toeplitz
MIDAPACK_LIB=$(MIDAPACK_ROOT)/lib

.PHONY: lib
all : 
	@echo "starts compiling ........"
	mkdir -p $(MIDAPACK_ROOT) $(MIDAPACK_OBJ) $(MAPMAT_OBJ) $(TOEPLITZ_OBJ) $(MIDAPACK_LIB)
	make mapmat
	make toeplitz 
	make lib

example :
	make mapmat_example
	make toeplitz_example

test :
	make mapmat_test
	make toeplitz_test

toeplitz: ./src/toeplitz/ 
	make -C ./src/toeplitz/

toeplitz_test: ./test/toeplitz/ ./src/toeplitz/
	make test -C ./test/toeplitz/

toeplitz_example: ./lib
	make example -C test/toeplitz/


mapmat: ./src/mapmat/ 
	make -C ./src/mapmat/

mapmat_test: ./test/mapmat/ ./src/mapmat/
	make test -C ./test/mapmat/

mapmat_example: ./lib
	make example -C test/mapmat/


lib:
	ar r $(MIDAPACK_LIB)/$(LIBNAME).a $(MAPMAT_OBJ)/mapmat.o $(MAPMAT_OBJ)/mapmatc.o \
		$(MAPMAT_OBJ)/bitop.o $(MAPMAT_OBJ)/als.o $(MAPMAT_OBJ)/alm.o \
        $(MAPMAT_OBJ)/csort.o $(MAPMAT_OBJ)/cindex.o $(MAPMAT_OBJ)/ring.o $(MAPMAT_OBJ)/butterfly.o \
        $(TOEPLITZ_OBJ)/toeplitz.o $(TOEPLITZ_OBJ)/toeplitz_seq.o $(TOEPLITZ_OBJ)/toeplitz_block.o \
        $(TOEPLITZ_OBJ)/toeplitz_nofft.o $(TOEPLITZ_OBJ)/toeplitz_gappy.o $(TOEPLITZ_OBJ)/toeplitz_params.o \
        $(TOEPLITZ_OBJ)/toeplitz_rshp.o $(TOEPLITZ_OBJ)/toeplitz_utils.o $(TOEPLITZ_OBJ)/toeplitz_wizard.o
	ranlib $(MIDAPACK_LIB)/$(LIBNAME).a
	cc -qopenmp -shared $(MAPMAT_OBJ)/mapmat.o $(MAPMAT_OBJ)/mapmatc.o $(MAPMAT_OBJ)/bitop.o \
	 	$(MAPMAT_OBJ)/als.o $(MAPMAT_OBJ)/alm.o $(MAPMAT_OBJ)/csort.o $(MAPMAT_OBJ)/cindex.o \
		$(MAPMAT_OBJ)/ring.o $(MAPMAT_OBJ)/butterfly.o $(TOEPLITZ_OBJ)/toeplitz.o $(TOEPLITZ_OBJ)/toeplitz_seq.o \
		$(TOEPLITZ_OBJ)/toeplitz_block.o $(TOEPLITZ_OBJ)/toeplitz_nofft.o $(TOEPLITZ_OBJ)/toeplitz_gappy.o \
		$(TOEPLITZ_OBJ)/toeplitz_params.o $(TOEPLITZ_OBJ)/toeplitz_rshp.o $(TOEPLITZ_OBJ)/toeplitz_utils.o \
		$(TOEPLITZ_OBJ)/toeplitz_wizard.o $(FFTW_LIB) -o $(MIDAPACK_LIB)/$(LIBNAME).so


seq :
	make mapmat
	make toeplitz 
	make lib


clean:
	make clean -C ./src/toeplitz
	make clean -C ./src/mapmat
	make clean -C ./test/toeplitz
	make clean -C ./test/mapmat
	rm $(MIDAPACK_LIB)/*midapack.*
	rm -r $(MIDAPACK_ROOT)

doc : ./src ./src/mapmat ./src/toeplitz ./doc/mkdoc.dox
	pushd ./doc ; doxygen mkdoc.dox ; popd

tar :
	#rm -rf $(DIRTAR)
	#mkdir $(DIRTAR)
	tar -czvf $(PROJECT)$(VERSION).tgz ./* --directory=newdir --exclude '*.svn*' --exclude './src/lens2hat_repo/*'
	#mv $(PROJECT)$(VERSION).tgz $(DIRTAR)

release :
	./mkrelease



