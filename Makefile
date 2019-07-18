PROJECT = Midapack_
VERSION = 1.1b
DIR = /global/homes/e/elbouha/midapack
DIRTAR = export_tar
LIBNAME = libmidapack

all :
	@echo "starts compiling ........"
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


lib: ./src/mapmat/ ./src/toeplitz
	ar r ./lib/$(LIBNAME).a src/mapmat/mapmat.o src/mapmat/mapmatc.o src/mapmat/bitop.o src/mapmat/als.o src/mapmat/alm.o \
           src/mapmat/csort.o src/mapmat/cindex.o src/mapmat/ring.o src/mapmat/butterfly.o \
           src/toeplitz/toeplitz.o src/toeplitz/toeplitz_seq.o src/toeplitz/toeplitz_block.o \
           src/toeplitz/toeplitz_nofft.o src/toeplitz/toeplitz_gappy.o src/toeplitz/toeplitz_params.o \
           src/toeplitz/toeplitz_rshp.o src/toeplitz/toeplitz_utils.o src/toeplitz/toeplitz_wizard.o
	ranlib ./lib/$(LIBNAME).a
	cc -shared src/mapmat/mapmat.o src/mapmat/mapmatc.o src/mapmat/bitop.o src/mapmat/als.o src/mapmat/alm.o \
           src/mapmat/csort.o src/mapmat/cindex.o src/mapmat/ring.o src/mapmat/butterfly.o \
           src/toeplitz/toeplitz.o src/toeplitz/toeplitz_seq.o src/toeplitz/toeplitz_block.o \
           src/toeplitz/toeplitz_nofft.o src/toeplitz/toeplitz_gappy.o src/toeplitz/toeplitz_params.o \
           src/toeplitz/toeplitz_rshp.o src/toeplitz/toeplitz_utils.o src/toeplitz/toeplitz_wizard.o \
					 -L/opt/cray/pe/fftw/3.3.6.3/haswell/lib -lfftw3 -lfftw3_threads -L/opt/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64 -liomp5 -o ./lib/$(LIBNAME).so

seq :
	make mapmat
	make toeplitz
	make lib

clean:
	make clean -C ./src/toeplitz
	make clean -C ./src/mapmat
	make clean -C ./test/toeplitz
	make clean -C ./test/mapmat
	rm lib/*midapack.a

doc : ./src ./src/mapmat ./src/toeplitz ./doc/mkdoc.dox
	pushd ./doc ; doxygen mkdoc.dox ; popd

tar :
	#rm -rf $(DIRTAR)
	#mkdir $(DIRTAR)
	tar -czvf $(PROJECT)$(VERSION).tgz ./* --directory=newdir --exclude '*.svn*' --exclude './src/lens2hat_repo/*'
	#mv $(PROJECT)$(VERSION).tgz $(DIRTAR)

release :
	./mkrelease
