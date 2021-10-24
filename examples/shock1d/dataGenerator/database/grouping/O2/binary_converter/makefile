#Fortran compiler given below
#FC=gfortran

#Compiler flags
FCFLAGS=-ffree-form -pipe -std=gnu -g -pg -O0 -Wall -fconvert=big-endian -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -Wsurprising -Wtabs -Wunderflow -Wintrinsic-shadow -ggdb -fdump-core -fno-unsafe-math-optimizations -finit-integer=n -finit-real=nan
#FCFLAGS=-O3 -fno-unsafe-math-optimizations -fconvert=big-endian

#Variables used to include Plato header/module files and link the library.
#Variables defining the directories storing header/module files and the shared library. These 
#are contained in the installation directory PLATO_INSTALL_DIR

default: all

convert.o: convert.F90
	$(FC) $(INC) $(FCFLAGS) -c convert.F90

convert: convert.o
	$(FC) $(INC) $(FCFLAGS) -o convert_ convert.o 

all: convert

clean:
	@echo 'Cleaning...'
	@rm -f *.o *_ *.mod
	@rm -f *.dat; rm -f species; rm -f components; rm -f elements; rm -f *.out
	@echo 'Done with clean!'
