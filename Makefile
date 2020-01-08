# ============================================================================
# Name        : Makefile
# Author      : Pierre-Carl Michaud
# Version     : 32.0
# Copyright   : Your copyright notice
# Description : Makefile pour Fonseca et al. JEEA
# Software used: Stata, Python (Anaconda) and Fortran (OpenMPI and gfortran)
# check README file for more details
# ============================================================================

.PHONY: all clean

# Change this line if you are using a different Fortran compiler (check flags, we are using openmpi with gfortran)
FORTRAN_COMPILER = mpif90
FFLAGS = -O3 -funroll-all-loops -finit-local-zero -mfpmath=sse -fstrength-reduce -fexpensive-optimizations  -fbounds-check   #-march=native
OBJ = obj

# this is to compile the statistical library
MAIN0 = libs/dcdflib.a

# this is the executable for estimation
MAIN1   =  ./runtime/estimate
SOURCE1 = $(OBJ)/newuoa.o $(OBJ)/core36.o   src/runestimation.f90 

# this is the executable to generate draws for MSM
MAIN2 = ./runtime/draws
SOURCE2 = $(OBJ)/newuoa.o $(OBJ)/core36.o  src/draws.f90

# thi is the executable to generate executable for simulation
MAIN3   =  ./runtime/generate
SOURCE3 = $(OBJ)/newuoa.o $(OBJ)/core36.o   src/generate.f90 

model: $(MAIN0) $(MAIN1) $(MAIN2) $(MAIN3) 
	@echo "Model compiled successfully"
input:
	stata-mp < do/master-dataprep.do > params/input/control.log
	@echo "*** Prepared inputs from data"
all : $(MAIN0)  $(MAIN1) $(MAIN2) $(MAIN3) 
	stata-mp < do/master-dataprep.do > params/input/control.log
	@echo "*** Sucessfully compiled model and prepared data"

$(MAIN0):
	gfortran -c libs/dcdflib.f/src/*.f
	ar r libs/dcdflib.a *.o
	rm *.o 

$(MAIN1): $(SOURCE1)
	$(FORTRAN_COMPILER) $(SOURCE1)  -J$(OBJ) libs/dcdflib.a    -o ./runtime/estimate  $(FFLAGS)  

$(MAIN2): $(SOURCE2)
	$(FORTRAN_COMPILER) $(SOURCE2)  -J$(OBJ) libs/dcdflib.a   -o ./runtime/draws $(FFLAGS)

$(MAIN3): $(SOURCE3)
	$(FORTRAN_COMPILER) $(SOURCE3)  -J$(OBJ) libs/dcdflib.a  -o ./runtime/generate $(FFLAGS)

$(OBJ)/core36.o: src/core36.f90 
	$(FORTRAN_COMPILER)  -J$(OBJ)  -c $< -o $@  $(FFLAGS) #-g

$(OBJ)/newuoa.o: src/newuoa.f90 
	gfortran  -J$(OBJ)  -c $< -o $@  $(FFLAGS) #-g


clean:
	rm -f obj/* libs/dcdflib.a
