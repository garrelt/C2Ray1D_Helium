# MAKROS----------------------------------------------------
# F90 compiler

#F90 = ifort
F90 = gfortran

# ifort options

#OPTIONS = -xW -O3 -vec_report -u -fpe0 -check all \
#-traceback -ipo -DIFORT

#OPTIONS = -pg -g  -O0 -u -fpe0 -check all \
#-traceback -ipo -DIFORT

GFORTFLAGS = -O
F90FLAGS1 = $(GFORTFLAGS) 
F90FLAGS = $(F90FLAGS1)
OPTIONS = $(F90FLAGS)

UTILS1=romberg.o string.o mrgrnk.o  ctrper.o

UTILS2= romberg.o string.o mrgrnk.o  ctrper.o

CONSTANTS1 = mathconstants.o cgsconstants.o  cgsphotoconstants.o \
 cgsastroconstants.o c2ray_parameters.o cosmoparms.o abundances.o \
atomic.o cosmoparms.o

CONSTANTS2 = mathconstants.o cgsconstants.o  cgsphotoconstants.o \
 cgsastroconstants.o c2ray_parameters.o cosmoparms.o abundances.o atomic.o \
cosmoparms.o

# there are two different versions of doric, one which solves the
# inhomogeneous system of 3 ode (doric2_july2010) and one which
# solves the homogeneous system of 5 ode (doric3b)

#DORIC= doric3b.o 
DORIC= doric.o # 

DORIC2= doric.o # 

# there are also some different versions of RADIATION.
# radiation_monocromatic takes different inputs, see input file
# inputs/isochromatic
RADIATION= radiation.o 
#RADIATION= radiation_monocromatic.o 
#RADIATION= ../radiation.o
#RADIATION2= radiation_monocromatic.o
RADIATION2= radiation.o
#RADIATION2= radiation.o
#-----------------------------------------------------------

# Building C2Ray_1D: 
# $@ means: name of the file to be made 

a: precision.o sizes.o $(CONSTANTS1) $(UTILS1) file_admin.o \
 no_mpi.o clocks.o grid.o tped.o  sourceprops_test_one_source.o  material.o cosmology.o\
cooling_h.o $(RADIATION) thermal.o time.o timeequation.o $(DORIC) \
photonstatistics.o evolve.o LTE.o timestep.o output.o C2RayA.o


	$(F90) $(OPTIONS) -o $@ precision.o $(CONSTANTS2) $(UTILS2) \
	file_admin.o sizes.o no_mpi.o  clocks.o grid.o \
	tped.o  sourceprops_test_one_source.o  cosmology.o  material.o cooling_h.o \
	$(RADIATION2) thermal.o time.o timeequation.o $(DORIC2) photonstatistics.o \
	evolve.o LTE.o timestep.o output.o C2RayA.o


clean : 
	rm -f *.o *.mod 

# Building object files: 
# $< means: name of related file that caused the action

.f90.o:
	$(F90) -c $(OPTIONS) $<
.F90.o:
	$(F90) -c $(OPTIONS) $<
f.mod:
	$(F90) -c $(OPTIONS) $<

# Suffix rules: List of significant suffixes
.SUFFIXES: .f90 .F90 .mod .o
