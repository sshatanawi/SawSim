# Makefile for compiling PETSc-based SawSim groundwater model
# Output: SawSim_GW.exe

#  The following variable must either be a path to petsc.pc or just "petsc" if petsc.pc
#  has been installed to a system location or can be found in PKG_CONFIG_PATH.
# PETSc config (export PETSC_DIR and PETSC_ARCH before building)
petsc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/petsc.pc
# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(petsc.pc)

FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES)) -cpp
CPPFLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm

# Object files
OBJS = SawSim_declaration.o SawSim_subroutines.o SawSim_main.o

# Final executable
SawSim_GW.exe: $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.F90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $<


clean:
	rm -f *.o *.mod SawSim_GW.exe

# Print PETSc build variables (optional debug)
print:
	@echo FC=$(FC)
	@echo FFLAGS=$(FFLAGS)
	@echo CPPFLAGS=$(CPPFLAGS)
	@echo LDFLAGS=$(LDFLAGS)
	@echo LDLIBS=$(LDLIBS)

# Calling libraries in some of subroutines should be in a form of #include <> at the beginning of the line
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscksp.h>
#		use petsc
#		use petscvec
#		use petscsys
#		use petscksp

