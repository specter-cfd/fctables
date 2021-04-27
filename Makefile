#****************************************************************************
# Makefile for compiling codes and linking with MPI, FFTP and FFTW libraries
# Mauro Fontana - 2/4/2019
#***************************************************************************

#****************************************************************************
# Edit library paths, compiler flags, and code options in Makefile.in
#****************************************************************************
include Makefile.in

#****************************************************************************
# Don't edit below this line
#****************************************************************************
CPP         := $(CPP_$(COMPILER))
FC          := $(FC_$(COMPILER))
CC          := $(CC_$(COMPILER))
FFLAGS      := $(FFLAGS_$(COMPILER))
TARGET_ARCH := $(TARGET_$(COMPILER))

# OpenMP-MPI  hybdridization
DO_HYBRIDyes:= $($(COMPILER)_OMP)
DO_HYBRIDno :=
LOPENMP     := $(DO_HYBRID$(OPENMP))

BLENDyes:= -DBLEND_
BLENDno :=

GMP_BUILD_DIR :=contrib/gmp/build
MPFR_BUILD_DIR:=contrib/mpfr/build

MPFUN_LINK_DEPS_fort:= 
MPFUN_LINK_DEPS_mpfr:= $(MPFR_BUILD_DIR)/lib/libmpfr.a
MPFUN_LINK_DEPS_mpfr:= $(MPFUN_LINK_DEPS_mpfr) $(GMP_BUILD_DIR)/lib/libgmp.a
MPFUN_LINK_DEPS     := $(MPFUN_LINK_DEPS_$(MPFUN_VERSION))

MPFUN_MAKE   := Makefile-$(MPFUN_VERSION)
MPFUN_BUILD  := contrib/mpfun/build-$(MPFUN_VERSION)
MPFUN_INCLUDE:= -I$(MPFUN_BUILD)

export FC CC

all: fc_tables

edit:
	$(CPP) -D_DIGITS=$(DIGITS) $(BLEND$(BLEND)) fc_tables.fpp \
		-o fc_tables.f90 

fc_tables: mpfun edit
	$(FC) $(FFLAGS) $(LOPENMP) mprlinalg_mod.f90 fc_tables.f90 \
		$(MPFUN_BUILD)/*.o $(MPFUN_LINK_DEPS)\
		$(MPFUN_INCLUDE) -o fc_tables

mpfun:
	$(MAKE) -C contrib/mpfun -f $(MPFUN_MAKE) COMPILER=$(COMPILER)\
		DIGITS=$(DIGITS)

test:
	$(FC) $(FFLAGS) fc_test.f08 -o fc_test

clean:
	rm -f *.mod fc_tables fc_test

distclean:
	$(MAKE) distclean -C contrib/mpfun -f $(MPFUN_MAKE)
	rm -f *.mod fc_tables fc_tables.f90 fc_test
