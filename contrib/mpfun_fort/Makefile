FFLAGS  = -O3

OBJS    = mpfuna.mod mpfunb.mod mpfunc.mod mpfund.mod mpfune.mod\
		 	 mpfunf.mod mpfung.mod second.o

REGEX   = s/parameter (mpipl = [0-9]\+)/parameter (mpipl = $(DIGITS))/

mpmodule.mod : $(OBJS)
	$(FC) $(FFLAGS) -c mpmodule.f90


# Explicit treatment of modules whose name differs from filename
mpfunb.mod : mpfunbq.f90
	$(FC) $(FFLAGS) -c mpfunbq.f90

# mpfung depends on mpfunf which gets recompiled always
mpfung.mod : mpfunf.mod mpfungq1.f90
	$(FC) $(FFLAGS) -c mpfungq1.f90

# Change precision
mpfunf.mod : mpfunf.f90
    sed -i '$(REGEX)' mpfunf.f90
	$(FC) $(FFLAGS) -c mpfunf.f90

%.mod: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

updatedigits :
	
clean :
	rm *.o *.mod

distclean : clean

.PHONY : mpfunf.mod
