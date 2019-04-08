UNDERS_INTEL = -assume nounderscore
UNDERS_GNU   = -fno-underscoring

FFLAGS       = -O3 $(UNDERS_$(COMPILER))

OBJS         = mpfuna.mod mpfunf.mod mpfung.mod second.o mpinterface.o

REGEX        = s/parameter (mpipl = [0-9]\+)/parameter (mpipl = $(DIGITS))/


mpmodule.mod : gmp mpfr $(OBJS)
	$(FC) $(FFLAGS) -c mpmodule.f90

# Explicit treatment of modules whose name differs from filename.
# mpfung depends on mpfunf which is always rebuilt.
mpfung.mod : mpfunf.mod mpfung1.f90
	$(FC) $(FFLAGS) -c mpfung1.f90

# This one and updatedigits recipes are a hack to recompile only mpfun
# if DIGITS changes
mpfunf.mod : mpfunf.f90
	sed -i '$(REGEX)' mpfunf.f90
	$(FC) $(FFLAGS) -c mpfunf.f90

%.mod: %.f90
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

mpinterface.o: mpinterface.c
	$(CC) -c -I$$PWD/mpfr/build/include -I$$PWD/gmp/build/include \
		-o mpinterface.o mpinterface.c

# Point gmp and mpfr at their static libs so that they get rebuilt
# only after a distclean. Note that a make distclean is run after
# install because the test folders get quite big otherwise.
gmp/build/lib/libgmp.a :
	cd gmp && CC=$(CC) ./configure --prefix=$$PWD/build
	CC=$(CC) $(MAKE) -C gmp
	CC=$(CC) $(MAKE) check -C gmp
	CC=$(CC) $(MAKE) install -C gmp
	CC=$(CC) $(MAKE) distclean -C gmp

mpfr/build/lib/libmpfr.a :
	cd mpfr && CC=$(CC) ./configure --with-gmp-lib=../gmp/build/lib \
		\-\-with-gmp-include=../gmp/build/include --prefix=$$PWD/build
	CC=$(CC) $(MAKE) -C mpfr
	CC=$(CC) $(MAKE) check -C mpfr
	CC=$(CC) $(MAKE) install -C mpfr
	CC=$(CC) $(MAKE) distclean -C mpfr

gmp : gmp/build/lib/libgmp.a
mpfr : mpfr/build/lib/libmpfr.a

clean :
	rm *.o *.mod

cleanlibs: 
	cd gmp && rm -r build/
	cd mpfr && rm -r build/

# Also uninstalls mpfr and gmp
distclean: cleanlibs clean

.PHONY: gmp mpfr mpfunf.mod
