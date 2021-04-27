!*****************************************************************************

!  MPFUN-MPFR: An MPFR-based Fortran arbitrary precision computation package
!  Global data definition module (module MPFUNA)

!  Revision date:  18 Feb 2021

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired) and University of California, Davis
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2021 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with several
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!    This version differs from the MPFUN20-Fort version by the same author in that
!    this is based on MPFR, which presently is the fastest available low-level
!    package for high-precision floating-point computation.  At the Fortran user
!    level, most application codes written for MPFUN20-Fort may be compiled with
!    MPFUN-MPFR -- i.e., MPFUN-MPFR is plug-compatible with MPFUN20-Fort.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package," 
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2020.pdf.

!  DESCRIPTION OF THIS MODULE (MPFUNA):
!    This module contains some global parameter definitions.

module mpfuna
implicit none

!----------------------------------------------------------------------------

!   Integer(4) constants:

!   Name     Default   Description
!   mpndpw       18    Largest n such that 10^n <= mpbdx (see below).
!   mpldb         6    Logical device number for output of error messages.
!   mpnbt        64    Number of significant bits in one mantissa word.
!   mpnstr     2048    Maximum length of certain input character strings.
!                        See usage in mpinp of module MPFUNC.
!   mpoutln      80    Length of output lines.  See usage in mpout of MPFUNC.
!   mprknd        8    Kind parameter for double precision.
!   mprnd         0    MPFR rounding mode.

!   Integer(8) constants

!   mpnan     -2^63+2  MPFR NAN flag.
!   mpzero    -2^63+1  MPFR zero flag.

!   Real (8) constants:

!   mpbdx      2^64    2^mpnbt, the radix for MP numbers.
!   mpbx2      2^128   Square of radix.
!   mpb13x     2^13    Constant for checking the 40-bit restriction.
!                        See usage in mpdmc40, mpmuld40 and mpdivd40 of MPFUNB.
!   mpdpw  Log10(2^64) DP approximation to number of digits per mantissa word.
!   mpexpmx   2^31     Largest permissible exponent.
!   mprdx   2^(-64)    Reciprocal of radix.
!   mprx2   2^(-128)   Reciprocal of square of radix.
!   mprandx  5501758857861179  Multiplier for pseudorandom number generator.

integer, public:: mpndpw, mpldb, mpnbt, mpnstr, mpoutln, mprknd, mprnd
integer(8), public:: mpnan, mpzero
parameter (mpndpw = 18, mpldb = 6, mpnbt = 64, mpnstr = 2048, mpoutln = 80, &
  mprknd = 8, mprnd = 0)
real (mprknd), public:: mpbdx, mpbx2, mpdpw, mpexpmx, mprdx, mprx2, mpb13x, &
  mprandx
parameter (mpnan = -9223372036854775806_8, mpzero = -9223372036854775807_8)
parameter (mpbdx = 2.d0 ** mpnbt, mpbx2 = mpbdx**2, mpdpw = 19.265919722494796d0, &
  mpexpmx = 2.d0**31, mprdx = 0.5d0 ** mpnbt, mprx2 = mprdx**2, &
  mpb13x = 2.d0 ** 13, mprandx = 5501758857861179.d0)

!  These are for debug and timing (for internal use only):

integer mpdeb(10)
real (mprknd) mptm(10), mpffterr
data mpffterr / 0.d0/

end module mpfuna


