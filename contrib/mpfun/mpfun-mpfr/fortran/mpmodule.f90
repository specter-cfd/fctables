
!  MPFUN-MPFR: An MPFR-based arbitrary precision computation package
!  Main module (module MPMODULE) -- references all other modules for user.

!  Revision date:  12 Feb 2021

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

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
    
!  DESCRIPTION OF THIS MODULE (MPMODULE):
!    This module links all lower-level modules and is the connection between
!    user codes and the lower modules.  See documentation for details.

module mpmodule
use mpfuna
use mpfunf
use mpfung
use mpfunh

end module mpmodule

