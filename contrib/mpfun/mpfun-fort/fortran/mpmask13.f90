!*****************************************************************************

function mpmask13 (b)

!  Revision date:  11 Feb 2021

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2021 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF THIS ROUTINE:
!    This convoluted-looking code tests whether B has more than 40
!    significant bits, returning the absolute value of B, with lower 13 
!    bits zeroed out. This function must be compiled separately, with lower
!    optimization, since -fast with gfortran, for instance, negates the test.

use mpfuna
implicit none
real (8) b, mpmask13, t1
t1 = mpb13x * abs (b)
mpmask13 = abs (abs (b) + t1) - abs (t1)
return
end
