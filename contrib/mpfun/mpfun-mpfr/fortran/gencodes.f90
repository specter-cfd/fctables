!  program gencodes

!  Revision date:  11 Mar 2021

!  AUTHOR:
!     David H. Bailey
!     Lawrence Berkeley National Lab (retired) and University of California, Davis
!     Email: dhbailey@lbl.gov

!  This reads an input source template from file fname (set in parameter statement
!  below) and generates the appropriate set of variant source files. Thus when one
!  of these files is updated, each of the variant files can be updated as well.

!  If fname = "mpfunb.f90" or "mpfunbq.f90" (for MPFUN20-Fort) then this program generates
!    mpfunb-new.f90 and mpfunbq-new.f90.
!  
!  If fname = "mpfung1.f90" or one of the other three "g" files, then it generates
!    mpfung1-new.f90, mpfung2-new.f90, mpfungq1-new.f90 and mpfungq2-new.f90.

!  If fname = "mpfunh1.f90" or one of the other three "h" files, then it generates
!    mpfunh1-new.f90, mpfunh2-new.f90, mpfunhq1-new.f90 and mpfunhq2-new.f90.

program gencodes
implicit none
character(32) fname
parameter (fname ='mpfung1.f90')

if (fname(6:6) == 'b') then
  call genb (fname)
  call genbq (fname)
elseif (fname(6:6) == 'g' .or. fname(6:6) == 'h') then
  call gen1 (fname)
  call gen2 (fname)
  call genq1 (fname)
  call genq2 (fname)
else
  write (6, *) 'filename not recognized'
endif

stop
end

subroutine genb (fname)
implicit none
integer i, l1, l2
character(32) fname, fname2
character(256) line1, line2

!   Generate mpfunb-new.f90:

open (11, file = fname)
rewind 11
fname2 = 'mpfunb-new.f90'
open (12, file = fname2)
rewind (12)

101 continue

read (11, '(a)', end=119) line1

if (line1(4:10) == 'Variant') then
  write (12, 11)
11 format (&
    '!  Variant standard: real(16) is NOT supported.')
  goto 101
elseif (line1 == '!>  If real(16) is supported, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
else
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)
  goto 101
endif

119 continue

close (12)
return

990 continue

write (12, *) 'unexpected end of file'
stop
end

subroutine genbq (fname)
implicit none
integer i, l1, l2
character(32) fname, fname2
character(256) line1, line2

!   Generate mpfunbq-new.f90:

open (11, file = fname)
rewind 11
fname2 = 'mpfunbq-new.f90'
open (12, file = fname2)
rewind (12)

101 continue

read (11, '(a)', end=119) line1

if (line1(4:10) == 'Variant') then
  write (12, 11)
11 format (&
    '!  Variant Q: real(16) IS supported.')
  goto 101
elseif (line1 == '!>  If real(16) is supported, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
else
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)
  goto 101
endif

119 continue

close (12)
return

990 continue

write (12, *) 'unexpected end of file'
stop
end

subroutine gen1 (fname)
implicit none
integer i, l1, l2
character(32) fname, fname2
character(256) line1, line2

!   Generate mpfung1-new.f90 or mpfunh1-new.f90:

open (11, file = fname)
rewind 11
write (fname2, '(a6,"1-new.f90")') fname(1:6)
open (12, file = fname2)
rewind (12)

101 continue

read (11, '(a)', end=119) line1

if (line1(4:10) == 'Variant') then
  write (12, 11)
11 format (&
    '!  Variant 1: Precision level arguments are NOT required; real(16) is NOT supported.')
  goto 101
elseif (line1 == '!>  In variant #1, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101 
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
elseif (line1 == '!>  If real(16) is supported, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
else
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)
  goto 101
endif

119 continue

close (12)
return

990 continue

write (12, *) 'unexpected end of file'
stop
end

subroutine gen2 (fname)
implicit none
integer i, l1, l2
character(32) fname, fname2
character(256) line1, line2

!   Generate mpfung2-new.f90 or mpfunh2-new.f90:

open (11, file = fname)
rewind 11
write (fname2, '(a6,"2-new.f90")') fname(1:6)
open (12, file = fname2)
rewind (12)

101 continue

read (11, '(a)', end=119) line1
! l1 = len_trim (line1)
! write (12, '(a)') line1(1:l1)

if (line1(4:10) == 'Variant') then
  write (12, 11)
11 format (&
    '!  Variant 2: Precision level arguments ARE required; real(16) is NOT supported.')
  goto 101
elseif (line1 == '!>  In variant #1, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101 
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
elseif (line1 == '!>  If real(16) is supported, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
else
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)
  goto 101
endif

119 continue

close (12)
return

990 continue

write (12, *) 'unexpected end of file'
stop
end

subroutine genq1 (fname)
implicit none
integer i, l1, l2
character(32) fname, fname2
character(256) line1, line2

!   Generate mpfungq1-new.f90 or mpfunhq1-new.f90:

open (11, file = fname)
rewind 11
write (fname2, '(a6,"q1-new.f90")') fname(1:6)
open (12, file = fname2)
rewind (12)

101 continue

read (11, '(a)', end=119) line1
! l1 = len_trim (line1)
! write (12, '(a)') line1(1:l1)

if (line1(4:10) == 'Variant') then
  write (12, 11)
11 format (&
    '!  Variant Q1: Precision level arguments are NOT required; real(16) IS supported.')
  goto 101
elseif (line1 == '!>  In variant #1, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101 
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
elseif (line1 == '!>  If real(16) is supported, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
else
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)
  goto 101
endif

119 continue

close (12)
return

990 continue

write (12, *) 'unexpected end of file'
stop
end

subroutine genq2 (fname)
implicit none
integer i, l1, l2
character(32) fname, fname2
character(256) line1, line2

!   Generate mpfungq2-new.f90 or mpfunhq2-new.f90:

open (11, file = fname)
rewind 11
write (fname2, '(a6,"q2-new.f90")') fname(1:6)
open (12, file = fname2)
rewind (12)

101 continue

read (11, '(a)', end=119) line1
! l1 = len_trim (line1)
! write (12, '(a)') line1(1:l1)

if (line1(4:10) == 'Variant') then
  write (12, 11)
11 format (&
    '!  Variant Q2: Precision level arguments ARE required; real(16) IS supported.')
  goto 101
elseif (line1 == '!>  In variant #1, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101 
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
elseif (line1 == '!>  If real(16) is supported, uncomment these lines:') then
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:12) /= '!  Otherwise' .and. line1(1:1) == '!') then
      line2 = line1(2:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') then
      goto 101
    elseif  (line1(1:12) == '!  Otherwise') then
      exit
    endif
  enddo

  do i = 1, 1000
    read (11, '(a)', end=990) line1
    l1 = len_trim (line1)
    if (line1(1:3) /= '!>>' .and. line1(1:1) /= '!') then
      line2 = '!' // line1(1:l1)
      line1 = line2
      l1 = len_trim (line1)
    endif
    write (12, '(a)') line1(1:l1)
    if (line1(1:3) == '!>>') goto 101
  enddo
else
  l1 = len_trim (line1)
  write (12, '(a)') line1(1:l1)
  goto 101
endif

119 continue

close (12)
return

990 continue

write (12, *) 'unexpected end of file'
stop
end
