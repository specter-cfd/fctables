! ---------------------------------------------------------------------
      MODULE fft_mod
! ---------------------------------------------------------------------
! Simple unoptimized fft borrowed from
! https://rosettacode.org/wiki/Fast_Fourier_transform#Fortran.
!
      use, intrinsic :: iso_fortran_env

      IMPLICIT NONE
      INTEGER,       PARAMETER :: sp=REAL32
      INTEGER,       PARAMETER :: dp=2*REAL32

      INTEGER,       PARAMETER :: wp=sp
      REAL(KIND=wp), PARAMETER :: pi=3.141592653589793238460_wp
      CONTAINS
       
! ---------------------------------------------------------------------
      RECURSIVE SUBROUTINE fft(x)
! ---------------------------------------------------------------------
! In place Cooley-Tukey FFT
!
      COMPLEX(KIND=wp), DIMENSION(:), INTENT(INOUT)  :: x
      COMPLEX(KIND=wp)                               :: t
      INTEGER                                        :: N
      INTEGER                                        :: i
      COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE    :: even, odd
     
      N = SIZE(x)
     
      IF ( N .le. 1 ) RETURN
     
      ALLOCATE(odd((N+1)/2))
      ALLOCATE(even(N/2))
     
      ! Divide
      odd =x(1:N:2)
      even=x(2:N:2)
     
      ! Conquer
      CALL fft(odd)
      CALL fft(even)
     
      ! Combine
      DO i=1,N/2
         t = EXP(CMPLX(0.0_wp,-2.0_wp*pi*REAL(i-1,wp)/REAL(N,wp),KIND=wp))*even(i)
         x(i)     = odd(i) + t
         x(i+N/2) = odd(i) - t
      END DO
     
      DEALLOCATE(odd)
      DEALLOCATE(even)
     
      END SUBROUTINE fft
     
      END MODULE fft_mod



! =====================================================================
      PROGRAM FC_TEST
! =====================================================================
      USE fft_mod
    
      REAL (KIND=REAL64)  , DIMENSION(:,:), ALLOCATABLE  :: Q,A

      COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE    :: f
      REAL(KIND=wp)   , DIMENSION(:), ALLOCATABLE    :: fp,x,k

      INTEGER :: i, n, c, d

      CHARACTER (LEN=100) :: odir  ! Dummy variable
      CHARACTER (LEN=5)   :: cs,ds
   
      ! Total number of points, continuation points, and matching points
      NAMELIST / required / odir,d,c
      OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=required)
      CLOSE(1)

      CALL get_command_argument(1, cs)
      IF ( LEN_TRIM(cs) == 0 ) THEN
         PRINT*, "Usage: fc_test N, with N > c the total grid size.",&
            "N must be a power of 2."
         STOP
      ELSE
         READ(cs,*) n
         IF ( n <= c) THEN
            PRINT*, "Usage: fc_test N, with N > c the total grid ",&
               "size. N must be a power of 2."
            STOP
         ENDIF
      ENDIF

      n = n - c

      ALLOCATE ( A(c,d), Q(d,d) )
      ALLOCATE ( f(n+c), fp(n), x(n), k(n+c) )


      ! Load A and Q used for continuations
      WRITE(cs,'(I5)') c
      WRITE(ds,'(I5)') d
      OPEN(10, FILE=trim(odir) // '/A' // trim(adjustl(cs)) // '-' // &
              trim(adjustl(ds)) //  '.dat', &
              FORM='unformatted', ACCESS='stream')
         READ(10) A
      CLOSE(10)
      OPEN(10, FILE=trim(odir) // '/Q' // trim(adjustl(ds)) //  '.dat',&
              FORM='unformatted', ACCESS='stream')
         READ(10) Q
      CLOSE(10)
      Q = TRANSPOSE(Q)  ! Only need the transpose


      ! Grid from x=0 to x=10, function and derivative
      x      = (/( i*10.0_wp/(n-1), i=0,n-1 )/)
      f(1:n) = BESSEL_J0(20*x)*EXP(3*x/x(n))
      fp     = -20*BESSEL_J1(20*x)*EXP(3*x/x(n)) + &
                  BESSEL_J0(20*x)*EXP(3*x/x(n))*3/x(n)

      ! Wavenumbers
      DO i = 1,(n+c)/2
         k(i) = REAL(i-1,KIND=wp)
         k(i+(n+c)/2) = REAL(i-(n+c)/2-1,KIND=wp)
      END DO

      k = 2*pi*k/( 10.0_wp/(n-1)*(n+c) )  ! k = 2*pi*k_ind/L


      ! Perform continuation. First contribution of last d points,
      ! then add contribution of first d points.
      f(n+1:n+c) = MATMUL(A, MATMUL(Q, f(n-d+1:n)) )
      f(n+1:n+c) = f(n+1:n+c) + MATMUL(A(c:1:-1,:), MATMUL(Q, f(d:1:-1)) )

      ! Transform, derivate, antitransform
      ! ifft can be obtained as conj(fft(conj(f_k))).
      ! For real signals it is redundant anyway.
      CALL fft(f)
      f = CMPLX(0.0_wp, 1.0_wp)*k*f
      f = CONJG(f)
      CALL fft(f)
      f = CONJG(f)/(n+c)


      ! Print error comparison
      PRINT*, "Maximum relative error: ", &
         MAXVAL(ABS(REAL(f(1:n),kind=wp) -  fp)/ABS(fp))
      PRINT*, "Average relative error: ", &
         SUM(ABS(REAL(f(1:n),kind=wp) -  fp)/ABS(fp))/n

      END PROGRAM FC_TEST
