! ---------------------------------------------------------------------
      MODULE fft_mod
! ---------------------------------------------------------------------
! Simple unoptimized fft borrowed from
! https://rosettacode.org/wiki/Fast_Fourier_transform#Fortran.
!
      IMPLICIT NONE
      INTEGER,       PARAMETER :: wp=selected_real_kind(15,300)
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
    
      REAL (KIND=wp)  , DIMENSION(:,:), ALLOCATABLE  :: Q,A,Qder

      COMPLEX(KIND=wp), DIMENSION(:), ALLOCATABLE    :: f
      REAL(KIND=wp)   , DIMENSION(:), ALLOCATABLE    :: fp,x,k

      REAL(KIND=wp)  :: dxp

      INTEGER :: i, n, c, d, tkind

      CHARACTER (LEN=100) :: odir  ! Dummy variable
      CHARACTER (LEN=5)   :: Cstr,dstr
   
      ! Total number of points, continuation points, and matching points
      ! Neum is currently unused
      NAMELIST / required / odir,d,c,tkind
      OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=required)
      CLOSE(1)

      CALL get_command_argument(1, Cstr)
      IF ( LEN_TRIM(Cstr) == 0 ) THEN
         PRINT*, "Usage: fc_test N, with N > c the total grid size.",&
            "N must be a power of 2."
         STOP
      ELSE
         READ(Cstr,*) n
         IF ( n <= c) THEN
            PRINT*, "Usage: fc_test N, with N > c the total grid ",&
               "size. N must be a power of 2."
            STOP
         ENDIF
      ENDIF

      n = n - c


      ALLOCATE ( f(n+c), fp(n), x(n), k(n+c) )
      ALLOCATE ( Q(d,d) )

      ! Grid from x=0 to x=10, function and derivative
      x      = (/( i*1.0_wp/(n-1), i=0,n-1 )/)
      f(1:n) = BESSEL_J0(20*x)*EXP(3*x/x(n))
      fp     = -20*BESSEL_J1(20*x)*EXP(3*x/x(n)) + &
                  BESSEL_J0(20*x)*EXP(3*x/x(n))*3/x(n)
  
      ! Load Q 
      WRITE(dstr,'(I5)') d
      OPEN(10, FILE=trim(odir) // '/Q' // trim(adjustl(dstr)) //  '.dat',&
              FORM='unformatted', ACCESS='stream')
         READ(10) Q
      CLOSE(10)

      IF (tkind .eq. 0) THEN
100      FORMAT( "Testing derivative estimation using ", i0, " matching ", &
                 "points and ", i0, " continuation points.")
         WRITE(*,100) d, c

         ! Load A
         ALLOCATE( A(c,d) )
         WRITE(Cstr,'(I5)') c
         OPEN(10, FILE=trim(odir) // '/A' // trim(adjustl(Cstr)) // '-' // &
                trim(adjustl(dstr)) //  '.dat', &
                FORM='unformatted', ACCESS='stream')
           READ(10) A
         CLOSE(10)
   
         ! Wavenumbers
         DO i = 1,(n+c)/2
            k(i) = REAL(i-1,KIND=wp)
            k(i+(n+c)/2) = REAL(i-(n+c)/2-1,KIND=wp)
         END DO
   
         k = 2*pi*k/( x(2)*(n+c) )  ! k = 2*pi*k_ind/L
   
         ! Perform continuation. First contribution of last d points,
         ! then add contribution of first d points.
         f(n+1:n+c) = MATMUL(A, MATMUL(TRANSPOSE(Q), f(n-d+1:n)) )
         f(n+1:n+c) = f(n+1:n+c) + &
                        MATMUL(A(c:1:-1,:), MATMUL(TRANSPOSE(Q), f(d:1:-1)) )
   
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

      ! Neumann reconstruction
      ELSEIF (tkind .eq. 1) THEN
101   FORMAT( "Testing endpoints reconstruction from the prescription of ", &
              "the normal derivative employing ", i0, " matching points.") 
      WRITE(*,101) d

         ! Load Q2n
         ALLOCATE(Qder(d,d))

         OPEN(10, FILE=trim(odir) // '/Qn' // trim(adjustl(dstr)) //  '.dat',&
                 FORM='unformatted', ACCESS='stream')
            READ(10) dxp
            READ(10) Qder
         CLOSE(10)

         fp(2:d) = f(2:d)
         fp(n-d+1:n-1) = f(n-d+1:n-1)
         fp(1) = fp(1)*(-x(2)/dxp)
         fp(n) = fp(n)*(x(2)/dxp)

         fp(1) = DOT_PRODUCT(Q(d,:), MATMUL(TRANSPOSE(Qder), fp(d:1:-1)))
         fp(n) = DOT_PRODUCT(Q(d,:), MATMUL(TRANSPOSE(Qder), fp(n-d+1:n)))

         ! Print error 
         PRINT*, "Relative reconstruction error at x=0: ", &
             ABS((fp(1)-f(1))/f(1))
         PRINT*, "Relative reconstruction error at x=L: ", &
             ABS((fp(n)-f(n))/f(n))

      ! Second normal derivative reconstruction
      ELSEIF (tkind .eq. 2) THEN
102   FORMAT( "Testing endpoints reconstruction from the prescription of ", &
              "the 2nd normal derivative employing ", i0, " matching points.")
      WRITE(*,102) d

         ! Load Q2n
         ALLOCATE(Qder(d,d))

         OPEN(10, FILE=trim(odir) // '/Q2n' // trim(adjustl(dstr)) //  '.dat',&
                 FORM='unformatted', ACCESS='stream')
            READ(10) dxp
            READ(10) Qder
         CLOSE(10)

         fp(2:d) = f(2:d)
         fp(n-d+1:n-1) = f(n-d+1:n-1)

         fp(1) = -20*20/2*(BESSEL_J0(20*x(1))-BESSEL_JN(2,20*x(1)))*EXP(3*x(1)/x(n)) &
                  -2*(20*BESSEL_J1(20*x(1))*EXP(3*x(1)/x(n))*3/x(n)) &
                  +BESSEL_J0(20*x(1))*EXP(3*x(1)/x(n))*(3/x(n))**2

         fp(n) = -20*20/2*(BESSEL_J0(20*x(n))-BESSEL_JN(2,20*x(n)))*EXP(3*x(n)/x(n)) &
                  -2*(20*BESSEL_J1(20*x(n))*EXP(3*x(n)/x(n))*3/x(n)) &
                  +BESSEL_J0(20*x(n))*EXP(3*x(n)/x(n))*(3/x(n))**2

         fp(1) = fp(1)*(x(2)/dxp)**2
         fp(n) = fp(n)*(x(2)/dxp)**2

         fp(1) = DOT_PRODUCT(Q(d,:), MATMUL(TRANSPOSE(Qder), fp(d:1:-1)))
         fp(n) = DOT_PRODUCT(Q(d,:), MATMUL(TRANSPOSE(Qder), fp(n-d+1:n)))

         ! Print error 
         PRINT*, "Relative reconstruction error at x=0: ", &
             ABS((fp(1)-f(1))/f(1))
         PRINT*, "Relative reconstruction error at x=L: ", &
             ABS((fp(n)-f(n))/f(n))
      ENDIF

      ! Clean up
      DEALLOCATE( x,f,fp )
      DEALLOCATE( Q )
      IF (tkind .eq. 0) THEN
         DEALLOCATE( A ) 
      ELSEIF ( tkind .ge. 0) THEN
         DEALLOCATE( Qder )
      ENDIF


      END PROGRAM FC_TEST
