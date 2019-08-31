!=================================================================
      PROGRAM FCTABLES
!=================================================================
! Program to compute FC-Gram tables for usage in NiFTy code.
! The related algorithm relies on arbitrary precision arithmetic
! provided here by the MPFUN library, by David H. Bailey and 
! available at: https://www.davidhbailey.com/dhbsoftware/.
!
! Details about FC-Gram can be found at
! https://doi.org/10.1016/j.jcp.2009.11.020
!
! IMPORTANT NOTE: THIS PROGRAM CAN ONLY RUN ON A SINGLE NODE.
!
! Required parameters:
!   - d  = matching points.
!   - c  = points in the continuating region.
!
! Optional parameters:
!   - z  = points in the "zero" region.
!   - e  = points in the "extended" region.
!   - o  = multiplying factor for the oversampled grid.
!   - bw = bandwith reduction in the trigonometric polynomial.
!
!      
! 2019 Mauro Fontana.
!      Department of Physics, 
!      Exact and Natural Sciences Faculty.
!      University of Buenos Aires.
!      e-mail: mfontana@df.uba.ar 
!=================================================================
   
      USE mpmodule
      use mprlinalg

      IMPLICIT NONE
   
      ! Vandermonde matrices, QR decomposition matrices,
      ! Trigonometric polynomial matrices and SVD matrices.
      TYPE (mp_real), DIMENSION(:,:), ALLOCATABLE :: P,Q,R
      TYPE (mp_real), DIMENSION(:,:), ALLOCATABLE :: oP,oQ
      TYPE (mp_real), DIMENSION(:,:), ALLOCATABLE :: T1,T2,A
      TYPE (mp_real), DIMENSION(:,:), ALLOCATABLE :: U,V

      ! Double precision versions of Q and A (OUTPUT)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Qd,Ad

      ! Spectra, coarse grids, oversampled grids, wavenumbers
      ! and temporal variables for solving the linear system.
      TYPE (mp_real), DIMENSION(:), ALLOCATABLE :: s
      TYPE (mp_real), DIMENSION(:), ALLOCATABLE :: grid_n,grid_d,grid_z
      TYPE (mp_real), DIMENSION(:), ALLOCATABLE :: ogrid_d,ogrid_z
      TYPE (mp_real), DIMENSION(:), ALLOCATABLE :: wn
      TYPE (mp_real), DIMENSION(:), ALLOCATABLE :: x,b

      ! Parameters and auxiliary variables
      TYPE (mp_real)   :: mprmp,mprmq

      INTEGER :: d,c,z,e,n
      INTEGER :: o,od,oc,oz,on
      INTEGER :: bw,nwn
      INTEGER :: i,j,k,l


      CHARACTER (LEN=100) :: odir
      CHARACTER (LEN=25)  :: ph
      CHARACTER (LEN=5)   :: cs,ds

      ! Namelists for parameters
      NAMELIST / required / odir,d,c
      NAMELIST / option / z,o,bw,e

      DIGS = _DIGITS
      eigenmin = mpreal('1.e-25')

      OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=required)
         READ(1,NML=option)
      CLOSE(1)

      ! Fill optional parameters if not defined.
      IF ( z  .eq. 0 )  z  = int(0.6*c)
      IF ( o  .eq. 0 )  o  = 10
      IF ( bw .eq. 0 )  bw = int(C/15.0 + 3)
      IF ( e  .eq. 0 )  e  = 25

      ! Assign derived parameters
      n = d+c+z+e

      od = o*(d-1) + 1
      oc = o*(c-1) + 1
      oz = o*(z-1) + 1
      on = o*n

      nwn = 2*(n/2-bw)+1


      ! Allocate first batch of matrices and all the vectors
      ALLOCATE( P(d,d), Q(d,d), R(d,d) )
      ALLOCATE( oP(od,d), oQ(od,d) )
      ALLOCATE( T1(od+oz,2*nwn), U(od+oz,2*nwn), V(2*nwn,2*nwn) )
      ALLOCATE( T2(c,2*nwn), A(c,d) )

      ALLOCATE( Qd(d,d), Ad(C,d) )

      ALLOCATE( s(2*nwn) )
      ALLOCATE( grid_n(n), grid_d(d), grid_z(z) )
      ALLOCATE( ogrid_d(od), ogrid_z(oz) )
      ALLOCATE( wn(nwn) )
      ALLOCATE( x(2*nwn), b(od+oz) )

      ! Populate grids
      WRITE(ph,'(i25)') n
      mprmp = 1/mpreal(ph//'.')

      WRITE(ph,'(i25)') on
      mprmq = 1/mpreal(ph//'.')
      
      grid_n  = (/( (i-1)*mprmp, i=1,n )/)
      grid_d  = (/( (i-1)*mprmp, i=1,d )/)
      grid_z  = (/( (d+C+i-1)*mprmp, i=1,z )/)

      ogrid_d = (/( (i-1)*mprmq, i=1,od )/)
      ogrid_z = (/( (d+C)*mprmp + (i-1)*mprmq, i=1,oz )/)

      wn      = (/( i*mpreal('1.'), i=-n/2+bw,n/2-bw )/)


      ! Construct Vandermonde matrices
      DO j=1,d
         DO i=1,d
            P(i,j) = grid_d(i)**(j-1)
         ENDDO
         DO i=1,od
            oP(i,j) = ogrid_d(i)**(j-1)
         ENDDO
      ENDDO
      ! Manual fix for mpreal's 0.0^0 == 0.0.
      P(1,1) = mpreal('1.'); oP(1,1) = mpreal('1.')


      ! Get QR decomposition of Vandermonde matrix in coarse grid.
      CALL mprqr(P,Q,R)

      ! Use R to "interpolate" Q in the finer grid.
      oQ = MPRSOLVEU(R,oP)



      ! Construct generic trigonometric polynomial in [0, 1) domain
      ! over oversampled matching and zero grids.
      DO j = 1,nwn
         DO i=1,od
            T1(i,j) = cos(2*mppi()*ogrid_d(i)*wn(j))
            T1(i,j+nwn) = sin(2*mppi()*ogrid_d(i)*wn(j))
         ENDDO
         DO i=1,oz
            T1(od+i,j) = cos(2*mppi()*ogrid_z(i)*wn(j))
            T1(od+i,j+nwn) = sin(2*mppi()*ogrid_z(i)*wn(j))
         ENDDO
      ENDDO

      ! Get SVD decomposition of T1
      CALL mprsvd(T1,U,s,V)


      ! Construct generic trigonometric polynomial in [0,1) domain
      ! over the continuation points
      DO j = 1,nwn
         DO i=1,C
            T2(i,j) = cos(2*mppi()*grid_n(d+i)*wn(j))
            T2(i,j+nwn) = sin(2*mppi()*grid_n(d+i)*wn(j))
         ENDDO
      ENDDO


      ! Solve system of equations for each orthogonal polynomial
      DO i=1,d
         ! RHS
         b(1:od) = (/( oQ(k,i), k=1,od )/)
         b(od+1:od+oz) = (/( mpreal('0.'), k=1,oz )/)

         ! Compute D*V^t*x
         x = MPRMATVEC(MPRTRANSPOSE(U), b)

         ! Compute V^t*x. (i.e. divide x_i by s_i)
         DO j = 1,2*nwn
            IF ( s(j) >= eigenmin ) x(j) = x(j)/s(j)
         ENDDO

         ! Compute x
         x = MPRMATVEC(V, x)

         ! Get continuation coefficients
         A(:,i) = MPRMATVEC(T2, x)
      ENDDO


      ! Convert to DOUBLE and save
      DO j=1,d
         DO i=1,C
            Ad(i,j) = DBLE(A(i,j))
         ENDDO
      ENDDO

      ! Convert to DOUBLE and save
      DO j=1,d
         DO i=1,d
            Qd(i,j) = DBLE(Q(i,j))
         ENDDO
      ENDDO

      WRITE(cs,'(I5)') c
      WRITE(ds,'(I5)') d

      OPEN(10, FILE=trim(odir) // '/A' // trim(adjustl(cs)) // '-' // &
              trim(adjustl(ds)) //  '.dat', &
              FORM='unformatted', ACCESS='stream')
         WRITE(10) Ad
      CLOSE(10)

      OPEN(10, FILE=trim(odir) // '/Q' // trim(adjustl(ds)) //  '.dat', &
              FORM='unformatted', ACCESS='stream')
         WRITE(10) Qd
      CLOSE(10)


      ! Finish
      DEALLOCATE( P, oP, R)
      DEALLOCATE( Q, oQ )
      DEALLOCATE( T1, U, V )
      DEALLOCATE( T2, A )

      DEALLOCATE( Qd, Ad )

      DEALLOCATE( s )
      DEALLOCATE( grid_n, grid_d, grid_z)
      DEALLOCATE( ogrid_d, ogrid_z )
      DEALLOCATE( wn )
      DEALLOCATE( x, b )



      END PROGRAM FCTABLES
