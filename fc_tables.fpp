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
      TYPE(mp_real), DIMENSION(:,:), ALLOCATABLE :: P,Q,R
      TYPE(mp_real), DIMENSION(:,:), ALLOCATABLE :: Qder,Rder
#ifdef BLEND_
      TYPE(mp_real), DIMENSION(:,:), ALLOCATABLE :: oP,oQ
      TYPE(mp_real), DIMENSION(:,:), ALLOCATABLE :: T1,T2,T3,A
#endif

      ! Double precision versions of Q and A (OUTPUT)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Qd,Ad

      ! Spectra, coarse grids, oversampled grids, wavenumbers
      ! and temporal variables for solving the linear system.
      TYPE(mp_real), DIMENSION(:), ALLOCATABLE :: grid_n,grid_d
#ifdef BLEND_
      TYPE(mp_real), DIMENSION(:), ALLOCATABLE :: ogrid_d,ogrid_z
      TYPE(mp_real), DIMENSION(:), ALLOCATABLE :: grid_z,wn
      TYPE(mp_real), DIMENSION(:), ALLOCATABLE :: x,b
#endif

      ! Parameters and auxiliary variables
      TYPE (mp_real)   :: L
      TYPE (mp_real)   :: mprmp,mprmq

      INTEGER :: d,c,z,e,n
      INTEGER :: o,od,oc,oz,on
      INTEGER :: bw,nwn
      INTEGER :: i,j,k
      INTEGER :: tkind

#ifdef BLEND_
      TYPE(MPRSVDPLAN) :: plansvd
      INTEGER :: resume,sstep,iters
#endif

      INTEGER :: nth
!$    INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS

      CHARACTER (LEN=100) :: odir
      CHARACTER (LEN=25)  :: ph
      CHARACTER (LEN=5)   :: cs,ds

      ! Namelists for parameters
      NAMELIST / required / odir,d,c,tkind
      NAMELIST / option   / z,o,bw,e
#ifdef BLEND_
      NAMELIST / svd      / resume,sstep,iters
#endif

#ifdef BLEND_
      z=0; o=0; bw=0; e=0
      resume=0; iters=0; sstep=0
#endif

      OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=required)
#ifdef BLEND_
         READ(1,NML=svd)
         READ(1,NML=option)
#endif
      CLOSE(1)

#ifdef BLEND_
      ! Fill optional parameters if not defined.
      IF ( o  .le. 0 )  o  = 20
      IF ( bw .le. 0 )  bw = 4 
      IF ( z  .le. 0 )  z  = int(0.5*c)
      IF ( e  .le. 0 )  e  = c

      IF ( iters .le. 0 )  iters = 200
      IF ( sstep .le. 0 )  sstep = iters+1
#endif


      ! Construct strings with the resolutions (for saving)
      WRITE(cs,'(I5)') c
      WRITE(ds,'(I5)') d


      nth = 1
!$      nth = OMP_GET_MAX_THREADS()
200   FORMAT( "Obtaining operators for ", i0, " matching points and ", &
              i0, " continuation points") 
      WRITE(*,200) d, C
201   FORMAT( "Running in ", i0, " cores with ", i0, " digits of precision") 
      WRITE(*,201) nth, mpipl 
#ifdef BLEND_
202   FORMAT( "z = ", i0, ", o = ", i0, ", bw = ", i0, ", e = ", i0) 
      WRITE(*,202) z, o, bw, e
#endif

      ! Assign derived parameters
      n = d+c+z+e

      od = o*(d-1) + 1
      oc = o*(c-1) + 1
      oz = o*(z-1) + 1
      on = o*n

      nwn = 2*(n/2-bw)+1

      ! Allocate first batch of matrices and all the vectors
      ! NOTE: Der matrices are used only for Neumann tables, however
      ! they are small enough that allocating in common block should
      ! not incur in appreciable overhead.
      ALLOCATE( P(d,d), Q(d,d), R(d,d) )
      ALLOCATE( Qder(d,d), Rder(d,d) )
      ALLOCATE( Qd(d,d), Ad(C,d) )
      ALLOCATE( grid_n(n), grid_d(d))
#ifdef BLEND_
      ALLOCATE( oP(od,d), oQ(od,d) )
      ALLOCATE( T1(od+oz,2*nwn), T2(od+oz,2*nwn) )
      ALLOCATE( T3(c,2*nwn), A(c,d) )
      ALLOCATE( ogrid_d(od), ogrid_z(oz) )
      ALLOCATE( grid_z(z), wn(nwn) )
      ALLOCATE( x(2*nwn), b(od+oz) )
#endif

      ! Populate grids
      WRITE(ph,'(i25)') n
      ! Length of the domain
      L = mpreal('1.0')

      mprmp = L/mpreal(ph//'.')

      WRITE(ph,'(i25)') on
      mprmq = L/mpreal(ph//'.')
      
      grid_n  = (/( (i-1)*mprmp, i=1,n )/)
      grid_d  = (/( (i-1)*mprmp, i=1,d )/)

#ifdef BLEND_
      grid_z  = (/( (d+C+i-1)*mprmp, i=1,z )/)
      ogrid_d = (/( (i-1)*mprmq, i=1,od )/)
      ogrid_z = (/( (d+C)*mprmp + (i-1)*mprmq, i=1,oz )/)

      wn      = (/( i*mpreal('1.'), i=-n/2+bw,n/2-bw )/)
#endif

      ! Construct Vandermonde matrices
      DO j=1,d
         DO i=1,d
            P(i,j) = grid_d(i)**(j-1)
         ENDDO
      ENDDO
      ! Manual fix for mpreal's 0.0^0 == 0.0.
      P(1,1) = mpreal('1.')

#ifdef BLEND_
      DO j=1,d
         DO i=1,od
            oP(i,j) = ogrid_d(i)**(j-1)
         ENDDO
      ENDDO
      oP(1,1) = mpreal('1.')
#endif

      ! Get QR decomposition of Vandermonde matrix in coarse grid.
      CALL mprqr(P,Q,R)

      ! ---------------------
      ! Dirichlet Projector
      ! ---------------------
      IF ( tkind .eq. 0 ) THEN
         ! Convert to DOUBLE and save
         DO j=1,d
            DO i=1,d
               Qd(i,j) = DBLE(Q(i,j))
            ENDDO
         ENDDO

         OPEN(10, FILE=trim(odir) // '/Q' // trim(adjustl(ds)) //  '.dat', &
                 FORM='unformatted', ACCESS='stream')
            WRITE(10) Qd
         CLOSE(10)

      !---------------------
      ! Neumann Projector
      !---------------------
      ELSEIF ( tkind .eq. 1 ) THEN
         ! Modified Vandermonde matrix
         DO j=1,d
               P(d,j) = (j-1)*grid_d(d)**(j-2)
         ENDDO

         CALL mprqr(P,Qder,Rder)

         Q = MPRTRANSPOSE(Qder)
         ! Find R^-1 by using MPRSOLVU
         DO j=1,d
            DO i=1,d
               IF ( i .eq. j ) THEN
                  P(i,j) = mpreal('1.0')
               ELSE
                  P(i,j) = mpreal('0.0')
               ENDIF
            ENDDO
         ENDDO
         Qder = MPRSOLVEU(Rder,P)

         P = MPRMATMUL(Qder,Q)
         Rder = MPRMATMUL(R,P)
         Qder = MPRTRANSPOSE(Rder)

         ! Convert to DOUBLE and save both the grid spacing and Q to file
         DO j=1,d
            DO i=1,d
               Qd(i,j) = DBLE(Qder(i,j))
            ENDDO
         ENDDO
   
         OPEN(10, FILE=trim(odir) // '/Qn' // trim(adjustl(ds)) //  '.dat', &
                 FORM='unformatted', ACCESS='stream')
            WRITE(10) DBLE(grid_n(2))
            WRITE(10) Qd
         CLOSE(10)

      !---------------------
      ! Neumann^2 Projector
      !---------------------
      ELSE IF ( tkind .eq. 2 ) THEN
         ! Modified Vandermonde matrix
         DO j=1,d
               P(d,j) = (j-2)*(j-1)*grid_d(d)**(j-3)
         ENDDO

         CALL mprqr(P,Qder,Rder)

         Q = MPRTRANSPOSE(Qder)
         ! Find R^-1 by using MPRSOLVU
         DO j=1,d
            DO i=1,d
               IF ( i .eq. j ) THEN
                  P(i,j) = mpreal('1.0')
               ELSE
                  P(i,j) = mpreal('0.0')
               ENDIF
            ENDDO
         ENDDO
         Qder = MPRSOLVEU(Rder,P)

         P = MPRMATMUL(Qder,Q)
         Rder = MPRMATMUL(R,P)
         Qder = MPRTRANSPOSE(Rder)

         ! Convert to DOUBLE and save both the grid spacing and Q to file
         DO j=1,d
            DO i=1,d
               Qd(i,j) = DBLE(Qder(i,j))
            ENDDO
         ENDDO
   
         OPEN(10, FILE=trim(odir) // '/Q2n' // trim(adjustl(ds)) //  '.dat', &
                 FORM='unformatted', ACCESS='stream')
            WRITE(10) DBLE(grid_n(2))
            WRITE(10) Qd
         CLOSE(10)

      !---------------------
      ! Robin Projector
      !---------------------
      ELSE IF ( tkind .eq. 3 ) THEN  !Robin
         ! Modified Vandermonde matrix
         DO j=1,d
               P(d,j) = (j-1)*grid_d(d)**(j-2) + mpreal('1.0')*grid_d(d)**(j-1)
         ENDDO

         CALL mprqr(P,Qder,Rder)

         Q = MPRTRANSPOSE(Qder)
         ! Find R^-1 by using MPRSOLVU
         DO j=1,d
           DO i=1,d
               IF ( i .eq. j ) THEN
                  P(i,j) = mpreal('1.0')
               ELSE
                  P(i,j) = mpreal('0.0')
               ENDIF
            ENDDO
         ENDDO
         Qder = MPRSOLVEU(Rder,P)

         P = MPRMATMUL(Qder,Q)
         Rder = MPRMATMUL(R,P)
         Qder = MPRTRANSPOSE(Rder)

         ! Convert to DOUBLE and save both the grid spacing and Q to file
         DO j=1,d
            DO i=1,d
               Qd(i,j) = DBLE(Qder(i,j))
            ENDDO
         ENDDO
   
         OPEN(10, FILE=trim(odir) // '/Qr-' // trim(adjustl(ds)) //  'a1.dat', &
                 FORM='unformatted', ACCESS='stream')
            WRITE(10) DBLE(grid_n(2))
            WRITE(10) Qd
         CLOSE(10)
      ENDIF

#ifdef BLEND_
      !-------------------------------
      ! Blend-to-zero Operator
      !-------------------------------
      IF ( tkind .eq. 0 ) THEN
         ! Use R to "interpolate" Q in the finer grid.
         oQ = MPRSOLVEU(R,oP)

         ! Construct generic trigonometric polynomial in [0, 1) domain
         ! over oversampled matching and zero grids.
         DO j = 1,nwn
            DO i=1,od
               T1(i,j)     = cos(2*mppi()*ogrid_d(i)*wn(j))
               T1(i,j+nwn) = sin(2*mppi()*ogrid_d(i)*wn(j))
            ENDDO
            DO i=1,oz
               T1(od+i,j)     = cos(2*mppi()*ogrid_z(i)*wn(j))
               T1(od+i,j+nwn) = sin(2*mppi()*ogrid_z(i)*wn(j))
            ENDDO
         ENDDO

         IF ( resume .eq. 0) THEN
            CALL plansvd%init_plan(T1)
         ELSEIF ( resume .eq. 1 ) THEN
            CALL plansvd%load_plan('svd_temp.dat')
         ENDIF
   
         ! Get orthogonal decomposition of T1
         IF ( plansvd%steps .lt. iters ) THEN
400         FORMAT( "Performing U*V^t decomposition (size: ", i0, "x", i0, ")") 
            WRITE(*,400) UBOUND(T1,1), UBOUND(T1,2) 
         ENDIF
         DO WHILE (plansvd%steps .lt. iters)
            CALL plansvd%iterate(min(sstep,iters))

            ! Check orthogonality 
            mprmp = check_orthogonal(plansvd%U)
            PRINT*,
            PRINT*, ACHAR(13), 'Max. inner product between columns of U: ',&
                    DBLE(mprmp)
            mprmp = check_orthogonal(plansvd%V)
            PRINT*, ACHAR(13), 'Max. inner product between columns of V: ',&
                    DBLE(mprmp)
   
            ! Check decomposition is accurate
            T2 = MPRMATMUL(plansvd%U, MPRTRANSPOSE(plansvd%V))
            mprmp = mpreal('0.')
            DO j=1,UBOUND(T2,2)
            DO i=1,UBOUND(T2,1)
               mprmp = max(mprmp,abs(T1(i,j) - T2(i,j)))
            ENDDO
            ENDDO
            PRINT*, ACHAR(13), 'Error in U*V^t decomposition: ', DBLE(mprmp)
            PRINT*,

            IF ( sstep .le. iters ) CALL plansvd%save_plan('svd_temp.dat')
         ENDDO

         ! Get SVD of T1 and check it
         CALL plansvd%to_svd()

         DO j=1,UBOUND(plansvd%V,2)
         DO i=1,UBOUND(T2,1)
            mprmp = mpreal('0.')
            DO k=1,UBOUND(plansvd%V,1)
               mprmp = mprmp + plansvd%U(i,k)*plansvd%sigma(k) * plansvd%V(j,k)
            ENDDO
            T2(i,j) = mprmp
         ENDDO
         ENDDO

         mprmp = mpreal('0.')
         DO j=1,UBOUND(T2,2)
         DO i=1,UBOUND(T2,1)
            mprmp = max(mprmp,abs(T1(i,j) - T2(i,j)))
         ENDDO
         ENDDO
         PRINT*, 'Error in SVD decomposition: ', DBLE(mprmp)
         PRINT*,


         ! Finally, solve the required system of equations! 
         ! Construct generic trigonometric polynomial in [0,1) domain
         ! over the continuation points
         DO j = 1,nwn
            DO i=1,C
               T3(i,j)     = cos(2*mppi()*grid_n(d+i)*wn(j))
               T3(i,j+nwn) = sin(2*mppi()*grid_n(d+i)*wn(j))
            ENDDO
         ENDDO
   
   
         ! Solve system of equations for each orthogonal polynomial
         DO i=1,d
            ! RHS
            b(1:od) = (/( oQ(k,i), k=1,od )/)
            b(od+1:od+oz) = (/( mpreal('0.'), k=1,oz )/)
   
            ! Compute D*V^t*x
            x = MPRMATVEC(MPRTRANSPOSE(plansvd%U), b)
   
            ! Compute V^t*x. (i.e. divide x_i by sigma_i)
            DO j = 1,2*nwn
               IF ( plansvd%sigma(j) >= plansvd%eigenmin ) THEN
                   x(j) = x(j)/plansvd%sigma(j)
               ENDIF
            ENDDO
   
            ! Compute x
            x = MPRMATVEC(plansvd%V, x)
   
            ! Get continuation coefficients
            A(:,i) = MPRMATVEC(T3, x)
         ENDDO
   
   
         ! Convert to DOUBLE and save
         DO j=1,d
            DO i=1,C
               Ad(i,j) = DBLE(A(i,j))
            ENDDO
         ENDDO
   
   
         OPEN(10, FILE=trim(odir) // '/A' // trim(adjustl(cs)) // '-' // &
                 trim(adjustl(ds)) //  '.dat', &
                 FORM='unformatted', ACCESS='stream')
            WRITE(10) Ad
         CLOSE(10)
   
      ENDIF
#endif

      ! Finish
      DEALLOCATE( P, Q, R)
      DEALLOCATE( Qder, Rder)
      DEALLOCATE( Qd, Ad )
      DEALLOCATE( grid_n, grid_d)
#ifdef BLEND_
      DEALLOCATE( oP, oQ )
      DEALLOCATE( T1, T2, T3, A )
      DEALLOCATE( ogrid_d, ogrid_z )
      DEALLOCATE( grid_z, wn )
      DEALLOCATE( x, b )
      CALL plansvd%destroy_plan()
#endif

      END PROGRAM FCTABLES




!            PRINT*, "m", plansvd%m, plansvd2%m
!            PRINT*, "n", plansvd%n, plansvd2%n
!            PRINT*, "digs", plansvd%digs, plansvd2%digs
!            PRINT*, "step", plansvd%step, plansvd2%step
!
!            PRINT*, "eigenmin", dble(plansvd%eigenmin), dble(plansvd2%eigenmin)
!            PRINT*, "eps", dble(plansvd%eps), dble(plansvd2%eps)
!            PRINT*, "gmax", dble(plansvd%gmax), dble(plansvd2%gmax)
!
!            PRINT*, dble(plansvd%sweep(1,od/2,od/2)), dble(plansvd2%sweep(1,od/2,od/2))
!            PRINT*, dble(plansvd%sweep(2,od/2,od/2)), dble(plansvd2%sweep(2,od/2,od/2))
!            PRINT*, "sweep", MAXVAL(ABS(plansvd%sweep - plansvd%sweep))
!
!            PRINT*, dble(plansvd%U(od/2,nwn)), dble(plansvd2%U(od/2,nwn))
!            mprmp = mpreal('0.')
!            DO j=1,2*nwn
!            DO i=1,od
!               mprmp = max(mprmp, abs(plansvd%U(i,j)-plansvd2%U(i,j)))
!            ENDDO
!            ENDDO
!            PRINT*, "U", DBLE(mprmp)
!
!            PRINT*, dble(plansvd%sigma(nwn)), dble(plansvd2%sigma(nwn))
!            mprmp = mpreal('0.')
!            DO i=1,2*nwn
!               mprmp = max(mprmp, abs(plansvd%sigma(i)-plansvd2%sigma(i)))
!            ENDDO
!            PRINT*, "sigma", DBLE(mprmp)
!
!            PRINT*, dble(plansvd%V(nwn,nwn)), dble(plansvd2%V(nwn,nwn))
!            mprmp = mpreal('0.')
!            DO j=1,2*nwn
!            DO i=1,2*nwn
!               mprmp = max(mprmp, abs(plansvd%V(i,j)-plansvd2%V(i,j)))
!            ENDDO
!            ENDDO
!            PRINT*, "V", DBLE(mprmp)

