!=================================================================
      MODULE mprlinalg
!=================================================================
!     Module containing simple and non-optimized linear algebra
!     algorithms in arbitrary precision, compatible with MPFUN
!     library. Necessary to construct FC-Gram tables.
!     
! 2019 Mauro Fontana.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mfontana@df.uba.ar 
!=================================================================

      USE mpmodule

      IMPLICIT NONE

      ! Eigenvalues below eigenmin are treated as zero
      TYPE(mp_real) :: eigenmin

      ! Number of decimal precision digits
      INTEGER :: DIGS

      CONTAINS

!=================================================================
      FUNCTION MPRMATVEC(A,b) RESULT(out)
!-----------------------------------------------------------------
!     Compute matrix product of matrix A and column vector b
!     with arbitrary precision.
!
!     Input:
!           - A: matrix of mp_real type.
!           - b: vector of mp_real type.
!
!     Output:
!           - mp_real vector A x b.
!-----------------------------------------------------------------

      IMPLICIT NONE

      TYPE(mp_real), INTENT(IN), DIMENSION(:,:) :: A
      TYPE(mp_real), INTENT(IN), DIMENSION(:)   :: b
      TYPE(mp_real), DIMENSION(UBOUND(A,1))     :: out

      INTEGER :: m,n
      INTEGER :: i,j

      ! Error checking
      IF ( UBOUND(A,2) .ne. UBOUND(b,1) ) THEN
         PRINT*, "Non-conforming dimensions in MPRMATVEC."
         STOP
      ELSE
         m = UBOUND(A,1); n = UBOUND(A,2)
      ENDIF


      ! Initialize
      out = (/( mpreal('0.'), i=1,m )/)

      DO i=1,n
!$omp parallel do
         DO j=1,m
            out(j) = out(j) + A(j,i)*b(i)
         ENDDO
      ENDDO



      END FUNCTION MPRMATVEC


!=================================================================
      FUNCTION MPRMATMUL(A,B) RESULT(OUT)
!-----------------------------------------------------------------
!     Compute matrix product of matrices A and B with arbitrary
!     precision.
!
!     Input:
!           - A: matrix of np_real type.
!           - B: matrix of np_real type.
!
!     Output:
!           - mp_real matrix A x B.
!-----------------------------------------------------------------

      IMPLICIT NONE

      TYPE(mp_real), INTENT(IN), DIMENSION(:,:)         :: A, B
      TYPE(mp_real), DIMENSION(UBOUND(A,1),UBOUND(B,2)) :: OUT

      INTEGER :: i,j,k
      INTEGER :: l,m,n

      ! Error checking
      IF ( UBOUND(A,2) .ne. UBOUND(B,1) ) THEN
         PRINT*, "Non-conforming dimensions in MPRMATMUL."
         STOP
      ELSE
         m = UBOUND(A,1); n = UBOUND(B,2); l = UBOUND(A,2)
      ENDIF

      ! Initialize
      OUT = RESHAPE( (/( mpreal('0.'), i=1,m*n )/), SHAPE(OUT))

      ! Compute
!$omp parallel do
      DO j=1,n
         DO i=1,m
            DO k=1,l
               OUT(i,j) = OUT(i,j) + A(i,k)*B(k,j)
            ENDDO
         ENDDO
      ENDDO



      END FUNCTION MPRMATMUL


!=================================================================
      FUNCTION MPRTRANSPOSE(A) RESULT(OUT)
!-----------------------------------------------------------------
!     Transposes matrix A.
!
!     Input:
!           - A: matrix of np_real type.
!
!     Output:
!           - Transpose of A, of type mp_real.
!-----------------------------------------------------------------

      IMPLICIT NONE

      TYPE(mp_real), INTENT(IN), DIMENSION(:,:)         :: A
      TYPE(mp_real), DIMENSION(UBOUND(A,2),UBOUND(A,1)) :: OUT

      INTEGER :: m,n
      INTEGER :: i,j

      m = UBOUND(A,1); n = UBOUND(A,2)

      ! Initialize
      OUT = RESHAPE( (/( mpreal('0.'), i=1,n*m )/), SHAPE(OUT))

      ! Transpose
      DO j=1,n
         DO i=1,m
             OUT(j,i) = A(i,j)
         ENDDO
      ENDDO

      END FUNCTION MPRTRANSPOSE


!=================================================================
      FUNCTION MPRSOLVEU(U,A) RESULT(OUT)
!-----------------------------------------------------------------
!     Solves the upper triangular system X U = A for X by forward
!     substitution. This is problem is  equivalent to solving a 
!     vector equation for each row of X, A.
!
!     Note that, although possible, multithreading was not
!     implemented, as it results in slower execution for the 
!     sizes involved in the FC-Gram tables problem.
!
!     Input:
!           - U: matrix of np_real type.
!
!     Output:
!           - X: matrix of type mp_real.
!-----------------------------------------------------------------

      IMPLICIT NONE

      TYPE(mp_real), INTENT(IN), DIMENSION(:,:)         :: U, A
      TYPE(mp_real), DIMENSION(UBOUND(A,1),UBOUND(A,2)) :: OUT

      TYPE(mp_real) :: mprtmp

      INTEGER :: m,n
      INTEGER :: i,j,k

      ! Error checking
      IF ( UBOUND(A,2) .ne. UBOUND(U,1) .OR. &
             UBOUND(U,1) .ne. UBOUND(U,2) ) THEN
         PRINT*, "Non-conforming dimensions in MPRSOLVEU."
         STOP
      ELSE
         m = UBOUND(A,1); n = UBOUND(A,2)
      ENDIF

      ! Initialize
      OUT = RESHAPE( (/( mpreal('0.'), i=1,m*n )/), SHAPE(OUT))

      ! Solve
      DO i=1,m
         DO j=1,n
            mprtmp = mpreal('0.')
            DO k=1,j-1
               mprtmp = mprtmp + OUT(i,k)*U(k,j)
            ENDDO
            OUT(i,j) = (A(i,j) - mprtmp) * (1/U(j,j))
         ENDDO
      ENDDO


      END FUNCTION MPRSOLVEU


!=================================================================
      SUBROUTINE mprqr(A,Q,R)
!-----------------------------------------------------------------
!     Compute QR decomposition of matrix A with arbitrary
!     precision. The strategy employed is the simple
!     orthogonalization of columns by Gram-Schmidt's algorithm.
!     NOTE: A must be a square matrix.
!
!     Input:
!           - A: matrix of mp_real type.
!
!     Output:
!           - Q: The Q in QR, with type mp_real.
!           - R: The R in QR, with type mp_real.
!-----------------------------------------------------------------

      IMPLICIT NONE

      TYPE (mp_real), DIMENSION(:,:), INTENT(IN) :: A 
      TYPE (mp_real), DIMENSION(UBOUND(A,1),UBOUND(A,2)), INTENT(OUT) :: Q,R

      TYPE (mp_real), DIMENSION(UBOUND(A,1)) :: mprtmp
      INTEGER :: i,j,k,n

      IF ( UBOUND(A,1) .ne. UBOUND(A,2) ) THEN
         PRINT*, "mprqr only works with square matrices"
         STOP
      ELSE
         n = UBOUND(A,1)
      ENDIF

      ! Initialize arrays
      Q = RESHAPE( (/( mpreal('0.'), i=1,n*n )/), SHAPE(Q))
      R = RESHAPE( (/( mpreal('0.'), i=1,n*n )/), SHAPE(R))

      ! Orthogonalize first column
      Q(:,1) = A(:,1)

      R(1,1) = mpreal('0.')
      DO i=1,n
         R(1,1) = R(1,1) + Q(i,1)**2
      ENDDO
      R(1,1) = SQRT(R(1,1))

      Q(:,1) = (/( Q(i,1)* (1/R(1,1)), i=1,n )/)

      ! Gram-Schmidt for the remaining columns
      DO k=2,n
         Q(:,k) = (/( A(i,k), i=1,n )/)

         DO j=1,k-1
            R(j,k) = mpreal('0.')
            DO i=1,n
               R(j,k) = R(j,k) + Q(i,j)*Q(i,k)
            ENDDO
         ENDDO

         mprtmp = MPRMATVEC(Q, R(:,k))

         Q(:, k) = (/( Q(i,k) - mprtmp(i), i=1,n )/)


         R(k,k) = mpreal('0.')
         DO i=1,n
            R(k,k) = R(k,k) + Q(i,k)*Q(i,k)
         ENDDO
         R(k,k) = SQRT(R(k,k))

         Q(:,k) = (/( Q(i,k)*(1/R(k,k)), i=1,n )/)
      ENDDO

      END SUBROUTINE mprqr


!====================================================================
      SUBROUTINE mprsvd(A,U,l,V)
!--------------------------------------------------------------------
!     Compute thin SVD decomposition of matrix A with arbitrary
!     precision, employing a one sided Jacobi scheme.
!     Based on: https://doi.org/10.1137/0906007 and
!     https://doi.org/10.1137/0910023
!
!     Input:
!           - A: matrix of mp_real type. A has dimensions mxn.
!
!     Output:
!           - U: mxn matrix containing left eigenvectors of A,
!                of type mp_real.
!           - s: n-vector containing the singular values of A, of
!                type mp_real.
!           - V: nxn matrix containing right side eigenvectors of A,
!                of type mp_real.
!     
!      Note: The sweeping strategy requires n to be even, which 
!      always is in the generation of tables for FC-Gram.
!      For other uses, the sweeping topology must be changed.
!--------------------------------------------------------------------

      IMPLICIT NONE

      TYPE (mp_real), DIMENSION(:,:), INTENT(IN) :: A 
      TYPE (mp_real), DIMENSION(UBOUND(A,1),UBOUND(A,2)), INTENT(OUT) :: U
      TYPE (mp_real), DIMENSION(UBOUND(A,2),UBOUND(A,2)), INTENT(OUT) :: V
      TYPE (mp_real), DIMENSION(UBOUND(A,2)), INTENT(OUT) :: l

      INTEGER, DIMENSION(2,UBOUND(A,2)/2,UBOUND(A,2)-1) :: sweep


      TYPE (mp_real) :: c,s,t,alpha,beta,gama
      TYPE (mp_real) :: tol,eps,gmax,jgmax,igmax
      TYPE (mp_real) :: mins,tmp,tmp2

      INTEGER :: m,n
      INTEGER :: i,j,k,ii,jj
      INTEGER :: iter,stopp,istopp

      ! Columns with inner product below eps*gmax
      ! are treated as orthogonal for the current sweep.
      eps = mpreal('1.e-30')
      gmax = mpreal('0.')  ! 1st sweep uses all columns

      m = UBOUND(A,1)
      n = UBOUND(A,2)


      ! Plan for the SVD. Pag. 73 of the first paper.
      DO i=1,n/2
        sweep(1,i,1) = 2*i-1
        sweep(2,i,1) = 2*i
      ENDDO
      DO j=2,n-1
         sweep(1,1,j) = 1
         sweep(1,2,j) = sweep(2,1,j-1)
         DO i=2,n/2-1
            sweep(1,i+1,j) = sweep(1,i,j-1)
            sweep(2,i-1,j) = sweep(2,i,j-1)
         ENDDO
         sweep(2,n/2-1,j) = sweep(2,n/2,j-1)
         sweep(2,n/2,j) = sweep(1,n/2,j-1)
      ENDDO

      ! Initialization
      U = A
      DO j=1,n
         DO i=1,n
            IF ( i .eq. j ) THEN
               V(i,j) = mpreal('1.')
            ELSE
               V(i,j) = mpreal('0.')
            ENDIF
         ENDDO
      ENDDO



      ! Start sweeping. Maximum number of sweeps is detemined
      ! ad-hoc according to precision level. This isn't optimal
      ! but probably okay for the reasonable values of precision
      ! to be used.
!      DO iter=1,DIGS*25/35+50
      DO iter=1,200

      igmax = mpreal('0.')

100   FORMAT( a , 'SVD Sweep : ',i3,' out of a total of ',i3)
!      WRITE(*,100,advance="no") ACHAR(13), iter, DIGS*25/35+50
      WRITE(*,100,advance="no") ACHAR(13), iter, 200

      DO ii=1,n-1
         jgmax = mpreal('0.')
!$omp parallel do private(i,j,k,alpha,beta,gama,tmp,c,s,t)
         DO jj=1,n/2
            i = sweep(1,jj,ii)
            j = sweep(2,jj,ii)

            alpha = mpreal('0.')
            beta = mpreal('0.')
            gama = mpreal('0.')
   
            ! Solve rotation coefficients
            DO k=1,m
               alpha = alpha + U(k,i)**2
               beta = beta + U(k,j)**2
               gama = gama + U(k,i)*U(k,j)
            ENDDO
            IF ( ABS(gama) .le. eps*gmax ) THEN
               tmp = mpreal('0.')
               t = mpreal('0.')
            ELSE
               tmp = (beta - alpha)/(2*gama)
               t = SIGN(mpreal('1.'), tmp)/(ABS(tmp) + SQRT(1 + tmp**2))
            ENDIF
            c = 1/SQRT(1+t**2)
            s = t*c

            ! Update arrays
            DO k=1,m
               tmp = U(k,i)
               U(k,i) = c*tmp - s*U(k,j)
               U(k,j) = s*tmp + c*U(k,j)
            ENDDO
            DO k=1,n
               tmp = V(k,i)
               V(k,i) = c*tmp - s*V(k,j)
               V(k,j) = s*tmp + c*V(k,j)
            ENDDO
! Manual max reduction because OpenMP max reduction doesn't handle
! mpfun types.
!$omp critical            
            jgmax = MAX(jgmax, ABS(gama))
!$omp end critical
         ENDDO
         igmax = MAX(igmax,jgmax)
      ENDDO
      gmax = igmax
      ENDDO ! End sweep

      ! Compute spectrum
!$omp parallel do private(tmp,k)
      DO i=1,n
         tmp = mpreal('0.')
         DO k=1,m
            tmp = tmp + U(k,i)*U(k,i)
         ENDDO
         ! If actual eigenvalue, add it and normalize
         ! otherwise set to zero
         IF ( SQRT(tmp) >= eigenmin ) THEN
            l(i) = SQRT(tmp)
            U(:, i) = (/( U(k,i)/l(i), k=1,m )/)
         ELSE
            l(i) = mpreal('0.')
            U(:, i) = (/( mpreal('0.'), k=1,m )/)
         ENDIF
      ENDDO


      ! Check orthogonality of U
      tmp2 = mpreal('0.')
      DO i=1,n
      DO j=i+1,n
        tmp = mpreal('0.')
        DO k=1,m
           tmp = tmp + U(k,i)*U(k,j)
        ENDDO
        IF ( tmp > tmp2) THEN
        tmp2 = max(tmp2,tmp)
        ii = i
        jj = j
        ENDIF
      ENDDO
      ENDDO
      PRINT*,
      PRINT*, 'Max. inner product between columns of U: ', DBLE(tmp2)

      ! Check orthogonality of V
      tmp2 = mpreal('0.')
      DO i=1,n
        DO j=i+1,n
        tmp = mpreal('0.')
        DO k=1,n
           tmp = tmp + V(k,i)*V(k,j)
        ENDDO
        tmp2 = max(tmp2,tmp)
      ENDDO
      ENDDO
      PRINT*, 'Max. inner product between columns of V: ', DBLE(tmp2)


      END SUBROUTINE mprsvd


      END MODULE mprlinalg
