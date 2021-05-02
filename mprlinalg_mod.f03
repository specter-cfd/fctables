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

      TYPE, PUBLIC :: MPRSVDPLAN
          PRIVATE
          TYPE (mp_real), DIMENSION(:,:), ALLOCATABLE, PUBLIC   :: U,V
          TYPE (mp_real), DIMENSION(:)  , ALLOCATABLE, PUBLIC   :: sigma
          INTEGER, PUBLIC :: steps

          INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: sweep
          TYPE(mp_real)   :: eigenmin,eps,gmax
          INTEGER         :: m,n,digs

          CONTAINS
             PRIVATE
             PROCEDURE, PASS(this) :: svd_init_plan,svd_destroy_plan
             PROCEDURE, PASS(this) :: svd_load_plan,svd_save_plan
             PROCEDURE, PASS(this) :: svd_iterate,svd_to_svd
             PROCEDURE, PASS(this) :: svd_solve
             PROCEDURE, PASS(this) :: svd_get_eigenmin,svd_set_eigenmin
             PROCEDURE             :: svd_output_procedure,svd_input_procedure

             PROCEDURE, PUBLIC  :: init_plan    => svd_init_plan
             PROCEDURE, PUBLIC  :: destroy_plan => svd_destroy_plan
             PROCEDURE, PUBLIC  :: save_plan    => svd_save_plan
             PROCEDURE, PUBLIC  :: load_plan    => svd_load_plan
             PROCEDURE, PUBLIC  :: iterate      => svd_iterate
             PROCEDURE, PUBLIC  :: to_svd       => svd_to_svd
             PROCEDURE, PUBLIC  :: solve        => svd_solve
             PROCEDURE, PUBLIC  :: get_eigenmin => svd_get_eigenmin
             PROCEDURE, PUBLIC  :: set_eigenmin => svd_set_eigenmin
             GENERIC :: WRITE(unformatted)      => svd_output_procedure
             GENERIC :: READ(unformatted)       => svd_input_procedure
      END TYPE MPRSVDPLAN

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

      RETURN
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


      RETURN
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

      RETURN
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

      RETURN
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

      RETURN
      END SUBROUTINE mprqr

!======================================================================
      SUBROUTINE svd_init_plan(this,arr,eigenmin,eps)
!----------------------------------------------------------------------
!      Initialize a plan to compute a thin SVD decomposition of matrix
!     `arr` with arbitrary precision, employing a one sided Jacobi
!     scheme. Note that the transposed of V is computed, that is
!     A = U*sigma*V^t.
!     Based on: https://doi.org/10.1137/0906007 and
!     https://doi.org/10.1137/0910023.
!
!     Input:
!       - arr:       mxn array whose SVD wants to be computed.
!       - eigenmin:  eigenvalues under `eigenmin` are considered zero.
!       - eps     :  columns with inner product under eps*gmax are
!                    considered orthogonal for a given sweep.
!                    (details on gmax in `svd_iterate`)
!
!     Output:
!       - this: SVD plan
!
!      NOTE: The sweeping strategy requires n to be even, which always 
!     is in the generation of tables for FC-Gram. For other uses, the
!     sweeping topology must be changed.
!----------------------------------------------------------------------
      IMPLICIT NONE


      TYPE (mp_real), DIMENSION(:,:), INTENT(IN) :: arr
      CLASS(MPRSVDPLAN), INTENT(out)             :: this
      TYPE(mp_real), OPTIONAL                    :: eigenmin,eps

      INTEGER :: i,j
      INTEGER :: m,n
     
      ! Size
      m  = UBOUND(arr,1)
      n  = UBOUND(arr,2)
      this%m = m
      this%n = n

      ! Numerical parameters
      this%steps    = 0
      this%digs     = mpipl
      this%eigenmin = mpreal('1.e-50')
      this%eps      = mpreal('1.e-30')
      this%gmax     = mpreal('0.0')

      IF( PRESENT(eigenmin) ) this%eigenmin = eigenmin
      IF( PRESENT(eps) )           this%eps = eps

      ! Plan for the SVD. Pag. 73 of the first paper.
      ALLOCATE( this%sweep(2,n/2,n-1) )
      DO i=1,n/2
        this%sweep(1,i,1) = 2*i-1
        this%sweep(2,i,1) = 2*i
      ENDDO
      DO j=2,n-1
         this%sweep(1,1,j) = 1
         this%sweep(1,2,j) = this%sweep(2,1,j-1)
         DO i=2,n/2-1
            this%sweep(1,i+1,j) = this%sweep(1,i,j-1)
            this%sweep(2,i-1,j) = this%sweep(2,i,j-1)
         ENDDO
         this%sweep(2,n/2-1,j) = this%sweep(2,n/2,j-1)
         this%sweep(2,n/2,j)   = this%sweep(1,n/2,j-1)
      ENDDO

      ! Initialization of SVD
      ALLOCATE( this%U(m,n), this%sigma(n), this%V(n,n) )
      this%U     = arr
      DO j=1,n
         DO i=1,n
            IF ( i .eq. j ) THEN
               this%V(i,j) = mpreal('1.')
            ELSE
               this%V(i,j) = mpreal('0.')
            ENDIF
         ENDDO
      ENDDO

      DO i=1,n
         this%sigma = mpreal('0.')
      ENDDO

      RETURN
      END SUBROUTINE svd_init_plan

!======================================================================
      SUBROUTINE svd_destroy_plan(this)
!----------------------------------------------------------------------
!      Destroye an MPRSVD plan. 
!
!     Input:
!       - this: SVD plan instance.
!----------------------------------------------------------------------
      IMPLICIT NONE

      CLASS(MPRSVDPLAN), INTENT(INOUT) :: this

      DEALLOCATE(this%sweep)
      DEALLOCATE(this%sigma)
      DEALLOCATE(this%U)
      DEALLOCATE(this%V)

      RETURN
      END SUBROUTINE svd_destroy_plan

!======================================================================
      SUBROUTINE svd_output_procedure(this, unit, iostat, iomsg)
!----------------------------------------------------------------------
!       Defines the unformatted WRITE method for the derived datatype.
!----------------------------------------------------------------------

      IMPLICIT NONE

      CLASS(MPRSVDPLAN), INTENT(IN)   :: this
      INTEGER, INTENT(IN)             :: unit
      INTEGER, INTENT(OUT)            :: iostat
      CHARACTER(len=*), INTENT(INOUT) :: iomsg

      WRITE(unit, iostat=iostat, iomsg=iomsg) this%digs
      WRITE(unit, iostat=iostat, iomsg=iomsg) this%m,this%n,this%steps
      WRITE(unit, iostat=iostat, iomsg=iomsg) this%eigenmin,this%eps,this%gmax
      WRITE(unit, iostat=iostat, iomsg=iomsg) this%sweep,this%sigma 
      WRITE(unit, iostat=iostat, iomsg=iomsg) this%U,this%V

      RETURN
      END SUBROUTINE svd_output_procedure

!======================================================================
      SUBROUTINE svd_input_procedure(this, unit, iostat, iomsg)
!----------------------------------------------------------------------
!     Defines the unformatted READ method for the derived datatype.
!----------------------------------------------------------------------

      IMPLICIT NONE

      CLASS(MPRSVDPLAN), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)              :: unit
      INTEGER, INTENT(OUT)             :: iostat
      CHARACTER(len=*), INTENT(INOUT)  :: iomsg

      INTEGER :: m,n

      READ(unit, iostat=iostat, iomsg=iomsg) this%digs
      IF (this%digs .ne. mpipl ) THEN
         PRINT*, this%digs, mpipl
         PRINT*, "Attemped to load SVD plan saved with different "&
                 "numerical precision (i.e. DIGITS). Aborting..."
         STOP
      ENDIF

      READ(unit, iostat=iostat, iomsg=iomsg) this%m,this%n,this%steps
      READ(unit, iostat=iostat, iomsg=iomsg) this%eigenmin,this%eps,this%gmax

      m = this%m
      n = this%n

      ALLOCATE( this%sweep(2,n/2,n-1) )
      ALLOCATE( this%U(m,n), this%sigma(n), this%V(n,n) )

      READ(unit, iostat=iostat, iomsg=iomsg) this%sweep,this%sigma 
      READ(unit, iostat=iostat, iomsg=iomsg) this%U,this%V

      RETURN
      END SUBROUTINE svd_input_procedure


!======================================================================
      SUBROUTINE svd_save_plan(this, fname)
!----------------------------------------------------------------------
!      Save the current SVD plan to file `fname`.
!
!     Input:
!       - this:  SVD plan instance.
!       - fname: target filename to save the plan.   
!----------------------------------------------------------------------
      IMPLICIT NONE

      CLASS(MPRSVDPLAN), INTENT(IN) :: this
      CHARACTER(len=*),  INTENT(IN) :: fname

      OPEN(10, FILE=fname, FORM='unformatted', ACCESS='stream')
         WRITE(10) this
      CLOSE(10)
   
      RETURN
      END SUBROUTINE svd_save_plan

!======================================================================
      SUBROUTINE svd_load_plan(this, fname)
!----------------------------------------------------------------------
!      Load an SVD plan from file `fname`.
!
!     Input:
!       - this:  SVD plan instance.
!       - fname: target filename to save the plan.   
!----------------------------------------------------------------------
      IMPLICIT NONE

      CLASS(MPRSVDPLAN), INTENT(OUT) :: this
      CHARACTER(len=*),  INTENT(IN)  :: fname

      OPEN(10, FILE=fname, FORM='unformatted', ACCESS='stream')
         READ(10) this
      CLOSE(10)
   
      RETURN
      END SUBROUTINE svd_load_plan

!======================================================================
      SUBROUTINE svd_iterate(this,iters)
!----------------------------------------------------------------------
!      Compute the first part of a thin SVD decomposition for a matrix
!     with arbitrary precision, employing a one sided Jacobi scheme.
!     The target matrix A is defined on the call to `plan_init` as well
!     as other relevant parameters. Note that the transposed of V is
!     computed and that this subroutine iterates only the A = U'*V^t
!     decomposition. When the desired level of orthogonality is
!     achieved, A = U*sigma*V can be found calling `to_svd`.
!     Based on: https://doi.org/10.1137/0906007 and
!     https://doi.org/10.1137/0910023
!
!     Input:
!       - this:  SVD plan instance.
!       - iters: Number of iterations to perform.  
!
!     Output:
!       - plan%U: mxn matrix containing left eigenvectors and singular
!                 values of the target matrix.
!       - plan%V: nxn matrix containing right side eigenvectors of
!                 the target matrix.
!     
!      Note: The sweeping strategy requires n to be even, which 
!      always is in the generation of tables for FC-Gram.
!      For other uses, the sweeping topology must be changed.
!--------------------------------------------------------------------

      IMPLICIT NONE

      CLASS(MPRSVDPLAN), INTENT(INOUT) :: this
      INTEGER, INTENT(IN)              :: iters

      TYPE(mp_real) :: c,s,t
      TYPE(mp_real) :: alpha,beta,gama
      TYPE(mp_real) :: jgmax,igmax
      TYPE(mp_real) :: tmp

      INTEGER :: m,n
      INTEGER :: i,j,k,ii,jj
      INTEGER :: start
      INTEGER :: clock_start,clock_end,clock_rate

      m = UBOUND(this%U,1)
      n = UBOUND(this%U,2)

      start = this%steps

      CALL system_clock(clock_start,clock_rate)
      clock_end = clock_start

      ! Start iterating
      DO WHILE ( this%steps .lt. start+iters )

      ! Speed is averaged only on the current "batch" of iterations
100   FORMAT( a , 'SVD Sweep : ',i3,' of ',i3, ' (rate: ', f0.2, ' s/iter)')
      WRITE(*,100,advance="no") ACHAR(13), this%steps+1, start+iters, &
            dble(clock_end-clock_start)/((this%steps-start)*clock_rate)

      igmax = mpreal('0.')
      ! Start sweeping
      DO ii=1,n-1
         jgmax = mpreal('0.')
!$omp parallel do private(i,j,k,alpha,beta,gama,tmp,c,s,t)
         DO jj=1,n/2
            i = this%sweep(1,jj,ii)
            j = this%sweep(2,jj,ii)

            alpha = mpreal('0.')
            beta  = mpreal('0.')
            gama  = mpreal('0.')
   
            ! Solve rotation coefficients
            DO k=1,m
               alpha = alpha + this%U(k,i)**2
               beta  = beta  + this%U(k,j)**2
               gama  = gama  + this%U(k,i)*this%U(k,j)
            ENDDO

            ! Columns with inner product below eps*gmax
            ! are treated as orthogonal for the current sweep.
            IF ( ABS(gama) .le. this%eps*this%gmax ) THEN
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
               tmp = this%U(k,i)
               this%U(k,i) = c*tmp - s*this%U(k,j)
               this%U(k,j) = s*tmp + c*this%U(k,j)
            ENDDO
            DO k=1,n
               tmp = this%V(k,i)
               this%V(k,i) = c*tmp - s*this%V(k,j)
               this%V(k,j) = s*tmp + c*this%V(k,j)
            ENDDO
! Manual max reduction because OpenMP max reduction doesn't handle
! mpfun types.
!$omp critical            
            jgmax = MAX(jgmax, ABS(gama))
!$omp end critical
         ENDDO
         igmax = MAX(igmax,jgmax)
      ENDDO

      this%gmax = igmax
      CALL system_clock(clock_end)

      this%steps = this%steps+1
      ENDDO ! End of iteration

      RETURN
      END SUBROUTINE svd_iterate



!======================================================================
      SUBROUTINE svd_to_svd(this)
!----------------------------------------------------------------------
!      Computes the singular values and stores them in the vector sigma
!     while at the same time normalizing U. After calling this
!     procedure, one has A = U*sigma*V^t (i.e. the SVD decomposition).
!
!     Input:
!       - this:  SVD plan instance.
!
!     Output:
!       - this: SVD plan where ||U|| = 1 and sigma is now populated.
!--------------------------------------------------------------------

      IMPLICIT NONE

      CLASS(MPRSVDPLAN), INTENT(INOUT) :: this
      TYPE(mp_real) :: tmp

      INTEGER :: m,n
      INTEGER :: i,k

      m = UBOUND(this%U,1)
      n = UBOUND(this%U,2)

!$omp parallel do private(tmp,k)
      DO i=1,n
         tmp = mpreal('0.')
         DO k=1,m
            tmp = tmp + this%U(k,i)*this%U(k,i)
         ENDDO
         ! If actual eigenvalue, add it and normalize
         ! otherwise set to zero
         IF ( SQRT(tmp) .ge. this%eigenmin ) THEN
            this%sigma(i) = SQRT(tmp)
            this%U(:, i)  = (/( this%U(k,i)/this%sigma(i), k=1,m )/)
         ELSE
            this%sigma(i) = mpreal('0.')
            this%U(:, i)  = (/( mpreal('0.'), k=1,m )/)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE svd_to_svd

!======================================================================
      FUNCTION svd_solve(this, b) RESULT(x)
!----------------------------------------------------------------------
!      Solves the system of equations A*x = b using the SVD of A.
!     Both x and b should be mpreal vectors.
!
!     Input:
!       - this:  SVD plan instance.
!       - b:     mpreal vector.          
!
!     Output:
!       - x: The solution to the system of equations.
!
!     TODO: implement A*X = B (solver for matrices)
!--------------------------------------------------------------------

      IMPLICIT NONE

      CLASS(MPRSVDPLAN), INTENT(IN)                          :: this
      TYPE(mp_real), DIMENSION(UBOUND(this%U,1)), INTENT(IN) :: b
      TYPE(mp_real), DIMENSION(UBOUND(this%U,2))             :: x

      INTEGER :: i

      ! Compute D*V^t*x  ( = U^t*b)
      x = MPRMATVEC(MPRTRANSPOSE(this%U), b)

      ! Compute V^t*x. (i.e. divide x_i by sigma_i)
!$omp parallel do
      DO i = 1,UBOUND(x,1)
         IF ( this%sigma(i) >= this%get_eigenmin() ) THEN
             x(i) = x(i)/this%sigma(i)
         ENDIF
      ENDDO

      ! Compute x
      x = MPRMATVEC(this%V, x)

      RETURN
      END FUNCTION svd_solve


!======================================================================
      FUNCTION svd_get_eigenmin(this) RESULT(out)
!----------------------------------------------------------------------
!      Returns the minimum singular value to consider in the SVD
!     decomposition. Singular values under eigenmin are considered 0.
!
!     Input:
!       - this: SVD plan
!
!     Output:
!       - out: Current minimum singular value to consider non 0
!
!     TODO: Sanity checks
!----------------------------------------------------------------------
      IMPLICIT NONE


      CLASS(MPRSVDPLAN), INTENT(IN) :: this
      TYPE(mp_real)                 :: out

      out = this%eigenmin
    
      RETURN
      END FUNCTION svd_get_eigenmin

!======================================================================
      FUNCTION svd_set_eigenmin(this, val) RESULT(stat)
!----------------------------------------------------------------------
!      Sets the minimum singular value to consider in the SVD
!     decomposition. Singular values under eigenmin are considered 0.
!
!     Input:
!       - this: SVD plan
!       - val: desired new minimum singular value to consider (MPREAL)
!
!     Output:
!       - this: SVD plan
!
!     TODO: Sanity checks
!----------------------------------------------------------------------
      IMPLICIT NONE


      CLASS(MPRSVDPLAN), INTENT(INOUT) :: this
      TYPE(mp_real), INTENT(IN)        :: val
      INTEGER                          :: stat

      this%eigenmin = val

      stat = 0   ! TODO change return value according to success or not
    
      RETURN
      END FUNCTION svd_set_eigenmin


!=====================================================================
      FUNCTION check_orthogonal(arr) RESULT(maxin)
!---------------------------------------------------------------------
!     Checks orthogonality between columns of MPREAL arrays.
!
!     Input:
!       - arr:  mpreal array.
!
!     Output:
!       - maxin: maximum inner product (mpreal).
!
!--------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(mp_real), DIMENSION(:,:), INTENT(IN) :: arr
      TYPE(mp_real) :: maxin

      TYPE(mp_real) :: locmax

      TYPE(mp_real) :: tmp
      INTEGER :: i,j,k

      maxin  = mpreal('0.')
      locmax = mpreal('0.')
!omp parallel do private(j,k,tmp,locmax)
      DO i=1,UBOUND(arr,2)
      DO j=i+1,UBOUND(arr,2)
         tmp = mpreal('0.')
         DO k=1,UBOUND(arr,1)
            tmp = tmp + arr(k,i)*arr(k,j)
         ENDDO
         locmax = MAX(locmax,tmp) 
      ENDDO
! Manual max reduction because OpenMP max reduction doesn't handle
! mpfun types.
!$omp critical            
         maxin = MAX(maxin, locmax)
!$omp end critical
      ENDDO

      RETURN
      END FUNCTION check_orthogonal

      END MODULE mprlinalg
