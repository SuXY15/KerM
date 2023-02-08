!=======================================================================
!  BLAS module library
!  Version 1.1
!
!  Contains only procedures used by PDFD code
!=======================================================================
!#1  subroutine daxpy
!#2  subroutine dcopy
!#3  function ddot
!#4  subroutine dscal
!#5  subroutine xerbla


!=======================================================================
!  LAPACK module, constructed from fortran 77 source code
!  99.04.05
!  Version 1.1
!
!  Contains only procedures used by PDFD code
!=======================================================================
!#6  subroutine dgemm
!#7  subroutine dgemv
!#8  subroutine dpbsv
!#9  subroutine dpbtf2
!#10 subroutine dpbtrf
!#11 subroutine dpbtrs
!#12 subroutine dpotf2
!#13 subroutine dsyr
!#14 subroutine dsyrk
!#15 subroutine dtbsv
!#16 subroutine dtrsm
!#17 function ilaenv
!#18 function lsame



      MODULE BLAS_LAPACK_MODULE


      REAL(KIND=8), PARAMETER :: ZERO = 0.0D0
      REAL(KIND=8), PARAMETER :: ONE  = 1.0D0
      REAL(KIND=8), PARAMETER :: TWO  = 2.0D0
      REAL(KIND=8), PARAMETER :: HALF = 0.5D0


      CONTAINS


!=======================================================================
!  BLAS module library
!=======================================================================


!#1
!***********************************************************************
      SUBROUTINE daxpy (n, da, dx, incx, dy, incy)
!***********************************************************************
!                b l a s  subprogram
!  Purpose: computation of y = a*x + y
!-----------------------------------------------------------------------
!  description of parameters
!
!  --input--
!        n  number of elements in input vector(s)
!       da  double precision scalar multiplier
!       dx  double precision vector with n elements
!     incx  storage spacing between elements of dx
!       dy  double precision vector with n elements
!     incy  storage spacing between elements of dy
!
!  --output--
!       dy  double precision result (unchanged if n .le. 0)
!-----------------------------------------------------------------------
!  Overwrite double precision dy with double precision da*dx + dy.
!  For i = 0 to n-1, replace  dy(ly+i*incy) with da*dx(lx+i*incx) +
!  dy(ly+i*incy), where lx = 1 if incx .ge. 0, else lx = (-incx)*n
!  and ly is defined in a similar way using incy.
!-----------------------------------------------------------------------
!  References: Lawson, C.L., Hanson, R.J., Kincaid, D.R., Krogh, F.T.,
!              "basic linear algebra subprograms for fortran usage",
!              Algorithm no. 539, transactions on mathematical
!              software, volume 5, number 3, September 1979, 308-323
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER (KIND=4), INTENT(IN)  :: n, incx, incy
      REAL    (KIND=8), INTENT(IN)  :: dx(*), da
      REAL    (KIND=8), INTENT(OUT) :: dy(*)

      INTEGER (KIND=4) :: m, ix, iy, i, ns

      IF (n >=0 .AND. da /= zero) THEN

!--------  code for nonequal or nonpositive increments.
         IF (incx /= incy .OR. incx < 1) THEN
            ix = 1
            iy = 1
            IF (incx < 0) ix = (-n+1)*incx + 1
            IF (incy < 0) iy = (-n+1)*incy + 1
            DO i = 1, n
               dy(iy) = dy(iy) + da*dx(ix)
               ix = ix + incx
               iy = iy + incy
            END DO

!--------   code for both increments equal to 1
         ELSE IF (incx == 1 .AND. incy == 1) THEN

!...........clean-up loop so remaining vector length is a multiple of 4.

            m = MOD(n,4)
            IF (m /= 0) THEN
               DO i = 1, m
                  dy(i) = dy(i) + da*dx(i)
               END DO
               IF (n < 4) RETURN
            END IF
            IF (n >= 4) THEN
               DO i = m+1, n, 4
                  dy(i)   = dy(i)   + da*dx(i)
                  dy(i+1) = dy(i+1) + da*dx(i+1)
                  dy(i+2) = dy(i+2) + da*dx(i+2)
                  dy(i+3) = dy(i+3) + da*dx(i+3)
               END DO
            END IF

!--------  code for equal, positive, nonunit increments.
         ELSE
            ns = n*incx
            DO i = 1, ns, incx
               dy(i) = da*dx(i) + dy(i)
            END DO
         END IF
      END IF

      END SUBROUTINE daxpy


!#2
!***********************************************************************
      SUBROUTINE dcopy (n, dx, incx, dy, incy)
!***********************************************************************
!                b l a s  subprogram
!  Purpose: vector copy y = x
!-----------------------------------------------------------------------
!  description of parameters
!
!  --input--
!        n  number of elements in input vector(s)
!       dx  double precision vector with n elements
!     incx  storage spacing between elements of dx
!       dy  double precision vector with n elements
!     incy  storage spacing between elements of dy
!
!  --output--
!       dy  copy of vector dx (unchanged if n .le. 0)
!-----------------------------------------------------------------------
!  Copy dx to dy.
!  For i = 0 to n-1, copy dx(lx+i*incx) to dy(ly+i*incy),
!  where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
!  defined in a similar way using incy.
!-----------------------------------------------------------------------
!  References: Lawson, C.L., Hanson, R.J., Kincaid, D.R., Krogh, F.T.,
!              "basic linear algebra subprograms for fortran usage",
!              Algorithm no. 539, transactions on mathematical
!              software, volume 5, number 3, September 1979, 308-323
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER (KIND=4), INTENT(IN)    :: n, incx, incy
      REAL    (KIND=8), INTENT(IN)    :: dx(*)
      REAL    (KIND=8), INTENT(INOUT) :: dy(*)

      INTEGER (KIND=4) :: ix, iy, i, m, ns

      IF (n > 0) THEN
 

!--------  code for both increments equal to 1

         IF (incx == incy .AND. incx == 1) THEN

!...........clean-up loop so remaining vector length is a multiple of 7.
            m = MOD (n,7)
            IF (m > 0) THEN
               DO i = 1, m
                  dy(i) = dx(i)
               END DO
            END IF
            IF (n > 6) THEN
               DO i = m+1, n, 7
                  dy(i)   = dx(i)
                  dy(i+1) = dx(i+1)
                  dy(i+2) = dx(i+2)
                  dy(i+3) = dx(i+3)
                  dy(i+4) = dx(i+4)
                  dy(i+5) = dx(i+5)
                  dy(i+6) = dx(i+6)
               END DO
            END IF

!--------  code for equal, positive, nonunit increments.

         ELSE IF (incx == incy .AND. incx > 1) THEN

            ns = n*incx
            DO i = 1, ns, incx
               dy(i) = dx(i)
            END DO

!--------  code for unequal or nonpositive increments.
 
         ELSE
 
            ix = 1
            iy = 1
            IF (incx < 0) ix = (-n+1)*incx + 1
            IF (incy < 0) iy = (-n+1)*incy + 1
            DO i = 1, n
              dy(iy) = dx(ix)
              ix = ix + incx
              iy = iy + incy
            END DO
      
         END IF
      END IF

      END SUBROUTINE dcopy


!#3
!***********************************************************************
      FUNCTION ddot(n, dx, incx, dy, incy)
!***********************************************************************
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
!=======================================================================

      IMPLICIT NONE

      REAL    (KIND=8) :: ddot
      INTEGER (KIND=4), INTENT(IN) :: incx, incy, n
      REAL    (KIND=8), INTENT(IN) :: dx(*), dy(*)

      INTEGER (KIND=4) :: i, ix, iy, m, mp1
      REAL    (KIND=8) :: dtemp

      ddot  = zero
      dtemp = zero
      IF (n <= 0) RETURN

      IF (incx == 1 .AND. incy == 1) THEN
!
!        code for both increments equal to 1
!
         m = MOD(n,5)
         IF (m /= 0) THEN
            DO i = 1, m
               dtemp = dtemp + dx(i)*dy(i)
            END DO
         END IF
         IF (n >= 5) THEN
            mp1 = m + 1
            DO i = mp1, n, 5
               dtemp = dtemp + dx(i)  *dy(i)   + &
                               dx(i+1)*dy(i+1) + &
                               dx(i+2)*dy(i+2) + &
                               dx(i+3)*dy(i+3) + &
                               dx(i+4)*dy(i+4)
            END DO
         END IF

      ELSE
!
!        code for unequal increments or equal increments not equal to 1
!
         ix = 1
         iy = 1
         IF (incx < 0) ix = (-n+1)*incx + 1
         IF (incy < 0) iy = (-n+1)*incy + 1
         DO i = 1, n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF

      ddot = dtemp

      END FUNCTION ddot


!#4
!***********************************************************************
      SUBROUTINE dscal(n, da, dx, incx)
!***********************************************************************
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!
!=======================================================================

      IMPLICIT NONE

      INTEGER (KIND=4), INTENT(IN)    :: incx, n
      REAL    (KIND=8), INTENT(IN)    :: da
      REAL    (KIND=8), INTENT(INOUT) :: dx(*)

      INTEGER (KIND=4) :: i, m, mp1, nincx

      IF ( n <= 0 .OR. incx <= 0 ) RETURN

      IF (incx == 1) THEN
!
!        code for increment equal to 1
!
         m = MOD(n,5)
         IF ( m /= 0 ) THEN
            DO i = 1, m
               dx(i) = da*dx(i)
            END DO
         END IF
         IF ( n >= 5 ) THEN
            mp1 = m + 1
            DO i = mp1, n, 5
               dx(i)   = da*dx(i)
               dx(i+1) = da*dx(i+1)
               dx(i+2) = da*dx(i+2)
               dx(i+3) = da*dx(i+3)
               dx(i+4) = da*dx(i+4)
            END DO
         END IF

      ELSE
!
!        code for increment not equal to 1
!
         nincx = n*incx
         DO i = 1, nincx, incx
            dx(i) = da*dx(i)
         END DO
      END IF

   END SUBROUTINE dscal


!#5
!***********************************************************************
   SUBROUTINE XERBLA( SRNAME, INFO )
!***********************************************************************
!
!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER(LEN=*)
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
!
!=======================================================================

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=*), INTENT(IN) :: SRNAME
      INTEGER (KIND=4), INTENT(IN) :: INFO

      WRITE( *, '(A,A,A,I2,A)') &
         ' ** On entry to ', TRIM(SRNAME), ' parameter number ', INFO, &
         ' had an illegal value'

      STOP

   END SUBROUTINE XERBLA



!=======================================================================
!  LAPACK module library
!=======================================================================


!#6
!*****************************************************************************
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
                         BETA, C, LDC )
!*****************************************************************************
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!========================================================================

!      USE BLAS_MODULE, ONLY : XERBLA

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN) :: TRANSA, TRANSB
      INTEGER (KIND=4), INTENT(IN) :: M, N, K, LDA, LDB, LDC
      REAL    (KIND=8), INTENT(IN) :: ALPHA, BETA

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(IN)    :: A( LDA, * )
      REAL    (KIND=8), INTENT(IN)    :: B( LDB, * )
      REAL    (KIND=8), INTENT(INOUT) :: C( LDC, * )

!     .. External Functions ..
!     LOGICAL :: LSAME
!     EXTERNAL LSAME

!     .. External Subroutines ..
!     EXTERNAL XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC MAX

!     .. Local Scalars ..
      LOGICAL :: NOTA, NOTB
      INTEGER (KIND=4) :: I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL    (KIND=8) :: TEMP

!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA = LSAME( TRANSA, 'N' )
      NOTB = LSAME( TRANSB, 'N' )
      IF ( NOTA ) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF ( NOTB ) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF (      ( .NOT.NOTA                 ).AND. &
                ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
                ( .NOT.LSAME( TRANSA, 'T' ) )      ) THEN
         INFO = 1
      ELSE IF ( ( .NOT.NOTB                 ).AND. &
                ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
                ( .NOT.LSAME( TRANSB, 'T' ) )      ) THEN
         INFO = 2
      ELSE IF ( M   < 0               ) THEN
         INFO = 3
      ELSE IF ( N   < 0               ) THEN
         INFO = 4
      ELSE IF ( K   < 0               ) THEN
         INFO = 5
      ELSE IF ( LDA < MAX( 1, NROWA ) ) THEN
         INFO = 8
      ELSE IF ( LDB < MAX( 1, NROWB ) ) THEN
         INFO = 10
      ELSE IF ( LDC < MAX( 1, M     ) ) THEN
         INFO = 13
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF ( ( M == 0 ).OR.( N == 0 ).OR. &
          ( ( ( ALPHA == ZERO ).OR.( K == 0 ) ).AND.( BETA == ONE ) ) ) THEN
         RETURN
      END IF
!
!     And if  alpha.eq.zero.
!
      IF ( ALPHA == ZERO ) THEN
         IF ( BETA == ZERO ) THEN
            DO J = 1, N
               DO I = 1, M
                  C( I, J ) = ZERO
               END DO
            END DO
         ELSE
            DO J = 1, N
               DO I = 1, M
                  C( I, J ) = BETA*C( I, J )
               END DO
            END DO
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF ( NOTB ) THEN
         IF ( NOTA ) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO J = 1, N
               IF ( BETA == ZERO ) THEN
                  DO I = 1, M
                     C( I, J ) = ZERO
                  END DO
               ELSE IF ( BETA /= ONE ) THEN
                  DO I = 1, M
                     C( I, J ) = BETA*C( I, J )
                  END DO
               END IF
               DO L = 1, K
                  IF ( B( L, J ) /= ZERO ) THEN
                     TEMP = ALPHA*B( L, J )
                     DO I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
                     END DO
                  END IF
               END DO
            END DO
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO J = 1, N
               DO I = 1, M
                  TEMP = ZERO
                  DO L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
                  END DO
                  IF ( BETA == ZERO ) THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
               END DO
            END DO
         END IF
      ELSE
         IF ( NOTA ) THEN
!
!           Form  C := alpha*A*B' + beta*C
!
            DO J = 1, N
               IF ( BETA == ZERO ) THEN
                  DO I = 1, M
                     C( I, J ) = ZERO
                  END DO
               ELSE IF ( BETA /= ONE ) THEN
                  DO I = 1, M
                     C( I, J ) = BETA*C( I, J )
                  END DO
               END IF
               DO L = 1, K
                  IF ( B( J, L ) /= ZERO ) THEN
                     TEMP = ALPHA*B( J, L )
                     DO I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
                     END DO
                  END IF
               END DO
            END DO
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO J = 1, N
               DO I = 1, M
                  TEMP = ZERO
                  DO L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
                  END DO
                  IF ( BETA == ZERO ) THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
               END DO
            END DO
         END IF
      END IF

      END SUBROUTINE DGEMM


!#7
!*****************************************************************************
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
                         BETA, Y, INCY )
!*****************************************************************************
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!===========================================================================

!      USE BLAS_MODULE, ONLY : XERBLA

      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER (KIND=4), INTENT(IN) :: INCX, INCY, LDA, M, N
      REAL    (KIND=8), INTENT(IN) :: ALPHA, BETA
      CHARACTER(LEN=1), INTENT(IN) :: TRANS

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(IN)    :: A( LDA, * ), X( * )
      REAL    (KIND=8), INTENT(INOUT) :: Y( * )

!     .. Local Scalars ..
      REAL    (KIND=8) :: TEMP
      INTEGER (KIND=4) :: I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY

!     .. External Functions ..
!     LOGICAL :: LSAME
!     EXTERNAL LSAME

!     .. External Subroutines ..
!     EXTERNAL XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC MAX
!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF      ( .NOT.LSAME( TRANS, 'N' ).AND. &
                .NOT.LSAME( TRANS, 'T' ).AND. &
                .NOT.LSAME( TRANS, 'C' )      ) THEN
         INFO = 1
      ELSE IF ( M < 0 ) THEN
         INFO = 2
      ELSE IF ( N < 0 ) THEN
         INFO = 3
      ELSE IF ( LDA < MAX( 1, M ) ) THEN
         INFO = 6
      ELSE IF ( INCX == 0 ) THEN
         INFO = 8
      ELSE IF ( INCY == 0 ) THEN
         INFO = 11
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF ( ( M == 0 ).OR.( N == 0 ).OR. &
         ( ( ALPHA == ZERO ).AND.( BETA == ONE ) ) ) THEN
         RETURN
      END IF
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF ( LSAME( TRANS, 'N' ) ) THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF ( INCX > 0 ) THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF ( INCY > 0 ) THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF ( BETA /= ONE ) THEN
         IF ( INCY == 1 ) THEN
            IF ( BETA == ZERO ) THEN
               DO I = 1, LENY
                  Y( I ) = ZERO
               END DO
            ELSE
               DO I = 1, LENY
                  Y( I ) = BETA*Y( I )
               END DO
            END IF
         ELSE
            IY = KY
            IF ( BETA == ZERO ) THEN
               DO I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
               END DO
            ELSE
               DO I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
               END DO
            END IF
         END IF
      END IF
      IF ( ALPHA == ZERO ) RETURN
      IF ( LSAME( TRANS, 'N' ) ) THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF ( INCY == 1 ) THEN
            DO J = 1, N
               IF ( X( JX ) /= ZERO ) THEN
                  TEMP = ALPHA*X( JX )
                  DO I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
                  END DO
               END IF
               JX = JX + INCX
            END DO
         ELSE
            DO J = 1, N
               IF ( X( JX ) /= ZERO ) THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
                  END DO
               END IF
               JX = JX + INCX
            END DO
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY
         IF ( INCX == 1 ) THEN
            DO J = 1, N
               TEMP = ZERO
               DO I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
               END DO
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
            END DO
         ELSE
            DO J = 1, N
               TEMP = ZERO
               IX   = KX
               DO I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
               END DO
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
            END DO
         END IF
      END IF

      END SUBROUTINE DGEMV


!#8
!*****************************************************************************
      SUBROUTINE DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!*****************************************************************************
!
!  -- LAPACK driver routine (version 1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!
!  Purpose
!  =======
!
!  DPBSV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N-by-N symmetric positive definite band matrix and X
!  and B are N-by-NRHS matrices.
!
!  The Cholesky decomposition is used to factor A as
!     A = U**T * U,  if UPLO = 'U', or
!     A = L * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix, and L is a lower triangular
!  matrix, with the same number of superdiagonals or subdiagonals as
!  A.  The factored form of A is then used to solve the system of
!  equations A * X = B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).
!          See below for further details.
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U!*T*U or A = L*L**T of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N-by-NRHS right hand side matrix B.
!          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i of A is not
!                positive definite, so the factorization could not be
!                completed, and the solution has not been computed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked ! are not used by the routine.
!
!=======================================================================

!      USE BLAS_MODULE, ONLY : XERBLA
      
      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN)  :: UPLO
      INTEGER (KIND=4), INTENT(IN)  :: KD, LDAB, LDB, N, NRHS
      INTEGER (KIND=4), INTENT(OUT) :: INFO

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(INOUT) :: AB( LDAB, * ), B( LDB, * )

!     .. Intrinsic Functions ..
      INTRINSIC MAX

!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF ( N < 0 ) THEN
         INFO = -2
      ELSE IF ( KD < 0 ) THEN
         INFO = -3
      ELSE IF ( NRHS < 0 ) THEN
         INFO = -4
      ELSE IF ( LDAB < KD+1 ) THEN
         INFO = -6
      ELSE IF ( LDB < MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DPBSV ', -INFO )
         RETURN
      END IF
!
!     Compute the Cholesky factorization A = U'*U or A = L*L'.
!
      CALL DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
      IF ( INFO == 0 ) THEN
!
!        Solve the system A*X = B, overwriting B with X.
!
         CALL DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!
      END IF

      END SUBROUTINE DPBSV


!#9
!*****************************************************************************
      SUBROUTINE DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
!*****************************************************************************
!
!  -- LAPACK routine (version 1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!  Purpose
!  =======
!
!  DPBTF2 computes the Cholesky factorization of a real symmetric
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix, U' is the transpose of U, and
!  L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of super-diagonals of the matrix A if UPLO = 'U',
!          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U'*U or A = L*L' of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  =====================================================================

!      USE BLAS_MODULE, ONLY : DSCAL, XERBLA

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN)  :: UPLO
      INTEGER (KIND=4), INTENT(IN)  :: KD, LDAB, N
      INTEGER (KIND=4), INTENT(OUT) :: INFO

!     .. Array Arguments ..
      REAL(KIND=8), INTENT(INOUT) :: AB( LDAB, * )

!     .. Local Scalars ..
      INTEGER (KIND=4) :: J, KLD, KN
      REAL    (KIND=8) :: AJJ
      LOGICAL :: UPPER

!     .. External Functions ..
!     LOGICAL :: LSAME
!     EXTERNAL LSAME

!     .. External Subroutines ..
!     EXTERNAL DSCAL, DSYR, XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC MAX, MIN, SQRT
!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF ( N < 0 ) THEN
         INFO = -2
      ELSE IF ( KD < 0 ) THEN
         INFO = -3
      ELSE IF ( LDAB < KD+1 ) THEN
         INFO = -5
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DPBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF ( N == 0 ) RETURN
!
      KLD = MAX( 1, LDAB-1 )
!
      IF ( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         DO J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = AB( KD+1, J )
            IF ( AJJ <= ZERO ) THEN
               INFO = J
               EXIT
            END IF
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
!
!           Compute elements J+1:J+KN of row J and update the
!           trailing submatrix within the band.
!
            KN = MIN( KD, N-J )
            IF ( KN > 0 ) THEN
               CALL DSCAL( KN, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL DSYR( 'Upper', KN, -ONE, AB( KD, J+1 ), KLD, &
                          AB( KD+1, J+1 ), KLD )
            END IF
         END DO
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         DO J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = AB( 1, J )
            IF ( AJJ <= ZERO ) THEN
               INFO = J
               EXIT
            END IF
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
!
!           Compute elements J+1:J+KN of column J and update the
!           trailing submatrix within the band.
!
            KN = MIN( KD, N-J )
            IF ( KN > 0 ) THEN
               CALL DSCAL( KN, ONE / AJJ, AB( 2, J ), 1 )
               CALL DSYR( 'Lower', KN, -ONE, AB( 2, J ), 1, &
                          AB( 1, J+1 ), KLD )
            END IF
         END DO
      END IF

      END SUBROUTINE DPBTF2


!#10
!*****************************************************************************
      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
!*****************************************************************************
!
!  -- LAPACK routine (version 1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!
!  Purpose
!  =======
!
!  DPBTRF computes the Cholesky factorization of a real symmetric
!  positive definite band matrix A.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER!1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
!          On entry, the upper or lower triangle of the symmetric band
!          matrix A, stored in the first KD+1 rows of the array.  The
!          j-th column of A is stored in the j-th column of the array AB
!          as follows:
!          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T of the band
!          matrix A, in the same storage format as A.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  N = 6, KD = 2, and UPLO = 'U':
!
!  On entry:                       On exit:
!
!      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!
!  Similarly, if UPLO = 'L' the format of A is as follows:
!
!  On entry:                       On exit:
!
!     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
!     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
!     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
!
!  Array elements marked * are not used by the routine.
!
!  Contributed by
!  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
!
!=======================================================================

!      USE BLAS_MODULE, ONLY : XERBLA

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN)  :: UPLO
      INTEGER (KIND=4), INTENT(IN)  :: KD, LDAB, N
      INTEGER (KIND=4), INTENT(OUT) :: INFO

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(INOUT) :: AB( LDAB, * )
!     ..
!     .. Parameters ..
      INTEGER (KIND=4),PARAMETER :: NBMAX = 32, LDWORK = NBMAX+1

!     .. Local Scalars ..
      INTEGER (KIND=4) :: I, I2, I3, IB, II, J, JJ, NB

!     .. Local Arrays ..
      REAL    (KIND=8) :: WORK( LDWORK, NBMAX )

!     .. External Functions ..
!     LOGICAL :: LSAME
!     INTEGER :: ILAENV
!     EXTERNAL LSAME, ILAENV

!     .. External Subroutines ..
!     EXTERNAL DGEMM, DPBTF2, DPOTF2, DSYRK, DTRSM, XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC MIN
!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF ( ( .NOT.LSAME( UPLO, 'U' ) ) .AND. &
           ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -1
      ELSE IF ( N < 0 ) THEN
         INFO = -2
      ELSE IF ( KD < 0 ) THEN
         INFO = -3
      ELSE IF ( LDAB < KD+1 ) THEN
         INFO = -5
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DPBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF ( N == 0 ) RETURN
!
!     Determine the block size for this environment
!
!     NB = ILAENV( 1, 'DPBTRF', UPLO, N, KD, -1, -1 )
      NB = ILAENV( 1, 'DPBTRF', N, KD, -1 )
!
!     The block size must not exceed the semi-bandwidth KD, and must not
!     exceed the limit set by the size of the local array WORK.
!
      NB = MIN( NB, NBMAX )
!
      IF ( NB <= 1 .OR. NB > KD ) THEN
!
!        Use unblocked code
!
         CALL DPBTF2( UPLO, N, KD, AB, LDAB, INFO )
      ELSE
!
!        Use blocked code
!
         IF ( LSAME( UPLO, 'U' ) ) THEN
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the upper triangle of the matrix in band
!           storage.
!
!           Zero the upper triangle of the work array.
!
            DO J = 1, NB
               DO I = 1, J - 1
                  WORK( I, J ) = ZERO
               END DO
            END DO
!
!           Process the band matrix one diagonal block at a time.
!
            DO I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL DPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
               IF ( II /= 0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF ( I+IB <= N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11   A12   A13
!                          A22   A23
!                                A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A12, A22 and
!                 A23 are empty if IB = KD. The upper triangle of A13
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF ( I2 > 0 ) THEN
!
!                    Update A12
!
                     CALL DTRSM( 'Left', 'Upper', 'Transpose', &
                                 'Non-unit', IB, I2, ONE, AB( KD+1, I ), &
                                 LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )
!
!                    Update A22
!
                     CALL DSYRK( 'Upper', 'Transpose', I2, IB, -ONE, &
                                 AB( KD+1-IB, I+IB ), LDAB-1, ONE, &
                                 AB( KD+1, I+IB ), LDAB-1 )
                  END IF
!
                  IF ( I3 > 0 ) THEN
!
!                    Copy the lower triangle of A13 into the work array.
!
                     DO JJ = 1, I3
                        DO II = JJ, IB
                           WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
                        END DO
                     END DO
!
!                    Update A13 (in the work array).
!
                     CALL DTRSM( 'Left', 'Upper', 'Transpose', &
                                 'Non-unit', IB, I3, ONE, AB( KD+1, I ), &
                                 LDAB-1, WORK, LDWORK )
!
!                    Update A23
!
                     IF ( I2 > 0 ) THEN
                        CALL DGEMM( 'Transpose', 'No Transpose', I2, I3, &
                                    IB, -ONE, AB( KD+1-IB, I+IB ), &
                                    LDAB-1, WORK, LDWORK, ONE, &
                                    AB( 1+IB, I+KD ), LDAB-1 )
                     END IF
!
!                    Update A33
!
                     CALL DSYRK( 'Upper', 'Transpose', I3, IB, -ONE, &
                                 WORK, LDWORK, ONE, AB( KD+1, I+KD ), &
                                 LDAB-1 )
!
!                    Copy the lower triangle of A13 back into place.
!
                     DO JJ = 1, I3
                        DO II = JJ, IB
                           AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
                        END DO
                     END DO
                  END IF
               END IF
            END DO
         ELSE
!
!           Compute the Cholesky factorization of a symmetric band
!           matrix, given the lower triangle of the matrix in band
!           storage.
!
!           Zero the lower triangle of the work array.
!
            DO J = 1, NB
               DO I = J + 1, NB
                  WORK( I, J ) = ZERO
               END DO
            END DO
!
!           Process the band matrix one diagonal block at a time.
!
            DO I = 1, N, NB
               IB = MIN( NB, N-I+1 )
!
!              Factorize the diagonal block
!
               CALL DPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
               IF ( II /= 0 ) THEN
                  INFO = I + II - 1
                  GO TO 150
               END IF
               IF ( I+IB <= N ) THEN
!
!                 Update the relevant part of the trailing submatrix.
!                 If A11 denotes the diagonal block which has just been
!                 factorized, then we need to update the remaining
!                 blocks in the diagram:
!
!                    A11
!                    A21   A22
!                    A31   A32   A33
!
!                 The numbers of rows and columns in the partitioning
!                 are IB, I2, I3 respectively. The blocks A21, A22 and
!                 A32 are empty if IB = KD. The lower triangle of A31
!                 lies outside the band.
!
                  I2 = MIN( KD-IB, N-I-IB+1 )
                  I3 = MIN( IB, N-I-KD+1 )
!
                  IF ( I2 > 0 ) THEN
!
!                    Update A21
!
                     CALL DTRSM( 'Right', 'Lower', 'Transpose', &
                                 'Non-unit', I2, IB, ONE, AB( 1, I ), &
                                 LDAB-1, AB( 1+IB, I ), LDAB-1 )
!
!                    Update A22
!
                     CALL DSYRK( 'Lower', 'No Transpose', I2, IB, -ONE, &
                                 AB( 1+IB, I ), LDAB-1, ONE, &
                                 AB( 1, I+IB ), LDAB-1 )
                  END IF
!
                  IF ( I3 > 0 ) THEN
!
!                    Copy the upper triangle of A31 into the work array.
!
                     DO JJ = 1, IB
                        DO II = 1, MIN( JJ, I3 )
                           WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
                        END DO
                     END DO
!
!                    Update A31 (in the work array).
!
                     CALL DTRSM( 'Right', 'Lower', 'Transpose', &
                                 'Non-unit', I3, IB, ONE, AB( 1, I ), &
                                 LDAB-1, WORK, LDWORK )
!
!                    Update A32
!
                     IF ( I2 > 0 ) THEN
                        CALL DGEMM( 'No transpose', 'Transpose', I3, I2, &
                                    IB, -ONE, WORK, LDWORK, &
                                    AB( 1+IB, I ), LDAB-1, ONE, &
                                    AB( 1+KD-IB, I+IB ), LDAB-1 )
                     END IF
!
!                    Update A33
!
                     CALL DSYRK( 'Lower', 'No Transpose', I3, IB, -ONE, &
                                 WORK, LDWORK, ONE, AB( 1, I+KD ), &
                                 LDAB-1 )
!
!                    Copy the upper triangle of A31 back into place.
!
                     DO JJ = 1, IB
                        DO II = 1, MIN( JJ, I3 )
                           AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
                        END DO
                     END DO
                  END IF
               END IF
            END DO
         END IF
      END IF
!
  150 CONTINUE

      END SUBROUTINE DPBTRF


!#11
!*****************************************************************************
      SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
!*****************************************************************************
!
!  -- LAPACK routine (version 1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!  Purpose
!  =======
!
!  DPBTRS solves a system of linear equations A*X = B with a symmetric
!  positive definite band matrix A using the Cholesky factorization
!  A = U**T*U or A = L*L**T computed by DPBTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangular factor stored in AB;
!          = 'L':  Lower triangular factor stored in AB.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KD      (input) INTEGER
!          The number of superdiagonals of the matrix A if UPLO = 'U',
!          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
!          The triangular factor U or L from the Cholesky factorization
!          A = U**T*U or A = L*L**T of the band matrix A, stored in the
!          first KD+1 rows of the array.  The j-th column of U or L is
!          stored in the array AB as follows:
!          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
!          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= KD+1.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================

!      USE BLAS_MODULE, ONLY : XERBLA

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN)  :: UPLO
      INTEGER (KIND=4), INTENT(IN)  :: KD, LDAB, LDB, N, NRHS
      INTEGER (KIND=4), INTENT(OUT) :: INFO
!     ..
!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(IN)    :: AB( LDAB, * )
      REAL    (KIND=8), INTENT(INOUT) :: B( LDB, * )
!     ..
!     .. Local Scalars ..
      LOGICAL :: UPPER
      INTEGER (KIND=4) :: J
!     ..
!     .. External Functions ..
!     LOGICAL :: LSAME
!     EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
!     EXTERNAL DTBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF ( N < 0 ) THEN
         INFO = -2
      ELSE IF ( KD < 0 ) THEN
         INFO = -3
      ELSE IF ( NRHS < 0 ) THEN
         INFO = -4
      ELSE IF ( LDAB < KD+1 ) THEN
         INFO = -6
      ELSE IF ( LDB < MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DPBTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF ( N == 0 .OR. NRHS == 0 ) RETURN
!
      IF ( UPPER ) THEN
!
!        Solve A*X = B where A = U'*U.
!
         DO J = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
!
!           Solve U*X = B, overwriting B with X.
!
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
         END DO
      ELSE
!
!        Solve A*X = B where A = L*L'.
!
         DO J = 1, NRHS
!
!           Solve L*X = B, overwriting B with X.
!
            CALL DTBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
!
!           Solve L'*X = B, overwriting B with X.
!
            CALL DTBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB, &
                        LDAB, B( 1, J ), 1 )
         END DO
      END IF

      END SUBROUTINE DPBTRS


!#12
!*****************************************************************************
      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO )
!*****************************************************************************
!
!  -- LAPACK routine (version 1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!
!  Purpose
!  =======
!
!  DPOTF2 computes the Cholesky factorization of a real symmetric
!  positive definite matrix A.
!
!  The factorization has the form
!     A = U' * U ,  if UPLO = 'U', or
!     A = L  * L',  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n by n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n by n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, if INFO = 0, the factor U or L from the Cholesky
!          factorization A = U'*U  or A = L*L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite, and the factorization could not be
!               completed.
!
!  =====================================================================

!      USE BLAS_MODULE, ONLY : DDOT, DSCAL, XERBLA

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN)  :: UPLO
      INTEGER (KIND=4), INTENT(IN)  :: LDA, N
      INTEGER (KIND=4), INTENT(OUT) :: INFO

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(INOUT) :: A( LDA, * )

!     .. Local Scalars ..
      INTEGER (KIND=4) :: J
      REAL    (KIND=8) :: AJJ
      LOGICAL :: UPPER

!     .. External Functions ..
!     LOGICAL :: LSAME
!     REAL    :: DDOT
!     EXTERNAL LSAME, DDOT

!     .. External Subroutines ..
!     EXTERNAL DGEMV, DSCAL, XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC MAX, SQRT
!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF ( N < 0 ) THEN
         INFO = -2
      ELSE IF ( LDA < MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DPOTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF ( N == 0 ) RETURN

      IF ( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U'*U.
!
         DO J = 1, N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF ( AJJ <= ZERO ) THEN
               A( J, J ) = AJJ
               INFO = J
               EXIT
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of row J.
!
            IF ( J < N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ), &
                           LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
         END DO
      ELSE
!
!        Compute the Cholesky factorization A = L*L'.
!
         DO J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), LDA )
            IF ( AJJ <= ZERO ) THEN
               A( J, J ) = AJJ
               INFO = J
               EXIT
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
!
!           Compute elements J+1:N of column J.
!
            IF ( J < N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ), &
                           LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
         END DO
      END IF

      END SUBROUTINE DPOTF2


!#13
!*****************************************************************************
      SUBROUTINE DSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
!*****************************************************************************
!
!  Purpose
!  =======
!
!  DSYR   performs the symmetric rank 1 operation
!
!     A := alpha*x*x' + A,
!
!  where alpha is a real scalar, x is an n element vector and A is an
!  n by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================

!      USE BLAS_MODULE, ONLY : XERBLA

      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER (KIND=4), INTENT(IN) :: INCX, LDA, N
      REAL    (KIND=8), INTENT(IN) :: ALPHA
      CHARACTER(LEN=1), INTENT(IN) :: UPLO

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(INOUT) :: A( LDA, * )
      REAL    (KIND=8), INTENT(IN)    :: X( * )

!     .. Local Scalars ..
      REAL    (KIND=8) :: TEMP
      INTEGER (KIND=4) :: I, INFO, IX, J, JX, KX

!     .. External Functions ..
!     LOGICAL :: LSAME
!     EXTERNAL LSAME

!     .. External Subroutines ..
!     EXTERNAL XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC MAX
!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF      ( .NOT.LSAME( UPLO, 'U' ).AND. &
                .NOT.LSAME( UPLO, 'L' )      ) THEN
         INFO = 1
      ELSE IF ( N < 0 ) THEN
         INFO = 2
      ELSE IF ( INCX == 0 ) THEN
         INFO = 5
      ELSE IF ( LDA < MAX( 1, N ) ) THEN
         INFO = 7
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DSYR  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF ( ( N == 0 ).OR.( ALPHA == ZERO ) ) RETURN
!
!     Set the start point in X if the increment is not unity.
!
      IF ( INCX <= 0 ) THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF ( INCX /= 1 ) THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      IF ( LSAME( UPLO, 'U' ) ) THEN
!
!        Form  A  when A is stored in upper triangle.
!
         IF ( INCX == 1 ) THEN
            DO J = 1, N
               IF ( X( J ) /= ZERO ) THEN
                  TEMP = ALPHA*X( J )
                  DO I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP
                  END DO
               END IF
            END DO
         ELSE
            JX = KX
            DO J = 1, N
               IF ( X( JX ) /= ZERO ) THEN
                  TEMP = ALPHA*X( JX )
                  IX   = KX
                  DO I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
                  END DO
               END IF
               JX = JX + INCX
            END DO
         END IF
      ELSE
!
!        Form  A  when A is stored in lower triangle.
!
         IF ( INCX == 1 ) THEN
            DO J = 1, N
               IF ( X( J ) /= ZERO ) THEN
                  TEMP = ALPHA*X( J )
                  DO I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP
                  END DO
               END IF
            END DO
         ELSE
            JX = KX
            DO J = 1, N
               IF ( X( JX ) /= ZERO ) THEN
                  TEMP = ALPHA*X( JX )
                  IX   = JX
                  DO I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP
                     IX        = IX        + INCX
                  END DO
               END IF
               JX = JX + INCX
            END DO
         END IF
      END IF

      END SUBROUTINE DSYR


!#14
!*****************************************************************************
      SUBROUTINE DSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA, &
                         BETA, C, LDC )
!*****************************************************************************
!
!  Purpose
!  =======
!
!  DSYRK  performs one of the symmetric rank k operations
!
!     C := alpha*A*A' + beta*C,
!
!  or
!
!     C := alpha*A'*A + beta*C,
!
!  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
!  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
!  in the second case.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of the  array  C  is to be  referenced  as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry,  TRANS  specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
!
!              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
!
!              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N specifies the order of the matrix C.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!           of  columns   of  the   matrix   A,   and  on   entry   with
!           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
!           of rows of the matrix  A.  K must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by n  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!           be at least  max( 1, k ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!           upper triangular part of the array C must contain the upper
!           triangular part  of the  symmetric matrix  and the strictly
!           lower triangular part of C is not referenced.  On exit, the
!           upper triangular part of the array  C is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!           lower triangular part of the array C must contain the lower
!           triangular part  of the  symmetric matrix  and the strictly
!           upper triangular part of C is not referenced.  On exit, the
!           lower triangular part of the array  C is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!=========================================================================

!      USE BLAS_MODULE, ONLY : XERBLA
      
      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN) :: UPLO, TRANS
      INTEGER (KIND=4), INTENT(IN) :: N, K, LDA, LDC
      REAL    (KIND=8), INTENT(IN) :: ALPHA, BETA

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(IN)    :: A( LDA, * )
      REAL    (KIND=8), INTENT(INOUT) :: C( LDC, * )

!     .. External Functions ..
!     LOGICAL :: LSAME
!     EXTERNAL LSAME

!     .. External Subroutines ..
!     EXTERNAL XERBLA

!     .. Intrinsic Functions ..
      INTRINSIC MAX

!     .. Local Scalars ..
      INTEGER (KIND=4) :: I, INFO, J, L, NROWA
      REAL    (KIND=8) :: TEMP
      LOGICAL :: UPPER

!     .. Executable Statements ..

!     Test the input parameters.

      IF ( LSAME( TRANS, 'N' ) ) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
!
      INFO = 0
      IF (      ( .NOT.UPPER               ).AND. &
                ( .NOT.LSAME( UPLO , 'L' ) )      ) THEN
         INFO = 1
      ELSE IF ( ( .NOT.LSAME( TRANS, 'N' ) ).AND. &
                ( .NOT.LSAME( TRANS, 'T' ) ).AND. &
                ( .NOT.LSAME( TRANS, 'C' ) )      ) THEN
         INFO = 2
      ELSE IF ( N   < 0               ) THEN
         INFO = 3
      ELSE IF ( K   < 0               ) THEN
         INFO = 4
      ELSE IF ( LDA < MAX( 1, NROWA ) ) THEN
         INFO = 7
      ELSE IF ( LDC < MAX( 1, N     ) ) THEN
         INFO = 10
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DSYRK ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF ( ( N == 0 ).OR. &
          ( ( ( ALPHA == ZERO ).OR.( K == 0 ) ).AND.( BETA == ONE ) ) ) THEN
         RETURN
      END IF
!
!     And when  alpha.eq.zero.
!
      IF ( ALPHA == ZERO ) THEN
         IF ( UPPER ) THEN
            IF ( BETA == ZERO ) THEN
               DO J = 1, N
                  DO I = 1, J
                     C( I, J ) = ZERO
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  DO I = 1, J
                     C( I, J ) = BETA*C( I, J )
                  END DO
               END DO
            END IF
         ELSE
            IF ( BETA == ZERO ) THEN
               DO J = 1, N
                  DO I = J, N
                     C( I, J ) = ZERO
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  DO I = J, N
                     C( I, J ) = BETA*C( I, J )
                  END DO
               END DO
            END IF
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF ( LSAME( TRANS, 'N' ) ) THEN
!
!        Form  C := alpha*A*A' + beta*C.
!
         IF ( UPPER ) THEN
            DO J = 1, N
               IF ( BETA == ZERO ) THEN
                  DO I = 1, J
                     C( I, J ) = ZERO
                  END DO
               ELSE IF ( BETA /= ONE ) THEN
                  DO I = 1, J
                     C( I, J ) = BETA*C( I, J )
                  END DO
               END IF
               DO L = 1, K
                  IF ( A( J, L ) /= ZERO ) THEN
                     TEMP = ALPHA*A( J, L )
                     DO I = 1, J
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
                     END DO
                  END IF
               END DO
            END DO
         ELSE
            DO J = 1, N
               IF ( BETA == ZERO ) THEN
                  DO I = J, N
                     C( I, J ) = ZERO
                  END DO
               ELSE IF ( BETA /= ONE ) THEN
                  DO I = J, N
                     C( I, J ) = BETA*C( I, J )
                  END DO
               END IF
               DO L = 1, K
                  IF ( A( J, L ) /= ZERO ) THEN
                     TEMP      = ALPHA*A( J, L )
                     DO I = J, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
                     END DO
                  END IF
               END DO
            END DO
         END IF
      ELSE
!
!        Form  C := alpha*A'*A + beta*C.
!
         IF ( UPPER ) THEN
            DO J = 1, N
               DO I = 1, J
                  TEMP = ZERO
                  DO L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
                  END DO
                  IF ( BETA == ZERO ) THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
               END DO
            END DO
         ELSE
            DO J = 1, N
               DO I = J, N
                  TEMP = ZERO
                  DO L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
                  END DO
                  IF ( BETA == ZERO ) THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
               END DO
            END DO
         END IF
      END IF

      END SUBROUTINE DSYRK


!#15
!*****************************************************************************
      SUBROUTINE DTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
!*****************************************************************************
!
!  Purpose
!  =======
!
!  DTBSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0 .le. K.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO J = 1, N
!                    M = K + 1 - J
!                    DO I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!                    END DO
!                 END DO
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO J = 1, N
!                    M = 1 - J
!                    DO I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!                    END DO
!                 END DO
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )!abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!  =====================================================================

!      USE BLAS_MODULE, ONLY : XERBLA
      
      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER (KIND=4), INTENT(IN) :: INCX, K, LDA, N
      CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(IN)    :: A( LDA, * )
      REAL    (KIND=8), INTENT(INOUT) :: X( * )

!     .. Local Scalars ..
      REAL    (KIND=8) :: TEMP
      INTEGER (KIND=4) :: I, INFO, IX, J, JX, KPLUS1, KX, L
      LOGICAL :: NOUNIT

!     .. Intrinsic Functions ..
      INTRINSIC MAX, MIN
!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF      ( .NOT.LSAME( UPLO , 'U' ).AND. &
                .NOT.LSAME( UPLO , 'L' )      ) THEN
         INFO = 1
      ELSE IF ( .NOT.LSAME( TRANS, 'N' ).AND. &
                .NOT.LSAME( TRANS, 'T' ).AND. &
                .NOT.LSAME( TRANS, 'C' )      ) THEN
         INFO = 2
      ELSE IF ( .NOT.LSAME( DIAG , 'U' ).AND. &
                .NOT.LSAME( DIAG , 'N' )      ) THEN
         INFO = 3
      ELSE IF ( N < 0 ) THEN
         INFO = 4
      ELSE IF ( K < 0 ) THEN
         INFO = 5
      ELSE IF ( LDA < ( K + 1 ) ) THEN
         INFO = 7
      ELSE IF ( INCX == 0 ) THEN
         INFO = 9
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DTBSV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF ( N == 0 ) RETURN
!
      NOUNIT = LSAME( DIAG, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF      ( INCX <= 0 ) THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF ( INCX /= 1 ) THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
      IF ( LSAME( TRANS, 'N' ) ) THEN
!
!        Form  x := inv( A )*x.
!
         IF ( LSAME( UPLO, 'U' ) ) THEN
            KPLUS1 = K + 1
            IF ( INCX == 1 ) THEN
               DO J = N, 1, -1
                  IF ( X( J ) /= ZERO ) THEN
                     L = KPLUS1 - J
                     IF ( NOUNIT ) X( J ) = X( J )/A( KPLUS1, J )
                     TEMP = X( J )
                     DO I = J - 1, MAX( 1, J - K ), -1
                        X( I ) = X( I ) - TEMP*A( L + I, J )
                     END DO
                  END IF
               END DO
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO J = N, 1, -1
                  KX = KX - INCX
                  IF ( X( JX ) /= ZERO ) THEN
                     IX = KX
                     L  = KPLUS1 - J
                     IF ( NOUNIT ) X( JX ) = X( JX )/A( KPLUS1, J )
                     TEMP = X( JX )
                     DO I = J - 1, MAX( 1, J - K ), -1
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      - INCX
                     END DO
                  END IF
                  JX = JX - INCX
               END DO
            END IF
         ELSE
            IF ( INCX == 1 ) THEN
               DO J = 1, N
                  IF ( X( J ) /= ZERO ) THEN
                     L = 1 - J
                     IF ( NOUNIT ) X( J ) = X( J )/A( 1, J )
                     TEMP = X( J )
                     DO I = J + 1, MIN( N, J + K )
                        X( I ) = X( I ) - TEMP*A( L + I, J )
                     END DO
                  END IF
               END DO
            ELSE
               JX = KX
               DO J = 1, N
                  KX = KX + INCX
                  IF ( X( JX ) /= ZERO ) THEN
                     IX = KX
                     L  = 1  - J
                     IF ( NOUNIT ) X( JX ) = X( JX )/A( 1, J )
                     TEMP = X( JX )
                     DO I = J + 1, MIN( N, J + K )
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      + INCX
                     END DO
                  END IF
                  JX = JX + INCX
               END DO
            END IF
         END IF
      ELSE
!
!        Form  x := inv( A')*x.
!
         IF ( LSAME( UPLO, 'U' ) ) THEN
            KPLUS1 = K + 1
            IF ( INCX == 1 ) THEN
               DO J = 1, N
                  TEMP = X( J )
                  L    = KPLUS1 - J
                  DO I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( I )
                  END DO
                  IF ( NOUNIT ) TEMP = TEMP/A( KPLUS1, J )
                  X( J ) = TEMP
               END DO
            ELSE
               JX = KX
               DO J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  L    = KPLUS1  - J
                  DO I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   + INCX
                  END DO
                  IF ( NOUNIT ) TEMP = TEMP/A( KPLUS1, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  IF ( J > K ) KX = KX + INCX
               END DO
            END IF
         ELSE
            IF ( INCX == 1 ) THEN
               DO J = N, 1, -1
                  TEMP = X( J )
                  L    = 1      - J
                  DO I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( I )
                  END DO
                  IF ( NOUNIT ) TEMP = TEMP/A( 1, J )
                  X( J ) = TEMP
               END DO
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  L    = 1       - J
                  DO I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   - INCX
                  END DO
                  IF ( NOUNIT ) TEMP = TEMP/A( 1, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  IF ( ( N - J ) >= K ) KX = KX - INCX
               END DO
            END IF
         END IF
      END IF

      END SUBROUTINE DTBSV


!#16
!*****************************************************************************
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
                         B, LDB )
!*****************************************************************************
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!  =====================================================================

!      USE BLAS_MODULE, ONLY : XERBLA

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN) :: SIDE, UPLO, TRANSA, DIAG
      INTEGER (KIND=4), INTENT(IN) :: M, N, LDA, LDB
      REAL    (KIND=8), INTENT(IN) :: ALPHA

!     .. Array Arguments ..
      REAL    (KIND=8), INTENT(IN)    :: A( LDA, * )
      REAL    (KIND=8), INTENT(INOUT) :: B( LDB, * )

!     .. Intrinsic Functions ..
      INTRINSIC MAX

!     .. Local Scalars ..
      INTEGER (KIND=4) :: I, INFO, J, K, NROWA
      REAL    (KIND=8) :: TEMP
      LOGICAL :: LSIDE, NOUNIT, UPPER

!
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF ( LSIDE ) THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )

      INFO   = 0
      IF (      ( .NOT.LSIDE                ).AND. &
                ( .NOT.LSAME( SIDE  , 'R' ) )      ) THEN
         INFO = 1
      ELSE IF ( ( .NOT.UPPER                ).AND. &
                ( .NOT.LSAME( UPLO  , 'L' ) )      ) THEN
         INFO = 2
      ELSE IF ( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
                ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
                ( .NOT.LSAME( TRANSA, 'C' ) )      ) THEN
         INFO = 3
      ELSE IF ( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
                ( .NOT.LSAME( DIAG  , 'N' ) )      ) THEN
         INFO = 4
      ELSE IF ( M   < 0               ) THEN
         INFO = 5
      ELSE IF ( N   < 0               ) THEN
         INFO = 6
      ELSE IF ( LDA < MAX( 1, NROWA ) ) THEN
         INFO = 9
      ELSE IF ( LDB < MAX( 1, M     ) ) THEN
         INFO = 11
      END IF
      IF ( INFO /= 0 ) THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF ( N == 0 ) RETURN
!
!     And when  alpha.eq.zero.
!
      IF ( ALPHA == ZERO ) THEN
         DO J = 1, N
            DO I = 1, M
               B( I, J ) = ZERO
            END DO
         END DO
         RETURN
      END IF
!
!     Start the operations.
!
      IF ( LSIDE ) THEN
         IF ( LSAME( TRANSA, 'N' ) ) THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF ( UPPER ) THEN
               DO J = 1, N
                  IF ( ALPHA /= ONE ) THEN
                     DO I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
                     END DO
                  END IF
                  DO K = M, 1, -1
                     IF ( B( K, J ) /= ZERO ) THEN
                        IF ( NOUNIT ) &
                           B( K, J ) = B( K, J )/A( K, K )
                        DO I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                        END DO
                     END IF
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  IF ( ALPHA /= ONE ) THEN
                     DO I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
                     END DO
                  END IF
                  DO K = 1, M
                     IF ( B( K, J ) /= ZERO ) THEN
                        IF ( NOUNIT ) &
                           B( K, J ) = B( K, J )/A( K, K )
                        DO I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
                        END DO
                     END IF
                  END DO
               END DO
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
            IF ( UPPER ) THEN
               DO J = 1, N
                  DO I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
                     END DO
                     IF ( NOUNIT ) TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
                  END DO
               END DO
            ELSE
               DO J = 1, N
                  DO I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
                     END DO
                     IF ( NOUNIT )  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
                  END DO
               END DO
            END IF
         END IF
      ELSE
         IF ( LSAME( TRANSA, 'N' ) ) THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF ( UPPER ) THEN
               DO J = 1, N
                  IF ( ALPHA /= ONE ) THEN
                     DO I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
                     END DO
                  END IF
                  DO K = 1, J - 1
                     IF ( A( K, J ) /= ZERO ) THEN
                        DO I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                        END DO
                     END IF
                  END DO
                  IF ( NOUNIT ) THEN
                     TEMP = ONE/A( J, J )
                     DO I = 1, M
                        B( I, J ) = TEMP*B( I, J )
                     END DO
                  END IF
               END DO
            ELSE
               DO J = N, 1, -1
                  IF ( ALPHA /= ONE ) THEN
                     DO I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
                     END DO
                  END IF
                  DO K = J + 1, N
                     IF ( A( K, J ) /= ZERO ) THEN
                        DO I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
                        END DO
                     END IF
                  END DO
                  IF ( NOUNIT ) THEN
                     TEMP = ONE/A( J, J )
                     DO I = 1, M
                       B( I, J ) = TEMP*B( I, J )
                     END DO
                  END IF
               END DO
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            IF ( UPPER ) THEN
               DO K = N, 1, -1
                  IF ( NOUNIT ) THEN
                     TEMP = ONE/A( K, K )
                     DO I = 1, M
                        B( I, K ) = TEMP*B( I, K )
                     END DO
                  END IF
                  DO J = 1, K - 1
                     IF ( A( J, K ) /= ZERO ) THEN
                        TEMP = A( J, K )
                        DO I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
                        END DO
                     END IF
                  END DO
                  IF ( ALPHA /= ONE ) THEN
                     DO I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
                     END DO
                  END IF
               END DO
            ELSE
               DO K = 1, N
                  IF ( NOUNIT ) THEN
                     TEMP = ONE/A( K, K )
                     DO I = 1, M
                        B( I, K ) = TEMP*B( I, K )
                     END DO
                  END IF
                  DO J = K + 1, N
                     IF ( A( J, K ) /= ZERO ) THEN
                        TEMP = A( J, K )
                        DO I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
                        END DO
                     END IF
                  END DO
                  IF ( ALPHA /= ONE ) THEN
                     DO I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
                     END DO
                  END IF
               END DO
            END IF
         END IF
      END IF

      END SUBROUTINE DTRSM


!#17
!*****************************************************************************
!     INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      INTEGER FUNCTION ILAENV( ISPEC, NAME, N1, N2, N4 )
!*****************************************************************************
!
!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 20, 1992
!
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF ( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=*), INTENT(IN) :: NAME !, OPTS
      INTEGER (KIND=4), INTENT(IN) :: ISPEC, N1, N2, N4
!     INTEGER (KIND=4), INTENT(IN) :: ISPEC, N1, N2, N3, N4

!     .. Local Scalars ..
      LOGICAL :: CNAME, SNAME
      CHARACTER(LEN=1) :: C1
      CHARACTER(LEN=2) :: C2, C4
      CHARACTER(LEN=3) :: C3
      CHARACTER(LEN=6) :: SUBNAM
      INTEGER (KIND=4) :: I, IC, IZ, NB, NBMIN, NX

!     .. Intrinsic Functions ..
      INTRINSIC CHAR, ICHAR, INT, MIN, REAL
!
!     .. Executable Statements ..
!
      SELECT CASE (ISPEC)
      
      CASE (1, 2, 3)
!
!     Convert NAME to upper case if the first character is lower case.
!
         ILAENV = 1
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1:1 ) )
         IZ = ICHAR( 'Z' )
         IF ( IZ == 90 .OR. IZ == 122 ) THEN
!
!           ASCII character set
!
            IF ( IC >= 97 .AND. IC <= 122 ) THEN
               SUBNAM( 1:1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I:I ) )
                  IF ( IC >= 97 .AND. IC <= 122 ) SUBNAM( I:I ) = CHAR( IC-32 )
               END DO
            END IF

         ELSE IF ( IZ == 233 .OR. IZ == 169 ) THEN
!
!           EBCDIC character set
!
            IF ( ( IC >= 129 .AND. IC <= 137 ) .OR. &
                 ( IC >= 145 .AND. IC <= 153 ) .OR. &
                 ( IC >= 162 .AND. IC <= 169 ) ) THEN
               SUBNAM( 1:1 ) = CHAR( IC+64 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I:I ) )
                  IF (      ( IC >= 129 .AND. IC <= 137 ) &
                       .OR. ( IC >= 145 .AND. IC <= 153 ) &
                       .OR. ( IC >= 162 .AND. IC <= 169 ) ) THEN
                     SUBNAM( I:I ) = CHAR( IC+64 )
                  END IF
               END DO
            END IF
 
         ELSE IF ( IZ == 218 .OR. IZ == 250 ) THEN
!
!           Prime machines:  ASCII+128
!
            IF ( IC >= 225 .AND. IC <= 250 ) THEN
               SUBNAM( 1:1 ) = CHAR( IC-32 )
               DO I = 2, 6
                  IC = ICHAR( SUBNAM( I:I ) )
                  IF ( IC >= 225 .AND. IC <= 250 ) SUBNAM( I:I ) = CHAR( IC-32 )
               END DO
            END IF
         END IF

         C1 = SUBNAM( 1:1 )
         SNAME = (C1 == 'S' .OR. C1 == 'D')
         CNAME = (C1 == 'C' .OR. C1 == 'Z')
         IF ( .NOT.( CNAME .OR. SNAME ) ) RETURN
         C2 = SUBNAM( 2:3 )
         C3 = SUBNAM( 4:6 )
         C4 = C3( 2:3 )
 
         SELECT CASE ( ISPEC )

         CASE (1)
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
            NB = 1
 
            IF ( C2 == 'GE' ) THEN
               IF ( C3 == 'TRF' ) THEN
                  IF ( SNAME ) THEN
                     NB = 64
                  ELSE
                     NB = 64
                  END IF
               ELSE IF ( C3 == 'QRF' .OR. &
                         C3 == 'RQF' .OR. &
                         C3 == 'LQF' .OR. &
                         C3 == 'QLF' ) THEN
                  IF ( SNAME ) THEN
                     NB = 32
                  ELSE
                     NB = 32
                  END IF
               ELSE IF ( C3 == 'HRD' ) THEN
                  IF ( SNAME ) THEN
                     NB = 32
                  ELSE
                     NB = 32
                  END IF
               ELSE IF ( C3 == 'BRD' ) THEN
                  IF ( SNAME ) THEN
                     NB = 32
                  ELSE
                     NB = 32
                  END IF
               ELSE IF ( C3 == 'TRI' ) THEN
                  IF ( SNAME ) THEN
                     NB = 64
                  ELSE
                     NB = 64
                  END IF
               END IF
            ELSE IF ( C2 == 'PO' ) THEN
               IF ( C3 == 'TRF' ) THEN
                  IF ( SNAME ) THEN
                     NB = 64
                  ELSE
                     NB = 64
                  END IF
               END IF
            ELSE IF ( C2 == 'SY' ) THEN
               IF ( C3 == 'TRF' ) THEN
                  IF ( SNAME ) THEN
                     NB = 64
                  ELSE
                     NB = 64
                  END IF
               ELSE IF ( SNAME .AND. C3 == 'TRD' ) THEN
                  NB = 1
               ELSE IF ( SNAME .AND. C3 == 'GST' ) THEN
                  NB = 64
               END IF
            ELSE IF ( CNAME .AND. C2 == 'HE' ) THEN
               IF ( C3 == 'TRF' ) THEN
                  NB = 64
               ELSE IF ( C3 == 'TRD' ) THEN
                  NB = 1
               ELSE IF ( C3 == 'GST' ) THEN
                  NB = 64
               END IF
            ELSE IF ( SNAME .AND. C2 == 'OR' ) THEN
               IF ( C3( 1:1 ).EQ.'G' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NB = 32
                  END IF
               ELSE IF ( C3( 1:1 ) == 'M' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NB = 32
                  END IF
               END IF
            ELSE IF ( CNAME .AND. C2 == 'UN' ) THEN
               IF ( C3( 1:1 ) == 'G' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NB = 32
                  END IF
               ELSE IF ( C3( 1:1 ) == 'M' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NB = 32
                  END IF
               END IF
            ELSE IF ( C2 == 'GB' ) THEN
               IF ( C3 == 'TRF' ) THEN
                  IF ( SNAME ) THEN
                     IF ( N4 <= 64 ) THEN
                        NB = 1
                     ELSE
                        NB = 32
                     END IF
                  ELSE
                     IF ( N4 <= 64 ) THEN
                        NB = 1
                     ELSE
                        NB = 32
                     END IF
                  END IF
               END IF
            ELSE IF ( C2 == 'PB' ) THEN
               IF ( C3 == 'TRF' ) THEN
                  IF ( SNAME ) THEN
                     IF ( N2 <= 64 ) THEN
                        NB = 1
                     ELSE
                        NB = 32
                     END IF
                  ELSE
                     IF ( N2 <= 64 ) THEN
                        NB = 1
                     ELSE
                        NB = 32
                     END IF
                  END IF
               END IF
            ELSE IF ( C2 == 'TR' ) THEN
               IF ( C3 == 'TRI' ) THEN
                  IF ( SNAME ) THEN
                     NB = 64
                  ELSE
                     NB = 64
                  END IF
               END IF
            ELSE IF ( C2 == 'LA' ) THEN
               IF ( C3 == 'UUM' ) THEN
                  IF ( SNAME ) THEN
                     NB = 64
                  ELSE
                     NB = 64
                  END IF
               END IF
            ELSE IF ( SNAME .AND. C2 == 'ST' ) THEN
               IF ( C3 == 'EBZ' ) THEN
                  NB = 1
               END IF
            END IF
            ILAENV = NB

         CASE (2)
!
!     ISPEC = 2:  minimum block size
!
            NBMIN = 2
            IF ( C2 == 'GE' ) THEN
               IF ( C3 == 'QRF' .OR. &
                    C3 == 'RQF' .OR. &
                    C3 == 'LQF' .OR. &
                    C3 == 'QLF' ) THEN
                  IF ( SNAME ) THEN
                     NBMIN = 2
                  ELSE
                      NBMIN = 2
                  END IF
               ELSE IF ( C3 == 'HRD' ) THEN
                  IF ( SNAME ) THEN
                     NBMIN = 2
                  ELSE
                     NBMIN = 2
                  END IF
               ELSE IF ( C3 == 'BRD' ) THEN
                  IF ( SNAME ) THEN
                     NBMIN = 2
                  ELSE
                     NBMIN = 2
                  END IF
               ELSE IF ( C3 == 'TRI' ) THEN
                  IF ( SNAME ) THEN
                     NBMIN = 2
                  ELSE
                     NBMIN = 2
                  END IF
               END IF
            ELSE IF ( C2 == 'SY' ) THEN
               IF ( C3 == 'TRF' ) THEN
                  IF ( SNAME ) THEN
                     NBMIN = 2
                  ELSE
                     NBMIN = 2
                  END IF
               ELSE IF ( SNAME .AND. C3 == 'TRD' ) THEN
                  NBMIN = 2
               END IF
            ELSE IF ( CNAME .AND. C2 == 'HE' ) THEN
               IF ( C3 == 'TRD' ) THEN
                  NBMIN = 2
               END IF
            ELSE IF ( SNAME .AND. C2 == 'OR' ) THEN
               IF ( C3( 1:1 ) == 'G' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NBMIN = 2
                  END IF
               ELSE IF ( C3( 1:1 ) == 'M' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NBMIN = 2
                  END IF
               END IF
            ELSE IF ( CNAME .AND. C2 == 'UN' ) THEN
               IF ( C3( 1:1 ) == 'G' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NBMIN = 2
                  END IF
               ELSE IF ( C3( 1:1 ) == 'M' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NBMIN = 2
                  END IF
               END IF
            END IF
            ILAENV = NBMIN

         CASE (3)
!
!     ISPEC = 3:  crossover point
!
            NX = 0
            IF ( C2 == 'GE' ) THEN
               IF ( C3 == 'QRF' .OR. &
                    C3 == 'RQF' .OR. &
                    C3 == 'LQF' .OR. &
                    C3 == 'QLF' ) THEN
                  IF ( SNAME ) THEN
                     NX = 128
                  ELSE
                     NX = 128
                  END IF
               ELSE IF ( C3 == 'HRD' ) THEN
                  IF ( SNAME ) THEN
                     NX = 128
                  ELSE
                     NX = 128
                  END IF
               ELSE IF ( C3 == 'BRD' ) THEN
                  IF ( SNAME ) THEN
                     NX = 128
                  ELSE
                     NX = 128
                  END IF
               END IF
            ELSE IF ( C2 == 'SY' ) THEN
               IF ( SNAME .AND. C3 == 'TRD' ) THEN
                  NX = 1
               END IF
            ELSE IF ( CNAME .AND. C2 == 'HE' ) THEN
               IF ( C3 == 'TRD' ) THEN
                  NX = 1
               END IF
            ELSE IF ( SNAME .AND. C2 == 'OR' ) THEN
               IF ( C3( 1:1 ) == 'G' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NX = 128
                  END IF
               END IF
            ELSE IF ( CNAME .AND. C2 == 'UN' ) THEN
               IF ( C3( 1:1 ) == 'G' ) THEN
                  IF ( C4 == 'QR' .OR. &
                       C4 == 'RQ' .OR. &
                       C4 == 'LQ' .OR. &
                       C4 == 'QL' .OR. &
                       C4 == 'HR' .OR. &
                       C4 == 'TR' .OR. &
                       C4 == 'BR' ) THEN
                     NX = 128
                  END IF
               END IF
            END IF
            ILAENV = NX

         END SELECT

      CASE (4)
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
         ILAENV = 6

      CASE (5)
!
!     ISPEC = 5:  minimum column dimension (not used)
!
         ILAENV = 2

      CASE (6)
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
         ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )

      CASE (7)
!
!     ISPEC = 7:  number of processors (not used)
!
         ILAENV = 1

      CASE (8)
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
         ILAENV = 50

      CASE DEFAULT
!
!     Invalid value for ISPEC
!
         ILAENV = -1
  
      END SELECT

      END FUNCTION ILAENV


!#18
!*****************************************************************************
      LOGICAL FUNCTION LSAME( CA, CB )
!*****************************************************************************
!
!  -- LAPACK auxiliary routine (version 1.1) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992 
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
!-----------------------------------------------------------------------------

      IMPLICIT NONE

!     .. Scalar Arguments ..
      CHARACTER(LEN=1), INTENT(IN) :: CA, CB

!     .. Intrinsic Functions ..
      INTRINSIC ICHAR

!     .. Local Scalars ..
      INTEGER (KIND=4) :: INTA, INTB, ZCODE
!
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = (CA == CB)
      IF ( LSAME ) RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF ( ZCODE == 90 .OR. ZCODE == 122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF ( INTA >= 97 .AND. INTA <= 122 ) INTA = INTA - 32
         IF ( INTB >= 97 .AND. INTB <= 122 ) INTB = INTB - 32
!
      ELSE IF ( ZCODE == 233 .OR. ZCODE == 169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF ( INTA >= 129 .AND. INTA <= 137 .OR. &
              INTA >= 145 .AND. INTA <= 153 .OR. &
              INTA >= 162 .AND. INTA <= 169 ) INTA = INTA + 64
         IF ( INTB >= 129 .AND. INTB <= 137 .OR. &
              INTB >= 145 .AND. INTB <= 153 .OR. &
              INTB >= 162 .AND. INTB <= 169 ) INTB = INTB + 64
!
      ELSE IF ( ZCODE == 218 .OR. ZCODE == 250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF ( INTA >= 225 .AND. INTA <= 250 ) INTA = INTA - 32
         IF ( INTB >= 225 .AND. INTB <= 250 ) INTB = INTB - 32
      END IF
      LSAME = (INTA == INTB)
!
!     End of LSAME
!
      END FUNCTION LSAME


      END MODULE BLAS_LAPACK_MODULE
