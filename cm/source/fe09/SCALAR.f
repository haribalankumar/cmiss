      REAL*8 FUNCTION SCALAR(N,A,B)

C#### Function: SCALAR
C###  Type: REAL*8
C###  Description:
C###    SCALAR calculates scalar product of two vectors A,B of length N.

      IMPLICIT NONE
!     Parameter List
      INTEGER N
      REAL*8 A(N),B(N)
!     Local Variables
      INTEGER i

      SCALAR=0.0D0
      DO i=1,n
        SCALAR=SCALAR+A(i)*B(i)
      ENDDO

      RETURN
      END


