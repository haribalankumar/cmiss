      INTEGER FUNCTION RIGHT_MULT(i,k,T)

C#### Function: RIGHT_MULT
C###  Type: INTEGER
C###  Description:
C###    RIGHT_MULT returns the right multiplicity of the ith knot in T.

      IMPLICIT NONE
!     Parameter List
      REAL*8 T(*)
!     Local Variables
      INTEGER i,ip,k

      RIGHT_MULT=0
      DO ip=1,k+1
        IF(T(i+ip-1).EQ.T(i)) RIGHT_MULT=RIGHT_MULT+1
      ENDDO

      RETURN
      END


