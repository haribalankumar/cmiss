      INTEGER FUNCTION LEFT_MULT(i,k,T)

C#### Function: LEFT_MULT
C###  Type: INTEGER
C###  Description:
C###    LEFT_MULT returns the left multiplicity of the ith knot in T.

      IMPLICIT NONE
!     Parameter List
      INTEGER i,k
!     Local Variables
      INTEGER ip
      REAL*8 T(*)

      LEFT_MULT=0
      DO ip=1,k+1
        IF(T(i-ip+1).EQ.T(i))LEFT_MULT=LEFT_MULT+1
      ENDDO

      RETURN
      END


