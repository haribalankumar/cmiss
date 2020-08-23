      REAL*8 FUNCTION ALFA(A,B,C)

C#### Function: ALFA
C###  Description:
C###    ALFA is gradient of characteristic through P.

      IMPLICIT NONE
!     Parameter List
      REAL*8 A,B,C

      ALFA=(B+DSQRT(B*B-4.0D0*A*C))/2.0D0*A

      RETURN
      END


