      REAL*8 FUNCTION BETA(A,B,C)

C#### Function:BETA
C###  Description:
C###    BETA is gradient of characteristic through Q.

      IMPLICIT NONE
!     Parameter List
      REAL*8 A,B,C

      BETA=(B-DSQRT(B*B-4.0D0*A*C))/2.0D0*A

      RETURN
      END


