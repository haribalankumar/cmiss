      REAL*8 FUNCTION ANG3P(A,B,C)

C#### Function: ANG3P
C###  Description:
C###    Calculates the angle between 2 vectors defined by the 3 points
C###    A,B,C.  The vectors are BA and BC.

      IMPLICIT NONE
!     Parameter List
      REAL*8 A(3),B(3),C(3)
!     Local variables
      INTEGER j
      REAL*8 DOT_PROD,LU,LV,U(3),V(3),ZERO_TOL

      ZERO_TOL=1.d-6 !error tolerance
      LU=0.d0
      LV=0.d0
      DO j=1,3
        U(j)=A(j)-B(j)
        V(j)=C(j)-B(j)
        LU=LU+U(j)**2.d0
        LV=LV+V(j)**2.d0
      ENDDO !j
      LU=DSQRT(LU)
      LV=DSQRT(LV)
      DO j=1,3
        U(j)=U(j)/LU
        V(j)=V(j)/LV
      ENDDO !j
      ANG3P=DOT_PROD(U,V)
      IF(ANG3P.GT.-1.d0-ZERO_TOL.AND.ANG3P.LT.-1.d0+ZERO_TOL) THEN
        ANG3P=DACOS(-1.d0) !fixes incase small rounding errors
      ELSE IF(ANG3P.GT.1.d0-ZERO_TOL.AND.ANG3P.LT.1.d0+ZERO_TOL) THEN
        ANG3P=DACOS(1.d0)
      ELSE
        ANG3P=DACOS(ANG3P)
      ENDIF

      RETURN
      END


