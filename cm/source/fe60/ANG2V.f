      REAL*8 FUNCTION ANG2V(U,V)

C#### Function: ANG2V
C###  Description:
C###    Calculates the angle between 2 vectors.

      IMPLICIT NONE
!     Parameter List
      INCLUDE 'tol00.cmn'
      REAL*8 U(3),V(3)
!     Local variables
      INTEGER j
      REAL*8 DOT_PROD,LU,LV

      LU=0.d0
      LV=0.d0
      DO j=1,3
        LU=LU+U(j)**2.d0
        LV=LV+V(j)**2.d0
      ENDDO !j
      LU=DSQRT(LU)
      LV=DSQRT(LV)
      DO j=1,3
        U(j)=U(j)/LU
        V(j)=V(j)/LV
      ENDDO !j
      ANG2V=DOT_PROD(U,V)
      IF(ANG2V.GT.-1.d0-LOOSE_TOL.AND.ANG2V.LT.-1.d0+LOOSE_TOL) THEN
        ANG2V=DACOS(-1.d0) !fixes incase small rounding errors
      ELSE IF(ANG2V.GT.1.d0-LOOSE_TOL.AND.ANG2V.LT.1.d0+LOOSE_TOL) THEN
        ANG2V=DACOS(1.d0)
      ELSE
        ANG2V=DACOS(ANG2V)
      ENDIF

      RETURN
      END


