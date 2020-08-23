      REAL*8 FUNCTION PAP4(K,XI)

C#### Function: PAP4
C###  Type: REAL*8
C###  Description:
C###    PAP4 evaluates 1D Pressure quartic (hat) auxiliary basis
C###    functionn at Xi.  Values and first derivatives are zero at
C###    Xi=0,1.

      IMPLICIT NONE
!     Parameter List
      INTEGER K
      REAL*8 XI

      GO TO (10,20,30),K
 10     PAP4=XI*XI*(XI-1.0D0)*(XI-1.0D0)
        RETURN
 20     PAP4=2.0D0*XI*(XI-1.0D0)*(2.0D0*XI-1.0D0)
        RETURN
 30     PAP4=12.0D0*XI*XI - 12.0D0*XI + 2.0D0
        RETURN
      END


