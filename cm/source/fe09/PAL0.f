      REAL*8 FUNCTION PAL0(K)

C#### Function: PAL0
C###  Type: REAL*8
C###  Description:
C###    PAL0 evaluates 1D constant auxiliary basis function.

      IMPLICIT NONE
!     Parameter List
      INTEGER K

      GO TO (10,20,30),K
 10     PAL0=1.0d0
        RETURN
 20     PAL0=0.0d0
        RETURN
 30     PAL0=0.0d0
        RETURN
      END


