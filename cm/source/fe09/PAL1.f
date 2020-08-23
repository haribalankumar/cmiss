      REAL*8 FUNCTION PAL1(K,XI)

C#### Function: PAL1
C###  Type: REAL*8
C###    PAL1 evaluates 1D linear Legendre auxiliary basis function
C###    at Xi.

      IMPLICIT NONE
!     Parameter List
      INTEGER K
      REAL*8 XI

      GO TO (10,20,30),K
 10     PAL1=XI
        RETURN
 20     PAL1=1.0d0
        RETURN
 30     PAL1=0.0d0
        RETURN
      END


