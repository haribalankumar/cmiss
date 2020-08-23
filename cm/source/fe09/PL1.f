      REAL*8 FUNCTION PL1(I,K,XI)

C#### Function: PL1
C###  Type: REAL*8
C###  Description:
C###    PL1 evaluates 1D linear Lagrange basis function at XI.

      IMPLICIT NONE
!     Parameter List
      INTEGER I,K
      REAL*8 XI

      GO TO (10,20,30),K
 10     GO TO (11,12),I
 11       PL1=1.0d0-XI
          RETURN
 12       PL1=XI
          RETURN
 20     GO TO (21,22),I
 21       PL1=-1.0d0
          RETURN
 22       PL1=1.0d0
          RETURN
 30     PL1=0.0d0
        RETURN
      END


