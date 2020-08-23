      REAL*8 FUNCTION PL2(I,K,XI)

C#### Function: PL2
C###  Type: REAL*8
C###  Description:
C###    PL2 evaluates 1D quadratic Lagrange basis function at XI.

      IMPLICIT NONE
!     Parameter List
      INTEGER I,K
      REAL*8 XI

      GO TO (10,20,30),K
 10     GO TO (11,12,13),I
 11       PL2=1.0d0-3.0d0*XI+2.0d0*XI*XI
          RETURN
 12       PL2=4.0d0*XI*(1.0d0-XI)
          RETURN
 13       PL2=XI*(XI+XI-1.0d0)
          RETURN
 20     GO TO (21,22,23),I
 21       PL2=4.0d0*XI-3.0d0
          RETURN
 22       PL2=4.0d0-8.0d0*XI
          RETURN
 23       PL2=4.0d0*XI-1.0d0
          RETURN
 30     GO TO (31,32,33),I
 31       PL2=4.0d0
          RETURN
 32       PL2=-8.0d0
          RETURN
 33       PL2=4.0d0
          RETURN
      END


