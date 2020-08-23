      REAL*8 FUNCTION RESIDFUN(M,N,B,BETA)

C#### Function: RESIDFUN
C###  Type: REAL*8
C###  Description:
C###    Calculates the intrisic residual.
CC JMB 10-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER M, N
      REAL*8 B(*), BETA(*)
!     Local Variables
      REAL*8 DELTA, BETA2, ZERO
      PARAMETER (ZERO = 0.0d0)
!     Functions
      REAL*8 DDOT

      BETA2 = DDOT(M, B, 1, B, 1) - DDOT(MIN(M,N), BETA, 1, BETA, 1)
      DELTA = ZERO
      IF( (M.GT.N).AND.(BETA2.GT.ZERO) ) THEN
        DELTA = BETA2
      ENDIF

      RESIDFUN = DELTA

      RETURN
      END

C---------------------------------------------------------------------
