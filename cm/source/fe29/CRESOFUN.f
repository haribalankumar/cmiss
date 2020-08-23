      REAL*8 FUNCTION CRESOFUN(LAMBDA,M,N,BETA,S,WORK)

C#### Function: CRESOFUN
C###  Type: REAL*8
C###  Description:
C###   Auxiliary routine for CRESO. Computes the -ve CRESO function.
CC JMB 13-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER M, N
      REAL*8 BETA(*), LAMBDA, S(*), WORK(*)

!     Local Variables
      INTEGER i
      REAL*8 F, ONE, PI, SUM, XI, ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)

      SUM = ZERO
      DO i = 1,MIN(M,N)
        F = S(i)**2 + LAMBDA**2
        XI = S(i)*BETA(i)/F
        PI = ONE -(4.0d0*LAMBDA**2)/F
        SUM = SUM + XI**2*PI
      ENDDO

      CRESOFUN = -SUM

      RETURN
      END

C---------------------------------------------------------------------
