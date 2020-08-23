      REAL*8 FUNCTION ZEROFUN(LAMBDA,M,N,BETA,S,WORK)

C#### Function: ZEROFUN
C###  Type: REAL*8
C###  Description:
C###    Auxiliary routine for ZEROCROSSING. Computes the -ve
C###  zero-crossing function.
CC JMB 13-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER M, N
      REAL*8 BETA(*), LAMBDA, S(*), WORK(*)
!     Local Variables
      INTEGER i, MN
      REAL*8 F
!     Functions
      REAL*8 DNRM2

      ! Initialisation
      MN = MIN(M,N)

      DO i = 1,MN
        F = S(i)**2/(S(i)**2 + LAMBDA**2)
        WORK(i + 1) = F*BETA(i)/S(i)
        WORK(MN + i + 1) = (1 - F)*BETA(i)
      ENDDO

      ZEROFUN = -LAMBDA**2*DNRM2(MN, WORK(2), 1)**2 +
     '  (DNRM2(MN, WORK(MN + 2), 1)**2 + WORK(1))

      RETURN
      END
