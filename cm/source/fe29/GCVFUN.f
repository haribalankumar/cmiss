      REAL*8 FUNCTION GCVFUN(LAMBDA,M,N,BETA,S,WORK)

C#### Function: GCVFUN
C###  Type: REAL*8
C###  Description:
C###    Auxiliary routine for GCV. Computes the general cross-
C###  validation (GCV) function.
CC JMB 13-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER M, N
      REAL*8 BETA(*), LAMBDA, S(*), WORK(*)
!     Local Variables
      INTEGER i, MN
      REAL*8 F, SUM, ZERO
      PARAMETER (ZERO = 0.0d0)
!     Functions
      REAL*8 DNRM2

      ! Initialisation
      MN = MIN(M,N)

      SUM = ZERO
      DO i = 1,MN
        F = LAMBDA**2/(S(i)**2 + LAMBDA**2)
        WORK(i + 1) = F*BETA(i)
        SUM = SUM + F
      ENDDO

      GCVFUN = (DNRM2(MN, WORK(2), 1)**2 + WORK(1))/(M - N + SUM)**2

      RETURN
      END

C---------------------------------------------------------------------
