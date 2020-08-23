      REAL*8 FUNCTION QUASIFUN(LAMBDA,M,N,BETA,S,WORK)

C#### Function: QUASIFUN
C###  Type: REAL*8
C###  Description:
C###   Auxiliary routine for QUASIOPT. Computes the quasi-optimality
C###  function.
CC JMB 13-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER M, N
      REAL*8 BETA(*), LAMBDA, S(*), WORK(*)
!     Local Variables
      INTEGER i, MN
      REAL*8 F, ONE, XI
      PARAMETER (ONE = 1.0d0)
!     Functions
      REAL*8 DNRM2

      ! Initialisation
      MN = MIN(M,N)

      DO i = 1,MN
        F = S(i)**2/(S(i)**2 + LAMBDA**2)
        XI = BETA(i)/S(i)
        WORK(i) = (ONE - F)*F*XI
      ENDDO

      QUASIFUN = DNRM2(MN, WORK(1), 1)

      RETURN
      END

C---------------------------------------------------------------------
