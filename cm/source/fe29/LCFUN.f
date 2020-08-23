      REAL*8 FUNCTION LCFUN(LAMBDA,M,N,BETA,S,WORK)

C#### Function: LCFUN
C###  Type: REAL*8
C###  Description:
C###    Auxililary routine for LCURVE. Computes the -ve of the
C###  curvature.
CC JMB 13-OCT-2000

      IMPLICIT NONE
!     Parameter List
      INTEGER M, N
      REAL*8 BETA(*), LAMBDA, S(*), WORK(*)
!     Local Variables
      INTEGER i, MN
      REAL*8 DDETA, DDLOGETA, DDLOGRHO, DDRHO, DETA, DLOGETA, DLOGRHO,
     '  DPHI, DPSI, DRHO, ETA, F, F1, F2, PHI, PSI, RHO, ZERO
      PARAMETER (ZERO = 0.0d0)
!     Functions
      REAL*8 DNRM2

      ! Initialisation
      MN = MIN(M,N)
      PHI = ZERO
      PSI = ZERO
      DPHI = ZERO
      DPSI = ZERO

      DO i = 1,MN
        F = S(i)**2/(S(i)**2 + LAMBDA**2)
        WORK(i + 1) = F*BETA(i)/S(i)
        WORK(MN + i + 1) = (1 - F)*BETA(i)
        F1 = -2.0d0*F*(1 - F)/LAMBDA
        F2 = -F1*(3.0d0 - 4.0d0*F)/LAMBDA
        PHI = PHI + F*F1*(BETA(i)/S(i))**2
        PSI = PSI + (1 - F)*F1*BETA(i)**2
        DPHI = DPHI + (F1**2 + F*F2)*(BETA(i)/S(i))**2
        DPSI = DPSI + (-F1**2 + (1 - F)*F2)*BETA(i)**2
      ENDDO

      ETA = DNRM2(MN, WORK(2), 1)
      RHO = DSQRT(DNRM2(MN, WORK(MN + 2), 1)**2 + WORK(1))

      ! Compute the first and second derivatives of ETA and RHO
      DETA = PHI/ETA
      DRHO = -PSI/RHO
      DDETA = DPHI/ETA - DETA*(DETA/ETA)
      DDRHO = -DPSI/RHO - DRHO*(DRHO/RHO)

      ! Convert derivatives to LOG(ETA) and LOG(RHO)
      DLOGETA = DETA/ETA
      DLOGRHO = DRHO/RHO
      DDLOGETA = DDETA/ETA - DLOGETA**2
      DDLOGRHO = DDRHO/RHO - DLOGRHO**2

      ! Compute -ve curvature
      LCFUN = -(DLOGRHO*DDLOGETA - DDLOGRHO*DLOGETA)/
     '  (DLOGRHO**2 + DLOGETA**2)**1.50d0

      RETURN
      END

C---------------------------------------------------------------------
