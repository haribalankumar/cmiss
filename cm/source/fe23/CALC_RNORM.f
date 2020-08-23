      REAL*8 FUNCTION CALC_RNORM(nBody,nHeart,NYNR,S_START,S_END,
     '  dPHIdTAU,d_TAU,PHIDIFF,YP)

C#### Function: CALC_RNORM
C###  Description:
C###    Computes the residual norm, which is the norm of
C###      dPhi/dTau.Tau_(k+1) - (Phi-Phi^) - dPhi/dTau.Tau_k

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER nBody,nHeart,NYNR(0:NY_R_M),S_START,S_END
      REAL*8 dPHIdTAU(nBody,nHeart,NTSM),d_TAU(nBody,NTSM),
     '  PHIDIFF(nBody,NTSM),YP(NYM,NIYM)

!     Local Variables
      INTEGER iB,nts,ny,iH
      REAL*8 dTAU,RNORM

      RNORM = 0.0d0
      CALC_RNORM = 0.0d0
      DO nts = S_START,S_END
        DO iB = 1,nBody
          dTAU = 0.0d0
          DO iH = 1,nHeart
            ny = NYNR(iH)
            dTAU = dTAU + dPHIdTAU(iB,iH,nts)*YP(ny,1)
          ENDDO
          RNORM = RNORM + ( dTAU - PHIDIFF(iB,nts) - d_TAU(iB,nts) )**2
        ENDDO
      ENDDO !nts

      CALC_RNORM = SQRT(RNORM)
      RETURN
      END


