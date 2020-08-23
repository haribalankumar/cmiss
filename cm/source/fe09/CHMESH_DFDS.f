      REAL*8 FUNCTION CHMESH_DFDS(ZETA,LENGTH,XP,np,nk)

C#### Function: CHMESH_DFDS
C###  Type: REAL*8
C###  Description:
C###    CHMESH_DFDS is one of J Crocombe's (temporary) change
C###    mesh functions.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      REAL*8 ZETA,LENGTH,XP(NKM,NVM,NJM,NPM)
      INTEGER nk,np
      REAL*8 CHMESH_FDERIV

      CHMESH_DFDS = (CHMESH_FDERIV(ZETA)/LENGTH)*XP(nk,1,3,np)
      RETURN
      END


