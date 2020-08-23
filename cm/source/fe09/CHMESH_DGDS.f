      REAL*8 FUNCTION CHMESH_DGDS(THETA,XP,np,nk,TEMP)

C#### Function: CHMESH_DGDS
C###  Type: REAL*8
C###  Description:
C###    CHMESH_DGDS is one of J Crocombe's (temporary) change
C###    mesh functions.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      REAL*8 THETA,TEMP,XP(NKM,NVM,NJM,NPM)
      INTEGER nk,np
      REAL*8 CHMESH_GDERIV
      CHMESH_DGDS =
     '  CHMESH_GDERIV(THETA)*
     '  ((XP(nk,1,2,np)*(XP(1,1,1,np)/(XP(1,1,1,np)**2
     '  +XP(1,1,2,np)**2)))-(TEMP*(XP(1,1,2,np)/(XP(1,1,1,np)**2
     '  +XP(1,1,2,np)**2))))
      RETURN
      END


