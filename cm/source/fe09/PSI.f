      REAL*8 FUNCTION PSI(IBT,IDO,INP,nb,nn,nk,nu,XI)

C#### Function: PSI
C###  Type: REAL*8
C###  Description:
C###    PSI evaluates a basis funtion at XI coords with basis
C###    functions as originally defined (jtyp5=1)

C**** Created by Carey Stevens 20 April 1998

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,nu
      REAL*8 XI(3)
!     Local Variables
      INTEGER ni,nk,nn
      REAL*8 PSI1,PSI2,PSI2_HERMITE,PSI2_XI,PSI5,XL(4)
      LOGICAL SECTOR

      SECTOR=.FALSE.
      DO ni=1,NIT(nb)
        IF(IBT(1,ni).EQ.5.OR.IBT(1,ni).EQ.6) SECTOR=.TRUE.
      ENDDO !ni
      IF(SECTOR) THEN
        PSI=PSI5(IBT,IDO,INP,nb,nu,nk,nn,XI)
      ELSE
        IF(IBT(1,1).LE.2.OR.IBT(1,1).EQ.9) THEN
          PSI=PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)
        ELSE IF(IBT(1,1).EQ.3.AND.IBT(2,1).EQ.4) THEN
          PSI=PSI2_HERMITE(IDO,INP,nb,nu,nk,nn,XI)
        ELSE IF(IBT(1,1).EQ.3.AND.IBT(3,1).EQ.2) THEN
          PSI=PSI2_XI(nb,nu,nn,XI)
        ELSE IF(IBT(1,1).EQ.3) THEN
          XL(2)=XI(1)
          XL(3)=XI(2)
          XL(4)=XI(3)
          XL(1)=1.0D0-XL(2)-XL(3)-XL(4)
          PSI=PSI2(IBT,INP,nb,nu,nn,XL)
        ENDIF
      ENDIF

      RETURN
      END


