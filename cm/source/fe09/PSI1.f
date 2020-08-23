      REAL*8 FUNCTION PSI1(IBT,IDO,INP,nb,nu,nk,nn,XI)

C#### Function: PSI1
C###  Type: REAL*8
C###  Description:
C###    PSI1 evaluates tensor product Lagrange/Hermite/Fourier basis
C###    functions at XI.

C**** IPU(nu,ni),nu=1,NUT(nb) identifies the complete set of partial
C**** derivatives with respect to Xi(ni).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,nk,nn,nu
      REAL*8 XI(3)
!     Local Variables
      INTEGER IPU(11,3),ni
      REAL*8 PF1,PH2,PH3,PL1,PL2,PL3

      DATA IPU/1,2,3,1,1,2,1,1,2,1,2,
     '         1,1,1,2,3,2,1,1,1,2,2,
     '         1,1,1,1,1,1,2,3,2,2,2/

      PSI1=1.0d0
      DO ni=1,NIT(nb)
        IF(IBT(1,ni).EQ.1) THEN
          IF(IBT(2,ni).EQ.1) PSI1=PSI1*PL1(INP(nn,ni),IPU(nu,ni),XI(ni))
          IF(IBT(2,ni).EQ.2) PSI1=PSI1*PL2(INP(nn,ni),IPU(nu,ni),XI(ni))
          IF(IBT(2,ni).EQ.3) PSI1=PSI1*PL3(INP(nn,ni),IPU(nu,ni),XI(ni))
        ELSE IF(IBT(1,ni).EQ.2) THEN
          IF(IBT(2,ni).EQ.1) THEN
            PSI1=PSI1*PH3(INP(nn,ni),IDO(nk,nn,ni),IPU(nu,ni),XI(ni))
          ELSE
            PSI1=PSI1*PH2(INP(nn,ni),IDO(nk,nn,ni),IPU(nu,ni),
     '        IBT(2,ni)-1,XI(ni))
          ENDIF
        ELSE IF(IBT(1,ni).EQ.9) THEN
          PSI1=PSI1*PF1(IDO(nk,nn,ni),IPU(nu,ni),XI(ni))
        ENDIF
      ENDDO

      RETURN
      END


