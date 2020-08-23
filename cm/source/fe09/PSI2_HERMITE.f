      REAL*8 FUNCTION PSI2_HERMITE(IDO,INP,nb,nu,nk,nn,XI)

C#### Function: PSI2_HERMITE
C###  Type: REAL*8
C###  Description:
C###    PSI2_HERMITE evaluates Hermite simplex basis function.
C**** IPU(nu,ni),nu=1,NUT(nb) identifies the complete set of partial
C**** derivatives with respect to Xi(ni) up to second order.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,nk,nn,nu
      REAL*8 XI(3)
!     Local Variables
      INTEGER IPU(11,3)
      REAL*8 PH3,PL2S1,PL2S3

      DATA IPU/1,2,3,1,1,2,1,1,2,1,2,
     '         1,1,1,2,3,2,1,1,1,2,2,
     '         1,1,1,1,1,1,2,3,2,2,2/

      PSI2_HERMITE=1.0d0
      IF(NIT(nb).EQ.2)THEN
!       !Basis function is hermite in xi1 and special quadratic in xi2
        IF(NKT(1,nb).EQ.1)THEN !Apex at node 1
          IF(nn.GT.1)THEN
            PSI2_HERMITE=PSI2_HERMITE*
     '         PH3(INP(nn,1),IDO(nk,nn,1),IPU(nu,1),XI(1))
          ENDIF
          IF((nn.EQ.1).AND.(IPU(nu,1).GT.1))THEN
            PSI2_HERMITE=0.0D0 !Deriv at node 1 in xi1 direction = 0.
          ELSE
            PSI2_HERMITE=PSI2_HERMITE*
     '         PL2S1(INP(nn,2),IDO(nk,nn,2),IPU(nu,2),XI(2))
          ENDIF
        ELSEIF(NKT(3,nb).EQ.1)THEN !Apex at node 3
          IF(nn.ne.3)THEN
            PSI2_HERMITE=PSI2_HERMITE*
     '         PH3(INP(nn,1),IDO(nk,nn,1),IPU(nu,1),XI(1))
          ENDIF
          IF((nn.EQ.3).AND.(IPU(nu,1).GT.1))THEN
            PSI2_HERMITE=0.0D0 !Deriv at node 3 in xi1 direction = 0.
          ELSE
            PSI2_HERMITE=PSI2_HERMITE*
     '         PL2S3(INP(nn,2),IDO(nk,nn,2),IPU(nu,2),XI(2))
          ENDIF
        ENDIF
      ENDIF
      IF(NIT(nb).EQ.3) THEN
      ENDIF

      RETURN
      END


