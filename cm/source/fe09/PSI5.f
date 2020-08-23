      REAL*8 FUNCTION PSI5(IBT,IDO,INP,nb,nu,nk,nn,XI)

C#### Function: PSI5
C###  Type: REAL*8
C###  Description:
C###    PSI5 evaluates Sector basis functions at coordinates XI.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,nk,nn,nu
      REAL*8 XI(3)
!     Local Variables
      INTEGER IPU(11,3),n,ni
      REAL*8 PH2,PH3,PL1,PL2,PL3,SUM
      LOGICAL COLLAPSEDNODE

      DATA IPU/1,2,3,1,1,2,1,1,2,1,2,
     '         1,1,1,2,3,2,1,1,1,2,2,
     '         1,1,1,1,1,1,2,3,2,2,2/

c cpb 20/7/95 Generalising for any sector

C***  Work out if local node nn is a collapsed node
      COLLAPSEDNODE=.FALSE.
      DO ni=1,NIT(nb)
        IF(IBT(1,ni).EQ.5) THEN
          IF(INP(nn,IBT(3,ni)).EQ.1) COLLAPSEDNODE=.TRUE.
        ELSE IF(IBT(1,ni).EQ.6) THEN
          IF(IBT(2,ni).EQ.4.AND.INP(nn,IBT(3,ni)).EQ.2) THEN
            COLLAPSEDNODE=.TRUE.
          ELSE IF(INP(nn,IBT(3,ni)).EQ.(IBT(2,IBT(3,ni))+1)) THEN
            COLLAPSEDNODE=.TRUE.
          ENDIF
        ENDIF
      ENDDO
C      IF(DOP) THEN
C        WRITE(OP_STRING,'('' nb='',I2,'' nu='',I2,'' nn='',I2,'
C     '    //''' nk='',I2)') nb,nu,nn,nk
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' COLLAPSED NODE='',L1)') COLLAPSEDNODE
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF
C**** Evaluate the basis functions
      PSI5=1.0d0
      DO ni=1,NIT(nb)
        IF(IBT(1,ni).EQ.1) THEN
          IF(IBT(2,ni).EQ.1) THEN
            PSI5=PSI5*PL1(INP(nn,ni),IPU(nu,ni),XI(ni))
          ELSE IF(IBT(2,ni).EQ.2) THEN
            PSI5=PSI5*PL2(INP(nn,ni),IPU(nu,ni),XI(ni))
          ELSE IF(IBT(2,ni).EQ.3) THEN
            PSI5=PSI5*PL3(INP(nn,ni),IPU(nu,ni),XI(ni))
          ENDIF
        ELSE IF(IBT(1,ni).EQ.2) THEN
C cpb 10/6/97 There are two ways to do a colapsed cubic Hermite.
C In both cases the product of the basis function times the derivative
C at the collapsed node is zero. The differences come about because
C either the basis function is zero or the derivative is zero. In the
C first case you get a "quadratic Hermite" basis and in the second case
C you get the normal "cubic Hermite" basis.

CC cpb 18/2/98 Use cubic Hermite
C cpb 10/6/97 Use quadratic for now
C          PSI5=PSI5*PH3(INP(nn,ni),IDO(nk,nn,ni),IPU(nu,ni),XI(ni))
          IF(IBT(2,ni).EQ.1) THEN
            PSI5=PSI5*PH3(INP(nn,ni),IDO(nk,nn,ni),IPU(nu,ni),XI(ni))
          ELSE
            PSI5=PSI5*PH2(INP(nn,ni),IDO(nk,nn,ni),IPU(nu,ni),
     '        IBT(2,ni)-1,XI(ni))
          ENDIF
        ELSE IF(IBT(1,ni).EQ.5.OR.IBT(1,ni).EQ.6) THEN
          IF(COLLAPSEDNODE) THEN
C****       If the local node is a collapsed node then we have to
C****       sum the basis functions in the collapsed xi direction
            SUM=0.0d0
            IF(IBT(2,ni).EQ.1) THEN
              DO n=1,2
                SUM=SUM+PL1(n,IPU(nu,ni),XI(ni))
              ENDDO !n
            ELSE IF(IBT(2,ni).EQ.2) THEN
              DO n=1,3
                SUM=SUM+PL2(n,IPU(nu,ni),XI(ni))
              ENDDO !n
            ELSE IF(IBT(2,ni).EQ.3) THEN
              DO n=1,4
                SUM=SUM+PL3(n,IPU(nu,ni),XI(ni))
              ENDDO !n
            ELSE IF(IBT(2,ni).EQ.4) THEN
              DO n=1,2
                SUM=SUM+PH3(n,IDO(nk,nn,ni),IPU(nu,ni),XI(ni))
              ENDDO !n
            ENDIF
            PSI5=PSI5*SUM
          ELSE
            IF(IBT(2,ni).EQ.1) THEN
              PSI5=PSI5*PL1(INP(nn,ni),IPU(nu,ni),XI(ni))
            ELSE IF(IBT(2,ni).EQ.2) THEN
              PSI5=PSI5*PL2(INP(nn,ni),IPU(nu,ni),XI(ni))
            ELSE IF(IBT(2,ni).EQ.3) THEN
              PSI5=PSI5*PL3(INP(nn,ni),IPU(nu,ni),XI(ni))
            ELSE IF(IBT(2,ni).EQ.4) THEN
              PSI5=PSI5*PH3(INP(nn,ni),IDO(nk,nn,ni),IPU(nu,ni),XI(ni))
            ENDIF
          ENDIF
        ENDIF
      ENDDO !ni

C      IF(DOP) THEN
C        WRITE(OP_STRING,'('' XI: '',3D12.4)') (XI(ni),ni=1,NIT(nb))
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,'('' PSI5: '',D12.4)') PSI5
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C      ENDIF

      RETURN
      END


