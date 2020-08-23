      SUBROUTINE FIBRE1(INDEX,IBT,IDO,INP,IW,NBJ,ne,NKJE,NPF,
     '  NPNE,NRE,NVJE,ALFA,DXIF,SE,XA,XE,XI,XP,ERROR,*)

C#### Subroutine: FIBRE1
C###  Description:
C###    FIBRE1 draws measured fibre field vector (DXIF) on curvilinear
C###    plane defined by constant Xi(3) in element ne at Xi(ni).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'scal01.cmn'
      INCLUDE 'trans00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),IW,NBJ(NJM,NEM),ne,
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 ALFA,DXIF,SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XI(NIM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER n,nb,ni,nj,njj,nk,nn,ns,nr
      REAL*8 DX1,DX2,FLENGTH,GL(3,3),
     '  GU(3,3),PXI,RG33,RL,RK1,RK2,X(3),XIP1,XIP2,
     '  XRC(2),YRC(2),Z(3),Z3D(3,2),ZRC(2)

      CALL ENTERS('FIBRE1',*9999)

      nr=NRE(ne)

      CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '  nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

      IF(IW.EQ.4.AND.ITYP10(nr).EQ.4) THEN !hammer proj with prol coords
        nb=NBJ(1,ne)
        DO nn=1,NNT(nb) !set lamda=1 and zero derivs
          XE(1+(nn-1)*NKT(0,nb),1)=1.d0
          DO nk=2,NKT(0,nb)
            XE(nk+(nn-1)*NKT(0,nb),1)=0.d0
          ENDDO
        ENDDO
      ENDIF
      IF(DOP) THEN
        DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
          nj=NJ_LOC(NJL_FIBR,njj,nr)
          nb=NBJ(nj,ne)
          WRITE(OP_STRING,'('' XE: '',8E11.3)')
     '      (XE(ns,nj),ns=1,NST(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'('' Xi:'',3E11.3)')
     '    (XI(ni),ni=1,NIT(NBJ(1,ne)))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL XMG(IBT,IDO,INP,NBJ(1,ne),nr,
     '  GL,GU,XE,XI,ERROR,*9999)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' GL:'',3E11.3)') GL(1,1),GL(2,2),GL(1,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      RL=DSQRT(GL(1,1)*GL(2,2)-GL(1,2)**2)
      RK1=DSQRT(GL(1,1))*DCOS(ALFA)
      RK2=(GL(1,2)*RK1+DSQRT(GL(1,1))*RL*DSIN(ALFA))/GL(1,1)
      IF(DOP) THEN
        WRITE(OP_STRING,'('' RL='',E11.3,'' RK1='',E11.3,'' RK2='','
     '    //'E11.3)') RL,RK1,RK2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(ITYP10(1).EQ.4) THEN
        RG33=DSQRT(GL(1,1))
      ELSE
        RG33=1.d0
      ENDIF
      DX1=SCALE*DXIF*RG33*(GL(2,2)*RK1-GL(1,2)*RK2)/RL**2
      DX2=SCALE*DXIF*RG33*(GL(1,1)*RK2-GL(1,2)*RK1)/RL**2
      IF(DOP) THEN
        WRITE(OP_STRING,'('' DXIF='',E11.3,'' DX1='',E11.3,'
     '    //''' DX2='',E11.3)') DXIF,DX1,DX2
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      XIP1=XI(1)
      XIP2=XI(2)
      XI(1)=XIP1-DX1/2.d0
      XI(2)=XIP2-DX2/2.d0
      DO nj=1,NJT
        nb=NBJ(nj,ne)
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '    XE(1,nj))
      ENDDO
      IF(IW.LT.3) THEN
        CALL XZ(ITYP10(1),X,Z)
      ELSE IF(IW.EQ.3) THEN
        CALL XZ(ITYP10(1),X,Z)
        CALL ZZ(Z,Z,TRANS)
      ELSE IF(IW.EQ.4) THEN
        Z(2)=X(2)
        Z(3)=X(3)
      ENDIF
      XRC(1)=Z(1)
      YRC(1)=Z(2)
      ZRC(1)=Z(3)
      XI(1)=XIP1+DX1/2.d0
      XI(2)=XIP2+DX2/2.d0
      DO nj=1,NJT
        nb=NBJ(nj,ne)
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '    XE(1,nj))
      ENDDO
      IF(IW.LT.3) THEN
        CALL XZ(ITYP10(1),X,Z)
      ELSE IF(IW.EQ.3) THEN
        CALL XZ(ITYP10(1),X,Z)
        CALL ZZ(Z,Z,TRANS)
      ELSE IF(IW.EQ.4) THEN
        Z(2)=X(2)
        Z(3)=X(3)
      ENDIF
      XRC(2)=Z(1)
      YRC(2)=Z(2)
      ZRC(2)=Z(3)
      DO n=1,2
        Z3D(1,n)=XRC(n)
        Z3D(2,n)=YRC(n)
        Z3D(3,n)=ZRC(n)
      ENDDO
      IF(IW.LE.3) THEN
        CALL POLYLINE(INDEX,IW,2,Z3D,ERROR,*9999)
      ELSE IF(IW.EQ.4) THEN
        FLENGTH=DSQRT((XRC(2)-XRC(1))**2+(YRC(2)-YRC(1))**2)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' FLENGTH='',E12.3)') FLENGTH
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        IF(FLENGTH.LT.0.2d0) THEN
          CALL POLYLINE(INDEX,IW,2,Z3D,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('FIBRE1')
      RETURN
 9999 CALL ERRORS('FIBRE1',ERROR)
      CALL EXITS('FIBRE1')
      RETURN 1
      END


