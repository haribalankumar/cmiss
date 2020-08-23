      SUBROUTINE FIBRE2(INDEX,IBT,IDO,INP,IW,NAN,NBJ,nr,XE,XG,
     '  ERROR,*)

C#### Subroutine: FIBRE2
C###  Description:
C###    FIBRE2 draws fitted fibre field vectors (DXIF) on curvilinear
C###    plane defined by constant Xi(3) in element ne at increments of
C###    DXI1 and DXI2 in Xi(1) and Xi(2) directions, respectively.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fibr00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER INDEX,IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IW,NAN(NIM,NAM,NBFM),NBJ(NJM),nr
      REAL*8 XE(NSM,NJM),XG(NJM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i1,i2,ITOT1,ITOT2,n,nb,nj,njj,ns
      REAL*8 DXINU(3,3),DX1,DX2,
     '  FLENGTH,GXL(3,3),GXU(3,3),PXI,RWX,SCALE,X(3),XI(3),XIP1,XIP2,
     '  XIP3,XX(2),YY(2),Z(3),Z3D(3,2),ZZ(2)

      CALL ENTERS('FIBRE2',*9999)

      IF(DOP) THEN
        DO njj=1,NJ_LOC(NJL_FIBR,0,nr)
          nj=NJ_LOC(NJL_FIBR,njj,nr)
          nb=NBJ(nj)
          WRITE(OP_STRING,'('' XE: '',8E11.3)')
     '      (XE(ns,nj),ns=1,NST(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF
      ITOT1=IDNINT(1.d0/DXI1+1.d-4)
      ITOT2=IDNINT(1.d0/DXI2+1.d-4)
      XIP3=XIF
      IF(DOP) THEN
        WRITE(OP_STRING,'('' XIP3='',E12.4)') XIP3
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      XIP1=-DXI1/2.d0
      DO i1=1,ITOT1
        XIP1=XIP1+DXI1
        XIP2=-DXI2/2.d0
        DO i2=1,ITOT2
          XIP2=XIP2+DXI2
          XI(1)=XIP1
          XI(2)=XIP2
          XI(3)=XIP3
C new MPN 30-Apr-96:
C         Interpolate geometric vars XG and derivs wrt Xi
          CALL XEXW(0,IBT,IDO,INP,NAN,NBJ,nr,XE,XG,XI,ERROR,*9999)
C         Calculate undeformed metric tensors wrt Nu (GXL,GXU) and
C         derivs (DXINU) of Xi wrt Nu  coords.
          CALL XGMG(1,0,NBJ(1),nr,DXINU,GXL,GXU,RWX,XG,ERROR,*9999)

C LC 24/2/97 archived section :
C      old MPN 30-Apr-96: old way of calculating material axes

C!!!      This scale factor seem dodgy!?
          SCALE=DSQRT(GXL(1,1))
          DX1=DXIF*SCALE*DXINU(1,1)
          DX2=DXIF*SCALE*DXINU(2,1)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' DXIF='',D11.3,'' DX1='',D11.3,'
     '        //''' DX2='',D11.3)') DXIF,DX1,DX2
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          XI(1)=XIP1-DX1/2.0d0
          XI(2)=XIP2-DX2/2.0d0
          IF(DOP) THEN
            WRITE(OP_STRING,'('' XI(1)='',D11.3,'' XI(2)='',D11.3)')
     '        XI(1),XI(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(PROJEC(1:2).EQ.'XI') THEN
            Z(1)=XI(1)
            Z(2)=XI(2)
          ELSE
            DO nj=1,NJT
              nb=NBJ(nj)
              X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '          XE(1,nj))
            ENDDO
            IF(DOP) THEN
              WRITE(OP_STRING,'('' X(1)='',D11.3,'' X(2)='',D11.3)')
     '          X(1),X(2)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(IW.EQ.4) THEN
              Z(2)=X(2)
              Z(3)=X(3)
            ELSE
              CALL XZ(ITYP10(1),X,Z)
            ENDIF
          ENDIF
          XX(1)=Z(1)
          YY(1)=Z(2)
          ZZ(1)=Z(3)
          XI(1)=XIP1+DX1/2.d0
          XI(2)=XIP2+DX2/2.d0

          IF(PROJEC(1:2).EQ.'XI') THEN
            Z(1)=XI(1)
            Z(2)=XI(2)
          ELSE
            DO nj=1,NJT
              nb=NBJ(nj)
              X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '          XE(1,nj))
            ENDDO
            IF(IW.EQ.4) THEN
              Z(2)=X(2)
              Z(3)=X(3)
            ELSE
              CALL XZ(ITYP10(1),X,Z)
            ENDIF
          ENDIF
          XX(2)=Z(1)
          YY(2)=Z(2)
          ZZ(2)=Z(3)
          DO n=1,2
            Z3D(1,n)=XX(n)
            Z3D(2,n)=YY(n)
            Z3D(3,n)=ZZ(n)
          ENDDO
          IF(IW.EQ.4.AND.PROJEC(1:6).EQ.'HAMMER') THEN
            FLENGTH=DSQRT((XX(2)-XX(1))**2+(YY(2)-YY(1))**2)
            IF(FLENGTH.LT.0.2D0) THEN
              CALL POLYLINE(INDEX,IW,2,Z3D,ERROR,*9999)
            ENDIF
          ELSE
            CALL POLYLINE(INDEX,IW,2,Z3D,ERROR,*9999)
          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('FIBRE2')
      RETURN
 9999 CALL ERRORS('FIBRE2',ERROR)
      CALL EXITS('FIBRE2')
      RETURN 1
      END


