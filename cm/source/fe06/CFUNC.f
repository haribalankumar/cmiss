      SUBROUTINE CFUNC(CONTYP,IBT,ICOORD,IDO,INP,NAN,NBH,NBJ,
     '  nh,NHE,NICONT,nj,nr,nx,
     '  PG,XE,XG,XIC,XICONT,ZE,ZG,U,DUXI,ERROR,*)

C#### Subroutine: CFUNC
C###  Description:
C###    CFUNC provides the field value U and its derivatives DUXI
C###    wrt the Xi coordinates.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),ICOORD,IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),
     '  nh,NHE,NICONT,nj,nr,nx
      REAL*8 DUXI(2),PG(NSM,NUM,NGM,NBM),U,XE(NSM,NJM),
     '  XG(NJM,NUM),XIC(4),XICONT,ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER CONTYP*(*),ERROR*(*)
!     Local Variables
      INTEGER i,ixi,j,j1,nb,NI(4),NKI(3)
      REAL*8 DX1XI(3,2),DX2XI(3,2),DXI,DXIXN(3,3),DXNXI(3,3),DXXI(3,2),
     '  DXZ,DZX,DZXI(3,2),EG(3,3),PHI(3),PST(3),PXI,
     '  RGX,RM(3,3),X(3),X1(3),X2(3),XI(3),
     '  Z(3)
      DATA NKI/2,4,7/

      CALL ENTERS('CFUNC',*9999)
      IF(NIT(NBH(nh)).EQ.2) THEN
        NI(1)=1
        NI(2)=2
        XI(1)=XIC(1)
        XI(2)=XIC(2)
      ELSE IF(NIT(NBH(nh)).EQ.3) THEN
        NI(1)=MOD(NICONT,3)+1
        NI(2)=MOD(NICONT+1,3)+1
        NI(3)=MOD(NICONT+2,3)+1
        XI(NI(1))=XIC(1)
        XI(NI(2))=XIC(2)
        XI(NI(3))=XICONT
      ENDIF
C     IF(DOP) WRITE(IO4,'('' NICONT='',I3,'' Xi: '',3F5.2)')
C    '  NICONT,XI(1),XI(2),XI(3)

      IF(CONTYP(1:5).EQ.'FIELD'.OR.CONTYP(1:9).EQ.'DEPENDENT') THEN
        nb=NBH(nh)
        U=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '    1,XI,ZE(1,nh))
        DUXI(1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '    NKI(NI(1)),XI,ZE(1,nh))
        DUXI(2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '    NKI(NI(2)),XI,ZE(1,nh))

      ELSE IF(CONTYP(1:9).EQ.'GEOMETRIC') THEN
        IF(ICOORD.EQ.ITYP10(1)) THEN
          nb=NBH(nh)
          U=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      1,XI,XE(1,nj))
          DUXI(1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      NKI(NI(1)),XI,XE(1,nj))
          DUXI(2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      NKI(NI(2)),XI,XE(1,nj))
        ELSE IF(ICOORD.EQ.1) THEN
          DO j=1,NJT
            nb=NBJ(j)
            X(j)     =PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        1,XI,XE(1,j))
            DXXI(j,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NKI(NI(1)),XI,XE(1,j))
            DXXI(j,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NKI(NI(2)),XI,XE(1,j))
          ENDDO
          CALL XZ(ITYP10(1),X,Z)
          DO i=1,2
            DZXI(nj,i)=0.d0
            DO j=1,NJT
              DZXI(nj,i)=DZXI(nj,i)+DZX(ITYP10(1),nj,j,X)*DXXI(j,i)
            ENDDO
          ENDDO
          U=Z(nj)
          DUXI(1)=DZXI(nj,1)
          DUXI(2)=DZXI(nj,2)
        ELSE IF(ITYP11(1).EQ.1) THEN
          DO j=1,NJT
            nb=NBJ(j)
            Z(j)     =PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        1,XI,XE(1,j))
            DZXI(j,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NKI(NI(1)),XI,XE(1,j))
            DZXI(j,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NKI(NI(2)),XI,XE(1,j))
          ENDDO
          CALL ZX(ICOORD,Z,X)
          DO i=1,2
            DXXI(nj,i)=0.d0
            DO j=1,NJT
              DXXI(nj,i)=DXXI(nj,i)+DXZ(ICOORD,nj,j,X)*DZXI(j,i)
            ENDDO
          ENDDO
          U=X(nj)
          DUXI(1)=DXXI(nj,1)
          DUXI(2)=DXXI(nj,2)
        ELSE
          DO j=1,NJT
            nb=NBJ(j)
            X1(j)     =PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        1,XI,XE(1,j))
            DX1XI(j,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NKI(NI(1)),XI,XE(1,j))
            DX1XI(j,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NKI(NI(2)),XI,XE(1,j))
          ENDDO
          CALL COORD(ITYP10(1),ICOORD,X1,X2,ERROR,*9999)
          DO i=1,2
            DX2XI(nj,i)=0.d0
            DO j=1,NJT
              DO j1=1,NJT
                DX2XI(nj,i)=DX2XI(nj,i)+DXZ(ICOORD,nj,j,X2)
     '            *DZX(ITYP10(1),j,j1,X1)*DX1XI(j1,i)
              ENDDO
            ENDDO
          ENDDO
          U=X2(nj)
          DUXI(1)=DX2XI(nj,1)
          DUXI(2)=DX2XI(nj,2)
        ENDIF

      ELSE IF(CONTYP(1:6).EQ.'STRAIN') THEN
         CALL ZEEX51(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '     DXIXN,DXNXI,EG,PG,PHI,PST,RGX,RM,XE,XG,XI,ZE,ZG,
     '    ERROR,*9999)
         U=PST(nh)
         DO ixi=1,2
           IF(XIC(ixi).LE.0.9999D0) THEN
             DXI= 0.0001D0
           ELSE
             DXI=-0.0001D0
           ENDIF
           XI(ixi)=XIC(ixi)+DXI
           CALL ZEEX51(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '       DXIXN,DXNXI,EG,PG,PHI,PST,RGX,RM,XE,XG,XI,ZE,ZG,
     '       ERROR,*9999)
           DUXI(ixi)=(PST(nh)-U)/DXI
           XI(ixi)=XIC(ixi)
         ENDDO

      ELSE IF(CONTYP(1:6).EQ.'CAUCHY') THEN

      ELSE IF(CONTYP(1:6).EQ.'VECTOR') THEN
        U=0.d0
        DUXI(1)=0.d0
        DUXI(2)=0.d0
        DO j=1,nj
          nb=NBJ(j)
          Z(j)     =PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      1,XI,ZE(1,nh+j))
          DZXI(j,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      NKI(NI(1)),XI,ZE(1,nh+j))
          DZXI(j,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '      NKI(NI(2)),XI,ZE(1,nh+j))
          U=U+Z(j)*Z(j)
          DUXI(1)=DUXI(1)+2.d0*Z(j)*DZXI(j,1)
          DUXI(2)=DUXI(2)+2.d0*Z(j)*DZXI(j,2)
        ENDDO
      ENDIF

      CALL EXITS('CFUNC')
      RETURN
 9999 CALL ERRORS('CFUNC',ERROR)
      CALL EXITS('CFUNC')
      RETURN 1
      END


