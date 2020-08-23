      SUBROUTINE STRAIN2(INDEX,IBT,IDO,INP,IPOINTTYP,IW,LD,NAN,
     '  NBH,NBJ,ne,NHE,nr,nx,PG,RG,XE,XG,XID,XIG,ZE,ZG,ERROR,*)

C#### Subroutine: STRAIN2
C###  Description:
C###    STRAIN2 is for principal strain output.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),IPOINTTYP,IW,LD(NDM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM),NBJ(NJM),ne,NHE,nr,nx
      REAL*8 PG(NSM,NUM,NGM,NBM),RG(NGM),XE(NSM,NJM),
     '  XG(NJM,NUM),XID(NIM,NDM),XIG(NIM,NGM,NBM),ZE(NSM,NHM),
     '  ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER INDEX_PL,j,nb,ni,NITB,nj,NU1(0:3),num_points,
     '  point_counter
      REAL*8 DP,DX(3),DXIXN(3,3),DXNXI(3,3),
     '  EG(3,3),PHI(3),PST(3),RM(3,3),SUM,X(3),XY(3,2)
      DATA NU1/1,2,4,7/

      CALL ENTERS('STRAIN2',*9999)
      nb=NBH(1)
      NITB=NIT(nb)
      IF(IPOINTTYP.EQ.1) THEN
        num_points=NGT(nb)
      ELSE
        num_points=NDT
      ENDIF
      DO point_counter=1,num_points
        IF((IPOINTTYP.EQ.1).OR.
     '    ((IPOINTTYP.EQ.2).AND.(LD(point_counter).EQ.ne))) THEN
          IF(DOP) THEN
            IF(IPOINTTYP.EQ.1) THEN
              WRITE(OP_STRING,'(/''Element '',I3,'' Gauss point '',I2)')
     '          ne,point_counter
            ELSE
              WRITE(OP_STRING,'(/''Element '',I3,'' Data point '',I2)')
     '          ne,point_counter
            ENDIF
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(IPOINTTYP.EQ.1) THEN
            CALL ZEEX51(IBT,IDO,INP,NAN,NBH,NBJ,point_counter,NHE,nr,nx,
     '        DXIXN,DXNXI,EG,PG,PHI,PST,RG(point_counter),RM,XE,XG,
     '        XIG(1,point_counter,nb),ZE,ZG,ERROR,*9999)
          ELSE
            CALL ZEEX51(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE,nr,nx,
     '        DXIXN,DXNXI,EG,PG,PHI,PST,RG(1),RM,XE,XG,
     '        XID(1,point_counter),ZE,ZG,ERROR,*9999)
          ENDIF

C ***   May want this later if we draw vectors on undeformed geometry
C       DO nj=1,NJT
C         DO mi=1,NITB
C           SUM=0.d0
C           DO K=1,NITB
C             SUM=SUM+XG(nj,NU1(K))*DXIXN(K,mi)
C           ENDDO
C           DZJXN(nj,mi)=SUM
C         ENDDO
C       ENDDO

C ***  Princ. strain vectors are drawn at Gauss point Xi coords XIG
C ***  and projected into deformed geometry
          DO j=1,NITB
            DP=PST(j)*SCALE
            DO nj=1,NJT
              SUM=0.d0
              DO ni=1,NITB
                SUM=SUM+RM(ni,j)*ZG(nj,NU1(ni))
              ENDDO
              DX(nj)=SUM*DP
            ENDDO
            DO nj=1,NJT
              X(nj)=ZG(nj,1)-DX(nj)/2.d0
            ENDDO
            CALL XZ(ITYP10(1),X,XY(1,1))
            DO nj=1,NJT
              X(nj)=ZG(nj,1)+DX(nj)/2.d0
            ENDDO
            CALL XZ(ITYP10(1),X,XY(1,2))
            IF(DP.LT.0) THEN
              INDEX_PL=INDEX+1
            ELSE
              INDEX_PL=INDEX
            ENDIF
            CALL POLYLINE(INDEX_PL,IW,2,XY,ERROR,*9999)
          ENDDO
        ENDIF
      ENDDO

      CALL EXITS('STRAIN2')
      RETURN
 9999 CALL ERRORS('STRAIN2',ERROR)
      CALL EXITS('STRAIN2')
      RETURN 1
      END


