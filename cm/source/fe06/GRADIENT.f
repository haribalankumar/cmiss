      SUBROUTINE GRADIENT(INDEX,NBH,NBJ,ne,NHE,nr,nx,
     '  PG,RG,XE,XG,ZE,ZG,ERROR,*)

C#### Subroutine: GRADIENT
C###  Description:
C###    GRADIENT draws gradient vectors of scalar field at Gauss points.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grad00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER INDEX,NBH(NHM,NCM),NBJ(NJM),ne,NHE,nr,nx
      REAL*8 PG(NSM,NUM,NGM,NBM),RG(NGM),XE(NSM,NJM),
     '  XG(NJM,NUM),ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nc,ng,ni,nj
      REAL*8 DXIX(3,3),GL(3,3),GU(3,3),PL(3,2)

      CALL ENTERS('GRADIENT',*9999)
      nc=1 ! Temporary
      IF(ITYP1(nr,nx).EQ.3) THEN
        nb=NBH(NH_LOC(1,nx),nc)
        DO ng=1,NGT(nb)
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Element '',I3,'' Gauss pt '',I2)')
     '        ne,ng
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
C Note use of IP=0 in call to XGMG (DXIX wrt Xi coords)
          CALL XGMG(0,0,NBJ(1),nr,DXIX,GL,GU,RG(ng),XG,
     '      ERROR,*9999)
          IF(DOP) THEN
            DO ni=1,NIT(nb)
              WRITE(OP_STRING,'('' DXIX('',I1,'',nj): '',3E11.3)')
     '          ni,(DXIX(ni,nj),nj=1,NJT)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDIF
          CALL ZEZG(1,NBH,ng,NHE,nx,DXIX,PG,ZE,ZG,ERROR,*9999)
          IF(DOP) THEN
            WRITE(OP_STRING,'('' ZG(1,2)='',E11.3,'
     '        //''' ZG(1,4)='',E11.3)') ZG(1,2),ZG(1,4)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          PL(1,1)=XG(1,1)-ZG(1,2)*SCALE
          PL(2,1)=XG(2,1)-ZG(1,4)*SCALE
          PL(1,2)=XG(1,1)+ZG(1,2)*SCALE
          PL(2,2)=XG(2,1)+ZG(1,4)*SCALE
          CALL POLYLINE(INDEX,1,2,PL,ERROR,*9999) !assume iw=1 AAY
        ENDDO
      ENDIF

      CALL EXITS('GRADIENT')
      RETURN
 9999 CALL ERRORS('GRADIENT',ERROR)
      CALL EXITS('GRADIENT')
      RETURN 1
      END


