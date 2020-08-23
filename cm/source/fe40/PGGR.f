      SUBROUTINE PGGR(NBH,NBJ,NW,GG,PG,RG,RR,XE,XG,X3G,ERROR,*)

C#### Subroutine: PGGR
C###  Description:
C###    PGGR finds the strain functions (GG) and curvature change
C###    function (RR) of the basis functions PG for all Gauss points
C###    of the current element ne.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NBH(NHM),NBJ(NJM),NW
      REAL*8 GG(2,2,48,*),PG(NSM,NUM,NGM,NBM),RG(NGM),RR(2,2,48,*),
     '  XE(NSM,NJM),XG(NJM,NUM),X3G(4,3,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ia,ib,il,im,nb,ng,nh,nhs,nhx,nr,ns,nu,NU1,NU2,
     '  nvar,NV1(3),NV2(3,3),nx
      REAL*8 BM(3,3,25),BL(3,3,25),CHG(3,3,3,25),CHTOFF(3,3,3),
     '  DBM(3,3,3,25),DXIX(3,3),GL(3,3),GU(3,3),GUG(3,3,25),
     '  SUM,SUM1,SUM2

      DATA NV1/2,4,7/
      DATA NV2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('PGGR',*9999)
C GBS 10-NOV-1994
      nr=1   !Temporary
      nx=1 !temporary
      DO NG=1,NGT(NBJ(1))
        CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
        CALL XGMG(0,NIT(NBJ(1)),NBJ(1),nr,DXIX,GL,GU,
     '    RG(ng),XG,ERROR,*9999)
        CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM(1,1,1,ng),GU,
     '    XG,X3G(1,1,ng),.TRUE.,ERROR,*9999)
        DO ia=1,2
          DO ib=1,2
            GUG(ia,ib,ng)=GU(ia,ib)
            BL(ia,ib,ng)=CHTOFF(3,ia,ib)
            SUM=0.0d0
            DO il=1,2
              SUM=SUM+GU(ia,il)*CHTOFF(3,il,ib)
              CHG(ia,ib,il,ng)=CHTOFF(ia,ib,il)
            ENDDO
            BM(ia,ib,ng)=SUM
            IF(DABS(BL(ia,ib,ng)).LT.1.0D-8) BL(ia,ib,ng)=0.d0
            IF(DABS(BM(ia,ib,ng)).LT.1.0D-8) BM(ia,ib,ng)=0.d0
          ENDDO
        ENDDO
      ENDDO

      nhs=0
      DO nvar=1,NVE(NW)
        nhx=NHV(nvar,NW)
        nh=NH_LOC(nhx,nx)
        nb=NBH(nh)
        DO ns=1,NST(nb)
          nhs=nhs+1
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,
     '        '('' nh='',I2,'' ns='',I3,'' nhs='',I3)') nh,ns,nhs
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          DO ng=1,NGT(nb)
            DO 50 ia=1,2
            DO 50 ib=1,2
              SUM1=0.0d0
              SUM2=0.0d0
              IF(nh.LE.2) THEN !in-plane displacements
                IF(nh.EQ.ib) THEN
                  DO il=1,2
                    nu=NV1(il)
                    SUM1=SUM1+0.50d0*GUG(ia,il,ng)*PG(ns,nu,ng,nb)
                  ENDDO
                ENDIF
                NU1=NV1(ib)
                SUM1=SUM1+0.50d0*GUG(ia,nh,ng)*PG(ns,NU1,ng,nb)
                DO il=1,2
                  SUM1=SUM1-0.50d0*GUG(ia,il,ng)*(CHG(nh,ib,il,ng)
     '                     +CHG(nh,il,ib,ng))*PG(ns,1,ng,nb)
                  NU2=NV1(il)
C                          DBM(nh,ib,il,ng)=0.d0
                  SUM2=SUM2+GUG(ia,il,ng)*(
     '                     +DBM(nh,ib,il,ng)*PG(ns,1,ng,nb)
     '                     +BM(nh,ib,ng)*PG(ns,NU2,ng,nb)
     '                     +BM(nh,il,ng)*PG(ns,NU1,ng,nb))
                  DO im=1,2
                    SUM2=SUM2+GUG(ia,il,ng)*(-BM(im,ib,ng)
     '                       *CHG(nh,im,il,ng)-BM(im,il,ng)
     '                       *CHG(nh,im,ib,ng))*PG(ns,1,ng,nb)
                  ENDDO
                ENDDO
              ELSE IF(nh.EQ.3) THEN !transverse displacement
                DO il=1,2
                  SUM1=SUM1-GUG(ia,il,ng)*BL(il,ib,ng)*PG(ns,1,ng,nb)
                  nu=NV2(il,ib)
                  SUM2=SUM2+GUG(ia,il,ng)*PG(ns,nu,ng,nb)
                  DO im=1,2
                    SUM2=SUM2-GUG(ia,il,ng)*CHG(im,il,ib,ng)
     '                *PG(ns,NV1(im),ng,nb)
                    SUM2=SUM2-GUG(ia,il,ng)*BM(im,il,ng)*BL(im,ib,ng)
     '                *PG(ns,1,ng,nb)
                  ENDDO
                ENDDO
              ENDIF
              GG(ia,ib,nhs,ng)=SUM1
              RR(ia,ib,nhs,ng)=SUM2
              IF(DABS(SUM1).LT.1.0D-8) GG(ia,ib,nhs,ng)=0.d0
              IF(DABS(SUM2).LT.1.0D-8) RR(ia,ib,nhs,ng)=0.d0
 50         CONTINUE
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            DO ia=1,2
              DO ib=1,2
                WRITE(OP_STRING,'('' GG('',I1,'','',I1,'',nhs,ng)='','
     '            //'10(1X,E11.4))')
     '            ia,ib,(GG(ia,ib,nhs,ng),ng=1,NGT(nb))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' RR('',I1,'','',I1,'',nhs,ng)='','
     '            //'10(1X,E11.4))')
     '            ia,ib,(RR(ia,ib,nhs,ng),ng=1,NGT(nb))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDDO
CC$          call mp_unsetlock()
          ENDIF
        ENDDO
      ENDDO

      CALL EXITS('PGGR')
      RETURN
 9999 CALL ERRORS('PGGR',ERROR)
      CALL EXITS('PGGR')
      RETURN 1
      END


