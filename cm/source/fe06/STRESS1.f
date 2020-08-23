      SUBROUTINE STRESS1(INDEX,IBT,IDO,INP,IW,NAN,NBH,NBJ,ne,NHE,
     '  NPNE,nr,NW,nx,CE,CG,CP,FEXT,PG,RG,XE,XG,XIG,YG,ZE,ZG,ERROR,*)

C#### Subroutine: STRESS1
C###  Description:
C###    STRESS1 is for principal stress output.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbst02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),IW,NAN(NIM,NAM,NBFM),NBH(NHM,NCM),NBJ(NJM),
     '  ne,NHE,NPNE(NNM,NBFM),nr,NW,nx
      REAL*8 CE(NMM),CG(NMM,NGM),CP(NMM,NPM),FEXT(NIFEXTM,NGM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),XE(NSM,NJM),
     '  XG(NJM,NUM),XIG(NIM,NGM,NBM),YG(NIYGM,NGM),
     '  ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,nb,nc,ng,ni,nj,NJ_LOC00_temp,ns,NU1(3)
      REAL*8 AXL(3,3),AZ,AZL(3,3),AZU(3,3),DNU,DX(3),
     '  DXINU(3,3),DXDNU(3,3),EG(3,3),EGNU(3,3),
     '  EVAL(3),EVEC(3),G,GL(3,3),GU(3,3),PGX,PHI,PHII(3),PPG(9,2,16),
     '  RG2D,RGZ,RGZ2D,RM(3,3),RT(3,3),
     '  TC(3,3),TG(3,3),TGNU(3,3),TN(3,3),TNA,X(3),XY(3,2)
      CHARACTER STRESSTYPE*17
      LOGICAL TEMP_FIBRE

      DATA NU1/2,4,7/
      DATA STRESSTYPE/'Total'/

      CALL ENTERS('STRESS1',*9999)
      nc=1 ! Temporary

      IF(ITYP1(nr,nx).GE.4) THEN
        DO ng=1,NGT(NBH(NH_LOC(1,nx),nc))
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Element '',I3,'' Gauss pt '',I2)')
     '        ne,ng
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(ITYP1(nr,nx).EQ.4.AND.(NW.eq.5.or.NW.eq.11.or.NW.EQ.12).
     '      OR.ITYP1(nr,nx).EQ.5) THEN
C ***       Calculate derivs of Xi wrt undef Nu (fibre) coords
C ***       (DXINU), & metrics (GL,GU) wrt undef Nu coords
            CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
            CALL XGMG(1,0,NBJ(1),nr,DXINU,GL,GU,RG(ng),XG,ERROR,*9999)
C ***       Calculate stress tensor TG & e.values & vectors
            IF(ITYP1(nr,nx).EQ.4) THEN
              IF(NW.EQ.5) THEN
                DO nb=1,NBFT
                  IF(nb.EQ.NBH(NH_LOC(1,nx),nc).or.
     '               nb.EQ.NBH(NH_LOC(2,nx),nc)
     '              .or.nb.EQ.NBH(3,nc)) THEN
                    DO ns=1,NST(nb)+NAT(nb)
                      DO ni=1,2
                        PPG(nb,ni,ns)=PGX(nb,ni,ns,DXINU,PG(1,1,ng,nb))
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
C ***           Calculate metric tensors wrt nu coords
                CALL ZEAZ45(NBJ,ng,nr,AXL,AZ,AZL,PG,PPG,
     '            XE,XG,ZE,ZG,ERROR,*9999)
                CALL AZTG45(NW,AXL,AZL,CG(1,ng),EG,TG,ERROR,*9999)
              ELSE IF(NW.eq.11.or.NW.EQ.12) THEN
                CALL ZEZG(1,NBH,ng,NHE,nx,DXINU,PG,ZE,ZG,ERROR,*9999)
C cpb 13/5/96 Old way
C                dX_dNu(1,1)=XG(1,NU1(1))*DXINU(1,1)+XG(1,NU1(2))
C     '            *DXINU(2,1)
C                dX_dNu(2,1)=XG(2,NU1(1))*DXINU(1,1)+XG(2,NU1(2))
C     '            *DXINU(2,1)
C                dX_dNu(1,2)=XG(1,NU1(1))*DXINU(1,2)+XG(1,NU1(2))
C     '            *DXINU(2,2)
C                dX_dNu(2,2)=XG(2,NU1(1))*DXINU(1,2)+XG(2,NU1(2))
C     '            *DXINU(2,2)
C                CALL CGS11(nr,NW,nx,CG(1,ng),
C     '            EG,TG,dX_dNu,ZG,ERROR,*9999)
                CALL MAT_VEC_NG(2,nr,RT(1,1),RT(1,2),RT(1,3),XG,
     '            ERROR,*9999)
                CALL CGS11(nr,nw,nx,CG(1,ng),EGNU,RT,TGNU,ZG,
     '            ERROR,*9999)
C cpb 23/5/97 not needed now
CC CPB 13/5/95 Rotate the stress tensor to get back into reference
CC coordinates
C                DO i=1,2
C                  DO j=1,2
CC cbp 22/5/97 Treat fibre/no fibre case the same.
CC                    IF(NJ_LOC(NJL_FIBR,0,nr).NE.0) THEN
C                    TG(i,j)=0.0D0
C                    EG(i,j)=0.0D0
C                    DO k=1,2
C                      DO l=1,2
C                        TG(i,j)=TG(i,j)+RT(k,i)*TGNU(k,l)*RT(l,j)
C                        EG(i,j)=EG(i,j)+RT(k,i)*EGNU(k,l)*RT(l,j)
C                      ENDDO !l
C                    ENDDO !k
CC                    ELSE
CC                      TG(i,j)=TGNU(i,j)
CC                      EG(i,j)=EGNU(i,j)
CC                    ENDIF
C                  ENDDO !j
C                ENDDO !i
C              ENDIF
C              CALL EVALUE(2,TG,EVAL,ERROR,*9999)
C              CALL EVECTR(2,TG,EVAL(1),EVEC,ERROR,*9999)
              ENDIF
              CALL EVALUE(2,TGNU,EVAL,ERROR,*9999)
              CALL EVECTR(2,TGNU,EVAL(1),EVEC,ERROR,*9999)
              PHI=DATAN2(EVEC(2),EVEC(1))

            ELSE IF(ITYP1(nr,nx).EQ.5) THEN
              CALL ZEZG(1,NBH,ng,NHE,nx,DXINU,PG,ZE,ZG,ERROR,*9999)
              CALL ZGMG(NBH(NH_LOC(1,nx),nc),nr,AZ,AZL,AZU,ZG,
     '          ERROR,*9999)
              nb=NBH(NH_LOC(1,nx),nc)
              CALL ZETX50('Fibre','Cauchy',STRESSTYPE,
     '          IBT,IDO,INP,NAN,NBH(1,nc),NBJ,ng,NHE,NPNE,nr,ne,nx,
     '          CE,CG,CP,FEXT(1,ng),PG,PHII,EVAL,RG(ng),RG2D,
     '          RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XIG(1,ng,nb),YG(1,ng),
     '          ZE,ZG,ERROR,*9999)
              PHI=PHII(1) !first principal angle
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' PHI = '',E12.4,'' degrees'')')
     '          PHI*180.0d0/PI
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

C           Rotate fibre axis into principal axis
            TEMP_FIBRE=.FALSE.
            IF(NJ_LOC(NJL_FIBR,0,nr).GE.1) THEN !fibres defined
              XG(NJ_LOC(NJL_FIBR,1,nr),1)=XG(NJ_LOC(NJL_FIBR,1,nr),1)+
     '          PHI
            ELSE
              NJ_LOC00_temp=NJ_LOC(0,0,0)
              NJ_LOC(0,0,0)=NJ_LOC(0,0,0)+1
              CALL ASSERT(NJ_LOC(0,0,0).LE.NJM,'>>Increase NJM',
     '          ERROR,*9999)
              NJ_LOC(NJL_FIBR,0,nr)=1
              NJ_LOC(NJL_FIBR,1,nr)=NJ_LOC(0,0,0)
              XG(NJ_LOC(NJL_FIBR,1,nr),1)=PHI
              TEMP_FIBRE=.TRUE.
            ENDIF
C ***       Calc derivs of Xi wrt coords now aligned with princ.axes
C ***       (DXINU), & metrics (GL,GU) wrt princ.axes
            CALL XEXG(NBJ,ng,nr,PG,XE,XG,ERROR,*9999)
            CALL XGMG(1,0,NBJ(1),nr,DXINU,GL,GU,G,XG,ERROR,*9999)
            IF(TEMP_FIBRE) THEN
C             if fibres have been temporarily defined undefine them
              XG(NJ_LOC(NJL_FIBR,1,nr),1)=0.0d0
              NJ_LOC(NJL_FIBR,1,nr)=0
              NJ_LOC(NJL_FIBR,0,nr)=0
              NJ_LOC(0,0,0)=NJ_LOC00_temp
              TEMP_FIBRE=.FALSE.
            ENDIF
CC ***       Calc derivs of Xi wrt nu coords now aligned with princ.axes
C            IF(NJ_LOC(NJL_FIBR,1,nr).GT.0) THEN !fibres defined
C              ETA=XG(NJ_LOC(NJL_FIBR,1,nr),1)
C            ELSE
C              ETA=0.d0
C            ENDIF
C            IF(JTYP12.LE.1) THEN
C              ETA1=ETA+PHI
C              ETA2=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA1
C            ELSE IF(JTYP12.EQ.2) THEN
C              ETA2=ETA-PHI
C              ETA1=DACOS(GL(1,2)/DSQRT(GL(1,1)*GL(2,2)))-ETA2
C            ENDIF
C            IF(DOP) THEN
C              WRITE(OP_STRING,'('' ETA1= '',F6.1,'' degrees'')')
C     '           ETA1*180.d0/PI
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C              WRITE(OP_STRING,'('' ETA2= '',F6.1,'' degrees'')')
C     '           ETA2*180.d0/PI
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
C            C1=DCOS(ETA1)
C            C2=DCOS(ETA2)
C            RK=DSQRT(GL(1,1)/GL(2,2))*C1/C2
C            G=GL(1,1)*GL(2,2)-GL(1,2)**2
C            DXINU(1,1)=( GL(2,2)*DSQRT(GL(1,1))*C1
C     '                  -GL(1,2)*DSQRT(GL(2,2))*C2)/G
C            DXINU(2,1)=(-GL(1,2)*DSQRT(GL(1,1))*C1
C     '                  +GL(1,1)*DSQRT(GL(2,2))*C2)/G
C            IF(DABS(C2).GT.1.E-5) THEN
C              RK=DSQRT(GL(1,1)/GL(2,2))*C1/C2
C              DXINU(1,2)=-1.d0/DSQRT(GL(1,1)+RK*(RK*GL(2,2)-
C     '          2.d0*GL(1,2)))
C              DXINU(2,2)=-RK*DXINU(1,2)
C            ELSE
C              DENOM=DSQRT(G*GL(1,1))
C              DXINU(1,2)=-GL(1,2)/DENOM
C              DXINU(2,2)= GL(1,1)/DENOM
C            ENDIF
C            IF(DOP) THEN
C              WRITE(OP_STRING,'('' DXINU: '',4E12.4)')
C     '         ((DXINU(I,J),J=1,2),I=1,2)
C              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C            ENDIF
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              DO ni=1,2
                DXDNU(nj,ni)=DXINU(1,ni)*XG(nj,NU1(1))
     '                      +DXINU(2,ni)*XG(nj,NU1(2))
              ENDDO !ni
            ENDDO !nj
            IF(DOP) THEN
              WRITE(OP_STRING,'('' DXDNU: '',6D12.4)')
     '          ((DXDNU(nj,i),i=1,2),nj=1,NJ_LOC(NJL_GEOM,0,nr))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF

C ***       Princ. stress vectors are drawn at Gauss point Xi coords XIG
C           XIP1=XIG(1,ng,NBH(NH_LOC(1,nx),nc))
C           XIP2=XIG(2,ng,NBH(NH_LOC(1,nx),nc))

            DO j=1,2
              IF(DOP) THEN
                WRITE(OP_STRING,'('' J='',I2)') j
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(PRSTMAX.GT.0.d0) THEN
                DNU=EVAL(j)/PRSTMAX*0.5D0*SCALE
              ELSE
                DNU=EVAL(j)*0.5D0*SCALE
              ENDIF
              DO nj=1,NJT
                DX(nj)=DNU*DXDNU(nj,j)
              ENDDO
              IF(DOP) THEN
C               WRITE(IO4,'('' XIP1='',E12.4,'' XIP2='',E12.4)') XIP1,
C    '            XIP2
                WRITE(OP_STRING,'('' DNU='',E12.4,'
     '            //''' DX(nj)='',3E12.4)') DNU,(DX(nj),nj=1,NJT)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
C             XI(1)=XIP1-DX1/2.d0
C             XI(2)=XIP2-DX2/2.d0
C             DO nj=1,NJT
C               IF(ITYP6(nr,nx).EQ.1) THEN
C                 nb=NBJ(nj)
C                 X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C    '              nb,1,XI,XE(1,nj))
C               ELSE IF(ITYP6(nr,nx).EQ.2) THEN
C                 nb=NBH(nj,nc)
C                 X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C    '              nb,1,XI,ZE(1,nj))
C               ENDIF
C             ENDDO
              DO nj=1,NJT
                IF(ITYP6(nr,nx).EQ.1) THEN          !linear problem
                  X(nj)=XG(nj,1)-DX(nj)/2.0d0
                ELSE IF(ITYP6(nr,nx).EQ.2) THEN     !nonlinear problem
                  X(nj)=ZG(nj,1)-DX(nj)/2.0d0
                ENDIF
              ENDDO

              CALL XZ(ITYP10(1),X,XY(1,1))
C             XI(1)=XIP1+DX1/2.d0
C             XI(2)=XIP2+DX2/2.d0
C             DO nj=1,NJT
C               IF(ITYP6(nr,nx).EQ.1) THEN
C                 nb=NBJ(nj)
C                 X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C    '              nb,1,XI,XE(1,nj))
C               ELSE IF(ITYP6(nr,nx).EQ.2) THEN
C                 nb=NBH(nj,nc)
C                 X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
C    '              nb,1,XI,ZE(1,nj))
C               ENDIF
C             ENDDO
              DO nj=1,NJT
                IF(ITYP6(nr,nx).EQ.1) THEN          !linear problem
                  X(nj)=XG(nj,1)+DX(nj)/2.0d0
                ELSE IF(ITYP6(nr,nx).EQ.2) THEN     !nonlinear problem
                  X(nj)=ZG(nj,1)+DX(nj)/2.0d0
                ENDIF
              ENDDO
              CALL XZ(ITYP10(1),X,XY(1,2))
              CALL POLYLINE(INDEX,IW,2,XY,ERROR,*9999)
            ENDDO
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('STRESS1')
      RETURN
 9999 IF(TEMP_FIBRE) THEN
C       If fibres have been temporarily defined undefine them
        XG(NJ_LOC(NJL_FIBR,1,nr),1)=0.0d0
        NJ_LOC(NJL_FIBR,1,nr)=0
        NJ_LOC(0,0,0)=NJ_LOC00_temp
        TEMP_FIBRE=.FALSE.
      ENDIF
      CALL ERRORS('STRESS1',ERROR)
      CALL EXITS('STRESS1')
      RETURN 1
      END
