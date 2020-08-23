      SUBROUTINE SSGGX(NBJ,GGM,GGPX,GL,GU,POIS,RRM,SSPX,
     '  X3G,XG,XI3,YMOD,ERROR,*)

C#### Subroutine: SSGGX
C###  Description:
C###    SSGGX calculates shell quantities.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER NBJ(NJM)
      REAL*8 GGM(2,*),GGPX(3,*),GL(3,*),GU(3,*),POIS,RRM(2,*),SSPX(3,*),
     '  X3G(4,*),XG(NJM,NUM),XI3,YMOD
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ia,ib,il,j,NITB,nr
      REAL*8 BL(3,3),BM(3,3),CHTOFF(3,3,3),DBM(3,3,3),GGX(3,3),
     '  RGU(3),SSX(3,3),SUM,SUM1,SUM2,SUM3,SUM4

      CALL ENTERS('SSGGX',*9999)
      nr=1 !temporary
      CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM,GU,XG,X3G,.TRUE.,
     '  ERROR,*9999)
      DO ia=1,2
        DO ib=1,2
          BL(ia,ib)=CHTOFF(3,ia,ib)
          SUM=0.d0
          DO il=1,2
            SUM=SUM+GU(ia,il)*CHTOFF(3,il,ib)
          ENDDO
          BM(ia,ib)=SUM
        ENDDO
      ENDDO

      NITB=NIT(NBJ(1))
      DO 60 ia=1,NITB
        GGX(ia,3)=0.d0
        GGX(3,ia)=0.d0
        SSX(ia,3)=0.d0
        SSX(3,ia)=0.d0
      DO 60  ib=1,NITB
        SUM1=GGM(ia,ib)-XI3*RRM(ia,ib)
        SUM2=GGM(ia,ib)-XI3*RRM(ia,ib)
        SUM3=0.d0
        DO 50 il=1,NITB
          SUM1=SUM1+2.d0*XI3*BM(ia,il)*GGM(il,ib)
          SUM2=SUM2+2.d0*XI3*BM(ia,il)*GGM(il,ib)
          IF(ia.EQ.ib) THEN
            SUM3=SUM3+GGM(il,il)-XI3*(RRM(il,il)-2.d0*BM(il,1)*GGM(1,il)
     '               -2.d0*BM(il,2)*GGM(2,il))
          ENDIF
 50     CONTINUE
        SUM3=SUM3*POIS/(1.d0-POIS)
        GGX(ia,ib)=SUM1
        SSX(ia,ib)=YMOD*(SUM2+SUM3)/(1.d0+POIS)
 60   CONTINUE
      SUM4=0.d0
      DO 70 ia=1,NITB
        SUM4=SUM4+GGX(ia,ia)
 70   CONTINUE
      GGX(3,3)=-SUM4*POIS/(1.d0-POIS)
      SSX(3,3)=0.d0
C
C  Transform into the (right-) physical components (strain and stress).
C
      DO 80 ia=1,NITB
        RGU(ia)=DSQRT(GL(ia,ia))*(1.d0-BL(ia,ia)*XI3/GL(ia,ia))
 80   CONTINUE
      RGU(3)=1.d0
C       WRITE(4,100)(RGU(i),i=1,3)
C100    FORMAT(/,'RGU: ',3E10.3,/)
      DO 90 i=1,NITB+1
      DO 90 j=1,NITB+1
        GGPX(i,j)=RGU(i)/RGU(j)*GGX(i,j)
        SSPX(i,j)=RGU(i)/RGU(j)*SSX(i,j)
 90   CONTINUE

      CALL EXITS('SSGGX')
      RETURN
 9999 CALL ERRORS('SSGGX',ERROR)
      CALL EXITS('SSGGX')
      RETURN 1
      END


