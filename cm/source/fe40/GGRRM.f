      SUBROUTINE GGRRM(NBJ,GGM,GU,RRM,X3G,XG,ZX,ERROR,*)

C#### Subroutine: GGRRM
C###  Description:

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER NBJ(NJM)
      REAL*8 GGM(2,*),GU(3,*),RRM(2,*),X3G(4,*),XG(NJM,NUM),ZX(10,*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ia,ib,id,il,im,NITB,nr,NU1,NU2,NU3,NV1(3),NV2(3,3)
      REAL*8 BL(3,3),BM(3,3),CHTOFF(3,3,3),DBM(3,3,3),SUM,SUM1,SUM2

      DATA NV1/2,4,7/
      DATA NV2/3,6,9,6,5,10,9,10,8/

      CALL ENTERS('GGRRM',*9999)
      nr=1 !temporary
      CALL TOFFEL(ITYP10(nr),NBJ(1),nr,CHTOFF,DBM,GU,XG,X3G,
     '  .TRUE.,ERROR,*9999)
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
      DO 60 ib=1,NITB
        SUM1=0.d0
        SUM2=0.d0
        NU2=NV1(ib)
        DO 50 il=1,NITB
          NU3=NV2(il,ib)
          SUM2=SUM2+GU(ia,il)*ZX(NU3,3)
          NU1=NV1(il)
          SUM1=SUM1+GU(ia,il)*(0.50d0*(ZX(NU1,ib)+ZX(NU2,il))-BL(il,ib)
     '             *ZX(1,3))
          DO 40 im=1,NITB
            SUM1=SUM1-(CHTOFF(im,ib,il)+CHTOFF(im,il,ib))*ZX(1,im)
     '               *0.5d0*GU(ia,il)
            SUM2=SUM2+GU(ia,il)*( -CHTOFF(im,il,ib)*ZX(NV1(im),3)+
     '               DBM(im,ib,il)
     '               *ZX(1,im)+BM(im,ib)*ZX(NU1,im)+BM(im,il)
     '               *ZX(NU2,im)-BM(im,il)*BL(im,ib)*ZX(1,3))
            DO 30 id=1,NITB
              SUM2=SUM2+GU(ia,il)*(-CHTOFF(im,id,il)
     '                 *BM(id,ib)-CHTOFF(im,id,ib)*BM(id,il))*ZX(1,im)
 30         CONTINUE
 40       CONTINUE
 50     CONTINUE
        GGM(ia,ib)=SUM1
        RRM(ia,ib)=SUM2
C       WRITE(3,590)ia,ib,SUM1,SUM2
C       WRITE(4,590)ia,ib,SUM1,SUM2
C590    FORMAT(' ia=',I2,' ib=',I2,' ggm(ia,ib)= ',E10.3,
C    '         ' rrm(ia,ib)= ',E10.3)
 60   CONTINUE

      CALL EXITS('GGRRM')
      RETURN
 9999 CALL ERRORS('GGRRM',ERROR)
      CALL EXITS('GGRRM')
      RETURN 1
      END


