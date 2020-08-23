      SUBROUTINE MIDMN(NBJ,GGM,GL,GU,POIS,RRM,SMX,SNX,THIC,
     '  YMOD,ERROR,*)

C#### Subroutine: MIDMN
C###  Description:
C###    MIDMN computes the midsurface stress resultants - SMX and SNX.
C###    These stresses are physical components. (The units are
C###    FORCE/UNIT LENGTH).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM)
      REAL*8 GGM(2,*),GL(3,*),GU(3,*),POIS,RRM(2,*),SMX(2,*),SNX(2,*),
     '  THIC,YMOD
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ia,ib,il,im,j,NITB
      REAL*8 D,GGL(2,2),RRL(2,2),SK,SM(2,2),SN(2,2),SUM1,SUM2

      CALL ENTERS('MIDMN',*9999)
      DO ia=1,2
        DO ib=1,2
          SUM1=0.0d0
          SUM2=0.0d0
          DO il=1,2
            SUM1=SUM1+GL(ia,il)*GGM(il,ib)
            SUM2=SUM2+GL(ia,il)*RRM(il,ib)
          ENDDO
          GGL(ia,ib)=SUM1
          RRL(ia,ib)=SUM2
        ENDDO
      ENDDO

C     DO ia=1,2
C       DO ib=1,2
C         RKL(ia,ib)=0.50d0*(RRL(ia,ib)+RRL(ib,ia))
C       ENDDO
C     ENDDO

      D=YMOD*THIC/(1.0d0-POIS**2)
      SK=YMOD*THIC**3/(12.0d0*(1.0d0-POIS**2))

      NITB=NIT(NBJ(1))
      DO ia=1,NITB
      DO ib=1,NITB
          SUM1=0.0d0
          SUM2=0.0d0
          DO il=1,NITB
            DO im=1,NITB
              SUM1=SUM1+((1.0d0-POIS)*GU(ia,il)*GU(ib,im)+POIS*GU(ia,ib)
     '          *GU(il,im))*GGL(il,im)
C             SUM2=SUM2+((1.0d0-POIS)*GU(ia,il)*GU(ib,im)+POIS*GU(ia,ib)
C    '          *GU(il,im))*RKL(il,im)
              SUM2=SUM2+((1.0d0-POIS)*GU(ia,il)*GU(ib,im)+POIS*GU(ia,ib)
     '          *GU(il,im))*RRL(il,im)
            ENDDO
          ENDDO
          SN(ia,ib)=D*SUM1
          SM(ia,ib)=SK*SUM2
        ENDDO
      ENDDO
C
C  Transform into the physical components (stress resultants).
C
      DO i=1,NITB
        DO j=1,NITB
          SMX(i,j)=SM(i,j)*DSQRT(GL(i,i)*GL(j,j))
          SNX(i,j)=SN(i,j)*DSQRT(GL(i,i)*GL(j,j))
        ENDDO
      ENDDO

      CALL EXITS('MIDMN')
      RETURN
 9999 CALL ERRORS('MIDMN',ERROR)
      CALL EXITS('MIDMN')
      RETURN 1
      END


