      SUBROUTINE CROOTS2(ITMAX,nh,NROOTS,ROOTS,TOL,ZE,ZVAL,ERROR,*)

C#### Subroutine: CROOTS2
C###  Description:
C###    CROOTS2 is a fast version of CROOTS for bicubic Hermite field.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ITMAX,nh,NROOTS
      REAL*8 ROOTS(2,12),TOL,ZE(NSM,NHM),ZVAL
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER it,n,ni1,NI2,nixi,nroot
      REAL*8 DXI,DXIMAX,DZXI(2),PH1111,PH1112,PH1121,PH1122,PH1211,
     '  PH1212,PH1221,PH1222,PH2111,PH2112,PH2121,PH2122,PH2211,PH2212,
     '  PH2221,PH2222,XI(2),XI1,XI2,Z,ZTOL,ZE01,ZE02,ZE03,ZE04,ZE05,
     '  ZE06,ZE07,ZE08,ZE09,ZE10,ZE11,ZE12,ZE13,ZE14,ZE15,ZE16

      CALL ENTERS('CROOTS2',*9999)
      NROOTS=0
      DXIMAX=0.25D0
      ZTOL=TOL*(1.d0+DABS(ZVAL))
      ZE01=ZE( 1,nh)
      ZE02=ZE( 2,nh)
      ZE03=ZE( 3,nh)
      ZE04=ZE( 4,nh)
      ZE05=ZE( 5,nh)
      ZE06=ZE( 6,nh)
      ZE07=ZE( 7,nh)
      ZE08=ZE( 8,nh)
      ZE09=ZE( 9,nh)
      ZE10=ZE(10,nh)
      ZE11=ZE(11,nh)
      ZE12=ZE(12,nh)
      ZE13=ZE(13,nh)
      ZE14=ZE(14,nh)
      ZE15=ZE(15,nh)
      ZE16=ZE(16,nh)
      DO ni1=1,2
        NI2=MOD(ni1,2)+1
        DO nixi=0,1
          XI(ni1)=DBLE(nixi)
          DO n=1,4
            XI(NI2)=DBLE(2*n-1)/8.d0
            DO it=1,ITMAX
              XI1=XI(1)
              XI2=XI(2)
              PH1111=1.d0+XI1*XI1*(2.d0*XI1-3.d0)
              PH1211=XI1*(XI1-1.d0)*(XI1-1.d0)
              PH2111=XI1*XI1*(3.d0-2.d0*XI1)
              PH2211=XI1*XI1*(XI1-1.d0)
              PH1121=6.d0*XI1*(XI1-1.d0)
              PH1221=(XI1-1.d0)*(3.d0*XI1-1.d0)
              PH2121=6.d0*XI1*(1.d0-XI1)
              PH2221=XI1*(3.d0*XI1-2.d0)
              PH1112=1.d0+XI2*XI2*(2.d0*XI2-3.d0)
              PH1212=XI2*(XI2-1.d0)*(XI2-1.d0)
              PH2112=XI2*XI2*(3.d0-2.d0*XI2)
              PH2212=XI2*XI2*(XI2-1.d0)
              PH1122=6.d0*XI2*(XI2-1.d0)
              PH1222=(XI2-1.d0)*(3.d0*XI2-1.d0)
              PH2122=6.d0*XI2*(1.d0-XI2)
              PH2222=XI2*(3.d0*XI2-2.d0)
              Z      =PH1111*PH1112*ZE01 + PH1211*PH1112*ZE02
     '               +PH1111*PH1212*ZE03 + PH1211*PH1212*ZE04
     '               +PH2111*PH1112*ZE05 + PH2211*PH1112*ZE06
     '               +PH2111*PH1212*ZE07 + PH2211*PH1212*ZE08
     '               +PH1111*PH2112*ZE09 + PH1211*PH2112*ZE10
     '               +PH1111*PH2212*ZE11 + PH1211*PH2212*ZE12
     '               +PH2111*PH2112*ZE13 + PH2211*PH2112*ZE14
     '               +PH2111*PH2212*ZE15 + PH2211*PH2212*ZE16
              DZXI(1)=PH1121*PH1112*ZE01 + PH1221*PH1112*ZE02
     '               +PH1121*PH1212*ZE03 + PH1221*PH1212*ZE04
     '               +PH2121*PH1112*ZE05 + PH2221*PH1112*ZE06
     '               +PH2121*PH1212*ZE07 + PH2221*PH1212*ZE08
     '               +PH1121*PH2112*ZE09 + PH1221*PH2112*ZE10
     '               +PH1121*PH2212*ZE11 + PH1221*PH2212*ZE12
     '               +PH2121*PH2112*ZE13 + PH2221*PH2112*ZE14
     '               +PH2121*PH2212*ZE15 + PH2221*PH2212*ZE16
              DZXI(2)=PH1111*PH1122*ZE01 + PH1211*PH1122*ZE02
     '               +PH1111*PH1222*ZE03 + PH1211*PH1222*ZE04
     '               +PH2111*PH1122*ZE05 + PH2211*PH1122*ZE06
     '               +PH2111*PH1222*ZE07 + PH2211*PH1222*ZE08
     '               +PH1111*PH2122*ZE09 + PH1211*PH2122*ZE10
     '               +PH1111*PH2222*ZE11 + PH1211*PH2222*ZE12
     '               +PH2111*PH2122*ZE13 + PH2211*PH2122*ZE14
     '               +PH2111*PH2222*ZE15 + PH2211*PH2222*ZE16
              IF(DZXI(NI2).EQ.0.d0) THEN
                IF((Z-ZVAL).EQ.0.d0) THEN
                  DXI=0.d0
                ELSE
                  GOTO 4
                ENDIF
              ELSE
                DXI=(ZVAL-Z)/DZXI(NI2)
              ENDIF
              IF(DABS(DXI).GT.DXIMAX) THEN
                DXI=DXIMAX*DSIGN(DXI,1.d0)
              ENDIF
              XI(NI2)=XI(NI2)+DXI
              IF((DABS(DXI).LT.TOL).AND.(DABS(Z-ZVAL).LT.ZTOL)) GOTO 2
              IF((XI(NI2).LT.0.d0).OR.(XI(NI2).GT.1.d0)) GOTO 4
            ENDDO
            GOTO 4
 2          IF((XI(NI2).LT.0.d0).OR.(XI(NI2).GT.1.d0)) GOTO 4
            DO nroot=1,NROOTS
              IF((DABS(XI(1)-ROOTS(1,nroot))
     '           +DABS(XI(2)-ROOTS(2,nroot))).LE.TOL) GOTO 4
            ENDDO
            NROOTS=NROOTS+1
            ROOTS(ni1,NROOTS)=XI(ni1)
            ROOTS(NI2,NROOTS)=XI(NI2)
 4          CONTINUE
          ENDDO
        ENDDO
      ENDDO

      CALL EXITS('CROOTS2')
      RETURN
 9999 CALL ERRORS('CROOTS2',ERROR)
      CALL EXITS('CROOTS2')
      RETURN 1
      END


