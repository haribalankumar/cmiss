      SUBROUTINE CROOTS1(ITMAX,nh,NROOTS,ROOTS,TOL,ZE,ZVAL,ERROR,*)

C#### Subroutine: CROOTS1
C###  Description:
C###    CROOTS1 is a fast version of CROOTS for bilinear field.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ITMAX,nh,NROOTS
      REAL*8 ROOTS(2,12),TOL,ZE(NSM,NHM),ZVAL
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER it,n,ni1,NI2,nixi,nroot
      REAL*8 DXI,DXIMAX,DZXI(2),PL111,PL112,PL121,PL122,PL211,PL212,
     '  PL221,PL222,XI(2),XI1,XI2,Z,ZTOL,ZE01,ZE02,ZE03,ZE04

      CALL ENTERS('CROOTS1',*9999)
      NROOTS=0
      DXIMAX=0.25D0
      ZTOL=TOL*(1.d0+DABS(ZVAL))
      ZE01=ZE( 1,nh)
      ZE02=ZE( 2,nh)
      ZE03=ZE( 3,nh)
      ZE04=ZE( 4,nh)
      DO ni1=1,2
        NI2=MOD(ni1,2)+1
        DO nixi=0,1
          XI(ni1)=DBLE(nixi)
          DO n=1,4
            XI(NI2)=DBLE(2*n-1)/8.d0
            DO it=1,ITMAX
              XI1=XI(1)
              XI2=XI(2)
              PL111=1.d0-XI1
              PL211=XI1
              PL121=-1.d0
              PL221= 1.d0
              PL112=1.d0-XI2
              PL212=XI2
              PL122=-1.d0
              PL222= 1.d0
              Z      =PL111*PL112*ZE01 + PL211*PL112*ZE02
     '               +PL111*PL212*ZE03 + PL211*PL212*ZE04
              DZXI(1)=PL121*PL112*ZE01 + PL221*PL112*ZE02
     '               +PL121*PL212*ZE03 + PL221*PL212*ZE04
              DZXI(2)=PL111*PL122*ZE01 + PL211*PL122*ZE02
     '               +PL111*PL222*ZE03 + PL211*PL222*ZE04
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

      CALL EXITS('CROOTS1')
      RETURN
 9999 CALL ERRORS('CROOTS1',ERROR)
      CALL EXITS('CROOTS1')
      RETURN 1
      END


