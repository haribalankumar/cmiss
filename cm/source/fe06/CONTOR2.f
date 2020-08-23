      SUBROUTINE CONTOR2(INDEX,IBT,IDO,INP,
     '  ITMAX,IW,NBJ,nh,nj,NROOTS,
     '  ROOTS,TOL,XE,XICONT,XISTEP,XLM,ZE,ZVAL,ERROR,*)

C#### Subroutine: CONTOR2
C###  Description:
C###    CONTOR2 is a fast version of CONTOR for bicubic Hermite field.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ITMAX,IW,NBJ(NJM),nh,nj,NROOTS
      REAL*8 ROOTS(2,12),TOL,XE(NSM,NJM),XICONT,XLM(3),ZE(NSM,NHM),ZVAL
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ISIGN,ISTEP,it,j,MSTEP,nb,nroo,nroot,nstep,NSTEPS
      REAL*8 DENOM,DXI(2),DZXI(2),PH1111,PH1112,PH1121,PH1122,PH1211,
     '  PH1212,PH1221,PH1222,PH2111,PH2112,PH2121,PH2122,PH2211,PH2212,
     '  PH2221,PH2222,PXI,XI(4),XI1,XI2,XINEW(4),XIOLD(4),XISTEP,
     '  XL(3,500),Z,ZE01,ZE02,ZE03,ZE04,ZE05,ZE06,ZE07,ZE08,
     '  ZE09,ZE10,ZE11,ZE12,ZE13,ZE14,ZE15,ZE16,ZTOL
      LOGICAL LROOT(12)

      CALL ENTERS('CONTOR2',*9999)
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

C *** LROOT is true if root lies between 0 & 1

      DO nroot=1,NROOTS
        IF ((DABS(ROOTS(1,nroot)).LE.TOL)
     '  .OR.(DABS(ROOTS(1,nroot)-1.d0).LE.TOL)
     '  .OR.(DABS(ROOTS(2,nroot)).LE.TOL)
     '  .OR.(DABS(ROOTS(2,nroot)-1.d0).LE.TOL)) THEN
          LROOT(nroot)=.TRUE.
        ELSE
          LROOT(nroot)=.FALSE.
        ENDIF
      ENDDO
      NSTEPS=4*IDNINT(1.d0/XISTEP)
      CALL ASSERT(NSTEPS.LT.500,'Too many steps in contour',ERROR,*9999)
      ISTEP=0
      DO nroot=1,NROOTS
        IF(LROOT(nroot)) THEN
          LROOT(nroot)=.FALSE.
          XI(1)=ROOTS(1,nroot)
          XI(2)=ROOTS(2,nroot)
          XI(3)=XICONT
          ISTEP=1
          IF(IW.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
            XL(1,ISTEP)=XI(1)
            XL(2,ISTEP)=XI(2)
          ELSE
            DO nj=1,NJT
              nb=NBJ(nj)
              XL(nj,ISTEP)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '          nb,1,XI,XE(1,nj))
            ENDDO
          ENDIF
          ISIGN=1
          XIOLD(1)=-1.d0
          XIOLD(2)=-1.d0
          DO nstep=1,NSTEPS
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
     '             +PH1111*PH1212*ZE03 + PH1211*PH1212*ZE04
     '             +PH2111*PH1112*ZE05 + PH2211*PH1112*ZE06
     '             +PH2111*PH1212*ZE07 + PH2211*PH1212*ZE08
     '             +PH1111*PH2112*ZE09 + PH1211*PH2112*ZE10
     '             +PH1111*PH2212*ZE11 + PH1211*PH2212*ZE12
     '             +PH2111*PH2112*ZE13 + PH2211*PH2112*ZE14
     '             +PH2111*PH2212*ZE15 + PH2211*PH2212*ZE16
            DZXI(1)=PH1121*PH1112*ZE01 + PH1221*PH1112*ZE02
     '             +PH1121*PH1212*ZE03 + PH1221*PH1212*ZE04
     '             +PH2121*PH1112*ZE05 + PH2221*PH1112*ZE06
     '             +PH2121*PH1212*ZE07 + PH2221*PH1212*ZE08
     '             +PH1121*PH2112*ZE09 + PH1221*PH2112*ZE10
     '             +PH1121*PH2212*ZE11 + PH1221*PH2212*ZE12
     '             +PH2121*PH2112*ZE13 + PH2221*PH2112*ZE14
     '             +PH2121*PH2212*ZE15 + PH2221*PH2212*ZE16
            DZXI(2)=PH1111*PH1122*ZE01 + PH1211*PH1122*ZE02
     '             +PH1111*PH1222*ZE03 + PH1211*PH1222*ZE04
     '             +PH2111*PH1122*ZE05 + PH2211*PH1122*ZE06
     '             +PH2111*PH1222*ZE07 + PH2211*PH1222*ZE08
     '             +PH1111*PH2122*ZE09 + PH1211*PH2122*ZE10
     '             +PH1111*PH2222*ZE11 + PH1211*PH2222*ZE12
     '             +PH2111*PH2122*ZE13 + PH2211*PH2122*ZE14
     '             +PH2111*PH2222*ZE15 + PH2211*PH2222*ZE16

            IF(DABS(Z-ZVAL).GT.ZTOL) THEN
              WRITE(OP_STRING,
     '          '('' !!!!!!Tolerance of contour exceeded at XI='','
     '          //'2(G12.4,4X)/''      Z='',G12.4,'
     '          //'/'' Z-ZVAL='',G12.4,/''   ZTOL='',G12.4)')
     '          XI(1),XI(2),Z,Z-ZVAL,ZTOL
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              GO TO 6
            ENDIF

C ***       Find increment normal to gradient & update Xi

            DENOM=DSQRT(DZXI(1)**2+DZXI(2)**2)*ISIGN
            IF(DABS(DENOM).GT.1.d-6) THEN
              DXI(1)=XISTEP*DZXI(2)/DENOM
              DXI(2)=XISTEP*DZXI(1)/DENOM
            ELSE
              DXI(1)=0.d0
              DXI(2)=0.d0
            ENDIF
            XINEW(1)=XI(1)+DXI(1)
            XINEW(2)=XI(2)-DXI(2)

C ***       Step in gradient direction until converged

            DO it=1,ITMAX
              XI1=XINEW(1)
              XI2=XINEW(2)
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

              DENOM=DZXI(1)**2+DZXI(2)**2
              IF(DABS(DENOM).GT.1.d-6) THEN
                DXI(1)=(Z-ZVAL)*DZXI(1)/DENOM
                DXI(2)=(Z-ZVAL)*DZXI(2)/DENOM
              ELSE
                DXI(1)=0.d0
                DXI(2)=0.d0
              ENDIF
              XINEW(1)=XINEW(1)-DXI(1)
              XINEW(2)=XINEW(2)-DXI(2)
              IF(((DXI(1)**2+DXI(2)**2).LE.TOL**2).AND.
     '          (DABS(Z-ZVAL).LE.ZTOL)) GO TO 3
            ENDDO

C ***       Avoid contour turning back on itself

 3          IF(((XINEW(1)-XIOLD(1))**2+
     '          (XINEW(2)-XIOLD(2))**2).LE.TOL**2) GO TO 6
            XIOLD(1)=XI(1)
            XIOLD(2)=XI(2)
            XI(1)=XINEW(1)
            XI(2)=XINEW(2)

            IF((XI(1).LE.0.d0).OR.(XI(1).GE.1.d0).OR.
     '         (XI(2).LE.0.d0).OR.(XI(2).GE.1.d0)) THEN

C ***         Reverse direction if 1st step moves outside element

              IF((nstep.EQ.1).AND.(ISIGN.EQ.1)) THEN
                ISIGN=-1
                XI(1)=ROOTS(1,nroot)
                XI(2)=ROOTS(2,nroot)
                GO TO 5
              ENDIF

C ***         When reach bdry join to that root (if it exists) with a
C ***         linear approx & then eliminate that root

              ISTEP=ISTEP+1
              DO nroo=nroot,NROOTS
                IF(LROOT(nroo)) THEN
                  IF(((XI(1)-ROOTS(1,nroo))**2+
     '                (XI(2)-ROOTS(2,nroo))**2).LE.4.d0*XISTEP**2) THEN
                    LROOT(nroo)=.FALSE.
                    XI(1)=ROOTS(1,nroo)
                    XI(2)=ROOTS(2,nroo)
                    IF(IW.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
                      XL(1,ISTEP)=XI(1)
                      XL(2,ISTEP)=XI(2)
                    ELSE
                      DO nj=1,NJT
                        nb=NBJ(nj)
                        XL(nj,ISTEP)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                    INP(1,1,nb),nb,1,XI,XE(1,nj))
                      ENDDO
                    ENDIF
                    GO TO 550
                  ENDIF
                ENDIF
              ENDDO
              IF(XI(1).LE.0.d0) THEN
                DENOM=XI(1)-XIOLD(1)
                IF(DABS(DENOM).GT.1.d-6) THEN
                  XI(2)=(XI(1)*XIOLD(2)-XI(2)*XIOLD(1))/DENOM
                ENDIF
                XI(1)=0.d0
              ELSE IF(XI(2).LE.0.d0) THEN
                DENOM=XI(2)-XIOLD(2)
                IF(DABS(DENOM).GT.1.d-6) THEN
                  XI(1)=(XI(2)*XIOLD(1)-XI(1)*XIOLD(2))/DENOM
                ENDIF
                XI(2)=0.d0
              ELSE IF(XI(1).GE.1.d0) THEN
                DENOM=XI(1)-XIOLD(1)
                IF(DABS(DENOM).GT.1.d-6) THEN
                  XI(2)=(XI(1)*XIOLD(2)-XI(2)*XIOLD(1)+XI(2)-XIOLD(2))
     '               /DENOM
                ENDIF
                XI(1)=1.d0
              ELSE IF(XI(2).GE.1.d0) THEN
                DENOM=XI(2)-XIOLD(2)
                IF(DABS(DENOM).GT.1.d-6) THEN
                  XI(1)=(XI(2)*XIOLD(1)-XI(1)*XIOLD(2)+XI(1)-XIOLD(1))
     '              /DENOM
                ENDIF
                XI(2)=1.d0
              ENDIF
              IF(IW.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
                XL(1,ISTEP)=XI(1)
                XL(2,ISTEP)=XI(2)
              ELSE
                DO nj=1,NJT
                  nb=NBJ(nj)
                  XL(nj,ISTEP)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,1,XI,XE(1,nj))
                ENDDO
              ENDIF
              GO TO 550
            ENDIF
            ISTEP=ISTEP+1
            IF(IW.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
              XL(1,ISTEP)=XI(1)
              XL(2,ISTEP)=XI(2)
            ELSE
              DO nj=1,NJT
                nb=NBJ(nj)
                XL(nj,ISTEP)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI,XE(1,nj))
              ENDDO
            ENDIF
 5          CONTINUE
          ENDDO

 550      CALL POLYLINE(INDEX,IW,ISTEP,XL,ERROR,*9999)

        ENDIF
 6      CONTINUE
      ENDDO
      IF(ISTEP.GT.0) THEN
        MSTEP=ISTEP/2
        DO j=1,NJT
          XLM(j)=XL(j,MSTEP)
        ENDDO
      ENDIF

      CALL EXITS('CONTOR2')
      RETURN
 9999 CALL ERRORS('CONTOR2',ERROR)
      CALL EXITS('CONTOR2')
      RETURN 1
      END


