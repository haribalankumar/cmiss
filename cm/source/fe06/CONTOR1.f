      SUBROUTINE CONTOR1(INDEX,IBT,IDO,INP,
     '  ITMAX,IW,NBJ,nh,nj,NROOTS,
     '  ROOTS,TOL,XE,XICONT,XISTEP,XLM,ZE,ZVAL,ERROR,*)

C#### Subroutine: CONTOR1
C###  Description:
C###    CONTOR1 is a fast version of CONTOR for bilinear field.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ITMAX,IW,NBJ(NJM),nh,nj,NROOTS
!     Local Variables
      REAL*8 ROOTS(2,12),TOL,XE(NSM,NJM),XICONT,XLM(3),ZE(NSM,NHM),ZVAL
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ISIGN,ISTEP,it,j,MSTEP,nb,nroo,nroot,nstep,NSTEPS
      REAL*8 DENOM,DXI(2),DZXI(2),PL111,PL112,PL121,PL122,PL211,PL212,
     '  PL221,PL222,PXI,XI(4),XI1,XI2,XINEW(4),XIOLD(4),XISTEP,
     '  XL(3,500),Z,ZE01,ZE02,ZE03,ZE04,ZTOL
      LOGICAL LROOT(12)

      CALL ENTERS('CONTOR1',*9999)
      ZTOL=TOL*(1.d0+DABS(ZVAL))
      ZE01=ZE( 1,nh)
      ZE02=ZE( 2,nh)
      ZE03=ZE( 3,nh)
      ZE04=ZE( 4,nh)

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
            PL111=1.d0-XI1
            PL211=XI1
            PL121=-1.d0
            PL221= 1.d0
            PL112=1.d0-XI2
            PL212=XI2
            PL122=-1.d0
            PL222= 1.d0
            Z      =PL111*PL112*ZE01 + PL211*PL112*ZE02
     '             +PL111*PL212*ZE03 + PL211*PL212*ZE04
            DZXI(1)=PL121*PL112*ZE01 + PL221*PL112*ZE02
     '             +PL121*PL212*ZE03 + PL221*PL212*ZE04
            DZXI(2)=PL111*PL122*ZE01 + PL211*PL122*ZE02
     '             +PL111*PL222*ZE03 + PL211*PL222*ZE04

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

      CALL EXITS('CONTOR1')
      RETURN
 9999 CALL ERRORS('CONTOR1',ERROR)
      CALL EXITS('CONTOR1')
      RETURN 1
      END


