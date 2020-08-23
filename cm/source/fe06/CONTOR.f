      SUBROUTINE CONTOR(INDEX,CONTYP,IBT,IDO,INP,ICOORD,
     '  ITMAX,IW,NICONT,NAN,NBH,NBJ,nh,NHE,nj,nr,NROOTS,nx,
     '  ROOTS,TOL,PG,XE,XG,XICONT,XISTEP,XLM,
     '  ZE,ZG,ZVAL,ERROR,*)

C#### Subroutine: CONTOR
C###  Description:
C###    CONTOR plots ZVAL contours of a field defined by
C###    CFUNC(XI,Z,DZXI).  On entry the NROOTS boundary XI(ni) values
C###    where CFUNC=ZVAL are stored in ROOTS(ni,nroot).

C**** EMAP is 'GEOMETRIC' if element to device map is defined by XE or
C**** EMAP is 'USER' if the element boundries are defined in device
C**** coordinates by ELIMIT(1)=min DC(1), ELIMIT(2)=min DC(2),
C****                ELIMIT(3)=max DC(1), ELIMIT(4)=max DC(2).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),ICOORD,IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ITMAX,IW,NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),
     '  nh,NHE,NICONT,nj,nr,NROOTS,nx
      REAL*8 PG(NSM,NUM,NGM,NBM),ROOTS(2,12),TOL,
     '  XE(NSM,NJM),XG(NJM,NUM),XICONT,XISTEP,XLM(3),ZE(NSM,NHM),
     '  ZG(NHM,NUM),ZVAL
      CHARACTER CONTYP*(*),ERROR*(*)
!     Local Variables
      INTEGER i,ISIGN,ISTEP,it,j,MSTEP,nroo,nroot,nstep,NSTEPS
      REAL*8 DENOM,DXI(2),DZXI(2),XI(4),XINEW(4),XIOLD(4),XL(3,500),
     '  Z,ZTOL
      LOGICAL LROOT(12)

      CALL ENTERS('CONTOR',*9999)
      ZTOL=TOL*(1.d0+DABS(ZVAL))

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
        IF(DOP) THEN
          WRITE(OP_STRING,'(//''  nroot='',I4,7X,''ROOT='','
     '      //'2(E12.3,4X))') nroot,ROOTS(1,nroot),ROOTS(2,nroot)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
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
            CALL XIXL(IBT,IDO,INP,ISTEP,IW,NBJ,XI,XE,XL,ERROR,*9999)
          ENDIF
          ISIGN=1
          XIOLD(1)=-1.d0
          XIOLD(2)=-1.d0
          DO nstep=1,NSTEPS
            CALL CFUNC(CONTYP,IBT,ICOORD,IDO,INP,NAN,NBH,NBJ,
     '        nh,NHE,NICONT,nj,nr,nx,
     '        PG,XE,XG,XI,XICONT,ZE,ZG,Z,DZXI,ERROR,*9999)
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
              CALL CFUNC(CONTYP,IBT,ICOORD,IDO,INP,NAN,NBH,NBJ,
     '          nh,NHE,NICONT,nj,nr,nx,
     '          PG,XE,XG,XINEW,XICONT,ZE,ZG,Z,DZXI,
     '          ERROR,*9999)
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
C             IF(DOP) THEN
C               CALL CFUNC(CONTYP,IBT,ICOORD,IDO,INP,NAN,NBH,NBJ,
C    '            nh,NHE,NICONT,nj,nr,nx,
C    '            PG,XE,XG,XI,XICONT,ZE,ZG,Z,DZXI,ERROR,*9999)
C               WRITE(IO4,2004) XINEW(1),XINEW(2),Z,DZXI(1),DZXI(2)
C2004           FORMAT(' Correction: XI=',2(E12.3,4X),'  Z=',E12.3,
C    '            3X,'DZXI=',2(E12.3,4X))
C             ENDIF
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
                      CALL XIXL(IBT,IDO,INP,ISTEP,IW,NBJ,XI,XE,XL,
     '                  ERROR,*9999)
                    ENDIF
                    IF(DOP) THEN
                      CALL CFUNC(CONTYP,IBT,ICOORD,IDO,INP,NAN,NBH,NBJ,
     '                  nh,NHE,NICONT,nj,nr,nx,
     '                  PG,XE,XG,XI,XICONT,ZE,ZG,Z,DZXI,ERROR,*9999)
                      WRITE(OP_STRING,'('' Estimation: XI='','
     '                  //'2(E12.3,4X),''Z='',E12.3,3X,''DZXI='','
     '                  //'2(E12.3,4X))')XI(1),XI(2),Z,DZXI(1),DZXI(2)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
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
     '            /DENOM
                ENDIF
                XI(2)=1.d0
              ENDIF
              IF(IW.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
                XL(1,ISTEP)=XI(1)
                XL(2,ISTEP)=XI(2)
              ELSE
                CALL XIXL(IBT,IDO,INP,ISTEP,IW,NBJ,XI,XE,XL,ERROR,*9999)
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' XL(nj,'',I3,''): '',3E10.3)')
     '            nstep+1,XL(1,nstep+1),XL(2,nstep+1),XL(3,nstep+1)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                CALL CFUNC(CONTYP,IBT,ICOORD,IDO,INP,NAN,NBH,NBJ,
     '            nh,NHE,NICONT,nj,nr,nx,
     '            PG,XE,XG,XI,XICONT,ZE,ZG,Z,DZXI,ERROR,*9999)
                WRITE(OP_STRING,'('' Estimation: XI='',2(E12.3,4X),'
     '            //'''Z='',E12.3,3X,''DZXI='',2(E12.3,4X))')
     '            XI(1),XI(2),Z,DZXI(1),DZXI(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              GO TO 550
            ENDIF
            ISTEP=ISTEP+1
            IF(IW.EQ.4.AND.PROJEC(1:2).EQ.'XI') THEN
              XL(1,ISTEP)=XI(1)
              XL(2,ISTEP)=XI(2)
            ELSE
              CALL XIXL(IBT,IDO,INP,ISTEP,IW,NBJ,XI,XE,XL,ERROR,*9999)
            ENDIF
 5          CONTINUE
          ENDDO

 550      IF(IW.EQ.5) THEN
C old MPN unused?
c           CALL FBSPLC(ZVAL)
            ERROR='>>Not Implemented'
            GO TO 9999
          ENDIF
          IF(IW.EQ.13) THEN !cross-section at constant theta
            DO i=1,ISTEP
              XL(2,i)=DSQRT(XL(2,i)**2+XL(3,i)**2)
            ENDDO
            CALL POLYLINE(INDEX,IW,ISTEP,XL,ERROR,*9999)
          ELSE
            CALL POLYLINE(INDEX,IW,ISTEP,XL,ERROR,*9999)
          ENDIF

        ENDIF
 6      CONTINUE
      ENDDO
      IF(ISTEP.GT.0) THEN
        MSTEP=ISTEP/2
        DO j=1,NJT
          XLM(j)=XL(j,MSTEP)
        ENDDO
      ENDIF

      CALL EXITS('CONTOR')
      RETURN
 9999 CALL ERRORS('CONTOR',ERROR)
      CALL EXITS('CONTOR')
      RETURN 1
      END


