      SUBROUTINE CROOTS(CONTYP,IBT,IDO,INP,ITMAX,NAN,NBH,NBJ,
     '  nh,NHE,nj,ICOORD,NICONT,nr,NROOTS,nx,ROOTS,TOL,XICONT,
     '  PG,XE,XG,ZE,ZG,ZVAL,ERROR,*)

C#### Subroutine: CROOTS
C###  Description:
C###    CROOTS calculates Xi-coordinates ROOTS(1,nroot),ROOTS(2,nroot),
C###    nroot=1,..NROOTS at which the contour Z=ZVAL crosses element
C###    boundaries.  CFUNC defines the field within the element.
C###    A Newton-Raphson iteration constrained to the boundaries is used
C###    to determine the roots.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),ICOORD,IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),ITMAX,NAN(NIM,NAM,NBFM),NBH(NHM),NBJ(NJM),nh,
     '  NHE,NICONT,nj,nr,NROOTS,nx
      REAL*8 PG(NSM,NUM,NGM,NBM),ROOTS(2,12),TOL,
     '  XE(NSM,NJM),XG(NJM,NUM),XICONT,ZE(NSM,NHM),ZG(NHM,NUM),ZVAL
      CHARACTER CONTYP*(*),ERROR*(*)
!     Local Variables
      INTEGER it,n,ni,ni1,NI2,nixi,nroot,nt
      REAL*8 DXI,DXIMAX,DZXI(2),XI(2),Z,ZTOL
      CHARACTER FORMAT*500

      CALL ENTERS('CROOTS',*9999)
      NROOTS=0
      nt=4
      DXIMAX=1.d0/DBLE(nt)
      ZTOL=TOL*(1.d0+DABS(ZVAL))
      DO ni1=1,2
        NI2=MOD(ni1,2)+1
        DO nixi=0,1
          XI(ni1)=DBLE(nixi)
          DO n=1,nt
            XI(NI2)=DBLE(2*n-1)/DBLE(2*nt)
            DO it=1,ITMAX
              CALL CFUNC(CONTYP,IBT,ICOORD,IDO,INP,NAN,NBH,NBJ,
     '          nh,NHE,NICONT,nj,nr,nx,
     '          PG,XE,XG,XI,XICONT,ZE,ZG,Z,DZXI,ERROR,*9999)
              IF(DOP) THEN
                WRITE(OP_STRING,*) ' Z DZXI(1) DZXI(2)=',Z,DZXI(1),
     '            DZXI(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
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
c17_sep_88    IF((XI(NI2).LT.-1.d0).OR.(XI(NI2).GT.2.d0)) GOTO 4
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
      IF(DOP) THEN
        FORMAT='(/'' ROOTS(ni,nroot)= '',2(G12.4,4X)'
     '       //' /''                  '',2(G12.4,4X))'
        WRITE(OP_STRING,FORMAT) ((ROOTS(ni,nroot),ni=1,2),
     '    nroot=1,NROOTS)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('CROOTS')
      RETURN
 9999 CALL ERRORS('CROOTS',ERROR)
      CALL EXITS('CROOTS')
      RETURN 1
      END


