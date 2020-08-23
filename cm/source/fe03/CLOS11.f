      SUBROUTINE CLOS11(IBT,IDO,INP,IT,ITMAX,NBJ,SQ,XE,XI,ZD,
     '  INELEM,ERROR,*)

C#### Subroutine: CLOS11
C###  Description:
C###    CLOS11 finds the Xi-coords at the closest approach of a 1D
C###    element to a data point with coordinates XD using a modified
C###    Newton algorithm.
C###    If INELEM is true then the closest point within the element is
C###    returned; if false then the projection from the element to the
C###    coordinate must be orthogonal.  If such a projection can't be
C###    obtained the xi position returned is outside the element.
C**** ITMAX is the maximum number of iterations.
C**** TOL is the required tolerance of the solution.
C**** VMAX denotes the maximum step length per iteration.
C**** NOTE: Works only for 1-D elements in rectangular cartesian coords.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IT,ITMAX,NBJ(NJM)
      REAL*8 SQ,XE(NSM,NJM),XI(1),ZD(NJM)
      LOGICAL INELEM
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER it1,it2,nb,nj
      REAL*8 D2SQXI,D2ZXI(3),DSQXI,DSQXIV,DZ(3),DZXI(3),PXI,SQLIN,
     '  SQOLD,TOL,V,VMAX,XILIN(1),Z(3)
      DATA VMAX /0.25D0/

      CALL ENTERS('CLOS11',*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(A)') ' >CLOS11 1-D rectangular cartesian'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      TOL=LOOSE_TOL
C      TOL2=TOL**2
      SQ=0.d0
      DO nj=1,NJT
        nb=NBJ(nj)
        Z(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '    XE(1,nj))
        DZ(nj)=Z(nj)-ZD(nj)
        SQ=SQ+DZ(nj)**2
      ENDDO
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(/'' ZD(nj)='',2(E12.6,4X),/'' XI='','
     '    //'(E12.6,4X),/''  Z(nj)='',2(E12.6,4X),/''  SQ='',E12.6)')
     '    (ZD(nj),nj=1,2),XI(1),(Z(nj),nj=1,2),SQ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      DO it1=1,ITMAX
        DSQXI =0.d0
        D2SQXI=0.d0
        DO nj=1,NJT
          nb=NBJ(nj)
          DZXI(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      nb,2,XI,XE(1,nj))
          DSQXI=DSQXI+DZXI(nj)*DZ(nj)
        ENDDO
        IF(INELEM) THEN
          IF(XI(1).EQ.0.0d0) THEN
            IF(DSQXI.GT.0.0d0) THEN
              DSQXI=0.0d0
              D2SQXI=1.0d0 !any positive number will do
            ENDIF
          ELSE IF(XI(1).EQ.1.0d0) THEN
            IF(DSQXI.LT.0.0d0) THEN
              DSQXI=0.0d0
              D2SQXI=1.0d0 !any positive number will do
            ENDIF
          ENDIF
        ENDIF
        IF(D2SQXI.EQ.0.0d0) THEN
          DO nj=1,NJT
            nb=NBJ(nj)
            D2ZXI(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,3,XI,XE(1,nj))
            D2SQXI=D2SQXI+DZXI(nj)*DZXI(nj)+D2ZXI(nj)*DZ(nj)
          ENDDO
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' DZXI(nj): '',2E11.4,'' D2ZXI(nj): '','
     '      //'2E11.4)') (DZXI(nj),nj=1,2),(D2ZXI(nj),nj=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(DABS(DSQXI).GT.D2SQXI*VMAX) THEN
          V=-DSIGN(VMAX,DSQXI)
        ELSE
          V=-DSQXI/D2SQXI
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(''    IT1='',I4,10X,''DSQXI='',E12.6,'
     '      //'/21X,''D2SQXI='',E12.6,/26X,''V='',E12.6)')
     '      it1,DSQXI,D2SQXI,V
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
C ***   Performs line search: fit a quadratic with the value & deriv at
C ***   current point IT1 & value at point IT2, then choose minimum point
        XILIN(1)=XI(1)+V
        IF(INELEM) THEN
          IF(XILIN(1).LT.0.0d0) THEN
            XILIN(1)=0.0d0
            V=-XI(1)
          ELSE IF(XILIN(1).GT.1.0d0) THEN
            XILIN(1)=1.0d0
            V=1.0d0-XI(1)
          ENDIF
        ENDIF
        SQLIN=0.d0
        DO nj=1,NJT
          nb=NBJ(nj)
          Z(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XILIN,
     '      XE(1,nj))
          DZ(nj)=Z(nj)-ZD(nj)
          SQLIN=SQLIN+DZ(nj)**2
        ENDDO
        DO it2=1,ITMAX
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(4X,''IT2='',I4,14X,''V='',E12.6,/21X,'
     '        //'''XILIN='',E12.6,/22X,''Z(nj)='',2(E12.6,4X)/22X,'
     '        //'''SQLIN='',E12.6)') it2,V,XILIN(1),(Z(nj),nj=1,2),SQLIN
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(SQLIN.LE.SQ) GOTO 5
          DSQXIV=DSQXI*V
C KAT 17Nov98: I think this is better as SUM is always -ve here
C         Added factor of half as this minimizes quadratic
          V=0.5d0*V*DSQXIV/(DSQXIV+SQ-SQLIN)
C          SUM=DSQXIV+SQ-SQLIN
C          IF(DABS(SUM).GT.1.d-5) V=V*DSQXIV/SUM
          XILIN(1)=XI(1)+V
          SQLIN=0.d0
          DO nj=1,NJT
            nb=NBJ(nj)
            Z(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XILIN,
     '        XE(1,nj))
            DZ(nj)=Z(nj)-ZD(nj)
            SQLIN=SQLIN+DZ(nj)**2
          ENDDO
        ENDDO

 5      SQOLD=SQ
        SQ=SQLIN
        XI(1)=XILIN(1)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(24X,''XI='',E12.6)') XI(1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IT=it1
        IF(.NOT.INELEM.AND.(XI(1).LT.-0.5D0.OR.XI(1).GT.1.5D0))
     '    GO TO 9998
C KAT 18Nov98: xi norm is already in unit sized dimensions
C        DXI=V
C        IF(((DXI**2)/(1.d0+XI(1)**2).LE.TOL2)
C     '      .AND.((SQOLD-SQ)/(1.d0+SQ).LE.TOL)) GO TO 9998
        IF(V.LE.TOL.AND.(SQOLD-SQ)/(1.d0+SQ).LE.TOL) GO TO 9998
      ENDDO

 9998 CALL EXITS('CLOS11')
      RETURN
 9999 CALL ERRORS('CLOS11',ERROR)
      CALL EXITS('CLOS11')
      RETURN 1
      END


