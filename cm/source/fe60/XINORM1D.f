      SUBROUTINE XINORM1D(IBT,IDO,INP,NBJ,XI2,X,XE,XI,XP_IB,ERROR,*)

C#### Subroutine: XINORM1D
C###  Description:
C###    XINORM1D calculates the XI coordinates of a point on a line.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM),XI2
      REAL*8 X(3),XE(NSM,NJM),XI(3),XP_IB(3)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nb,nj
      REAL*8 DIST,DIST_NEW,DIST1,DIST2,incr,PXI,xtemp1(3),xtemp2(3),
     '  xitemp1(3),xitemp2(3)
      LOGICAL FOUND

      CALL ENTERS('XINORM1D',*9999)

      INCR=0.01d0
      DIST=0.d0
      DIST1=0.d0
      DIST2=0.d0
      DO nj=1,3
        xitemp1(nj)=XI(nj)
        xitemp2(nj)=XI(nj)
      ENDDO !nj
      xitemp1(XI2)=xitemp1(XI2)+INCR
      xitemp2(XI2)=xitemp2(XI2)-INCR
      DO nj=1,3 !note that 3D elements assumed throughout
        nb=NBJ(nj) !basis fn
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI,XE(1,nj)) !global coordinates of surface point
        xtemp1(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    xitemp1,XE(1,nj)) !global coordinates of surface point
        xtemp2(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    xitemp2,XE(1,nj)) !global coordinates of surface point
        DIST=DIST+(X(nj)-XP_IB(nj))**2.d0
        DIST1=DIST1+(xtemp1(nj)-XP_IB(nj))**2.d0
        DIST2=DIST2+(xtemp2(nj)-XP_IB(nj))**2.d0
      ENDDO !nj
      DIST=DSQRT(DIST)
      DIST=DSQRT(DIST1)
      DIST=DSQRT(DIST2)
      FOUND=.FALSE.
      IF(DIST1.LT.DIST)THEN
        INCR=INCR
      ELSE IF(DIST2.LT.DIST)THEN
        INCR=-INCR
      ELSE
        FOUND=.TRUE.
      ENDIF
      DO WHILE(.NOT.FOUND)
        XI(XI2)=XI(XI2)+INCR
        DIST_NEW=0.d0
        DO nj=1,3 !note that 3D elements assumed throughout
          nb=NBJ(nj) !basis fn
          X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '      XI,XE(1,nj)) !global coordinates of surface point
          DIST_NEW=DIST_NEW+(X(nj)-XP_IB(nj))**2.d0
        ENDDO !nj
        DIST_NEW=DSQRT(DIST_NEW)
        IF(DIST_NEW.GT.DIST+LOOSE_TOL)THEN
          INCR=INCR*(-0.3d0)
        ENDIF
        IF(DABS(DIST-DIST_NEW).LE.LOOSE_TOL)THEN
          FOUND=.TRUE.
        ELSE
          DIST=DIST_NEW
        ENDIF
      ENDDO !WHILE

      CALL EXITS('XINORM1D')
      RETURN
 9999 CALL ERRORS('XINORM1D',ERROR)
      CALL EXITS('XINORM1D')
      RETURN 1
      END


