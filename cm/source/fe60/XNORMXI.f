      SUBROUTINE XNORMXI(IBT,IDO,INP,in1,in2,NBJ,X,XE,XI,XNORM,
     '  ERROR,*)

C#### Subroutine: XNORMXI
C###  Description:
C###    XNORMXI calculates global coordinates (X) and normal (XNORM)
C###    for a point (XI) on a face (in1-in2 face).

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBJ(NJM),in1,in2
      REAL*8 X(3),XE(NSM,NJM),XI(3),XNORM(3)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nb,nj
      REAL*8 DIST1,DIST2,DX(3,3),PXI,X05(3),XI05(3)

      CALL ENTERS('XNORMXI',*9999)

      DO nj=1,3
        XI05(nj)=0.5d0
      ENDDO !nj
      DO nj=1,3 !note that 3D elements assumed throughout
        nb=NBJ(nj) !basis fn
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI,XE(1,nj)) !global coordinates of surface point
        X05(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '    XI05,XE(1,nj)) !global coordinates of mid-point
        DX(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '    nb,2,XI,XE(1,nj)) !global Xi1 derivatives
        DX(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '    nb,4,XI,XE(1,nj)) !global Xi2 derivatives
        DX(nj,3)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '    nb,7,XI,XE(1,nj)) !global Xi3 derivatives
      ENDDO !nj
C***  Calculate unit outward normal to boundary surface
      XNORM(1)=(DX(2,in1)*DX(3,in2)-DX(3,in1)*DX(2,in2))
      XNORM(2)=(DX(3,in1)*DX(1,in2)-DX(1,in1)*DX(3,in2))
      XNORM(3)=(DX(1,in1)*DX(2,in2)-DX(2,in1)*DX(1,in2))
      CALL NORMALISE(3,XNORM,ERROR,*9999) !normalise XNORM
      DIST1=0.d0
      DIST2=0.d0
      DO nj=1,3
        DIST1=DIST1+(X(nj)+XNORM(nj)-X05(nj))**2.d0
        DIST2=DIST2+(X(nj)-XNORM(nj)-X05(nj))**2.d0
      ENDDO !nj
      DIST1=DSQRT(DIST1)
      DIST2=DSQRT(DIST2)
      IF(DIST1.LT.DIST2)THEN !determine outwards direction
        DO nj=1,3
          XNORM(nj)=-XNORM(nj)
        ENDDO !nj
      ENDIF

      CALL EXITS('XNORMXI')
      RETURN
 9999 CALL ERRORS('XNORMXI',ERROR)
      CALL EXITS('XNORMXI')
      RETURN 1
      END









