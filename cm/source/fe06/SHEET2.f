      SUBROUTINE SHEET2(INDEX,IW,ZD,ERROR,*)

C#### Subroutine: SHEET2
C###  Description:
C###    SHEET2 draws measured sheet angles on constant theta plane
C###    (CALC_SHEET=.FALSE.).

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER INDEX,IW
      REAL*8 ZD(NJM)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 ANGLE,DX,POINTS(3,2),RAD,X

      CALL ENTERS('SHEET2',*9999)
      CALL ASSERT(IW.EQ.13,
     '  ' Incorrect workstation ID: sb 13 for sheets',ERROR,*9999)
      RAD=DSQRT(ZD(2)**2+ZD(3)**2)
      X=ZD(1)
      ANGLE=ZD(6)
      WRITE(OP_STRING,'('' rad='',E12.3,'' x='',E12.3,'' angle='','
     '  //'E12.3)') RAD,X,ANGLE
      CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      DX=0.03D0*FOCUS
      POINTS(1,1)=X  -DX*DSIN(ANGLE)
      POINTS(2,1)=RAD+DX*DCOS(ANGLE)
      POINTS(1,2)=X  +DX*DSIN(ANGLE)
      POINTS(2,2)=RAD-DX*DCOS(ANGLE)
      CALL POLYLINE(INDEX,IW,2,POINTS,ERROR,*9999)

      CALL EXITS('SHEET2')
      RETURN
 9999 CALL ERRORS('SHEET2',ERROR)
      CALL EXITS('SHEET2')
      RETURN 1
      END


