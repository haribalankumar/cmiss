      SUBROUTINE COORD_TRANS(ANGLE_X,ANGLE_Y,ANGLE_Z,X_TRANS,
     '  TRANSLATE,X_INIT,ERROR,*)

C#### Subroutine: COORD_TRANS
C###  Description:
C###    COORD_TRANS rotates and translates a set of points.

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      REAL*8 ANGLE_X,ANGLE_Y,ANGLE_Z,X_TRANS(3),
     '  TRANSLATE(3),X_INIT(3)
      CHARACTER ERROR*(*)
!     Local variables
      REAL*8 CX,CY,CZ,SX,SY,SZ

      CALL ENTERS('COORD_TRANS',*9999)

      CX=DCOS(ANGLE_X)
      SX=DSIN(ANGLE_X)
      CY=DCOS(ANGLE_Y)
      SY=DSIN(ANGLE_Y)
      CZ=DCOS(ANGLE_Z)
      SZ=DSIN(ANGLE_Z)
      X_TRANS(1)=CZ*CY*X_INIT(1)-(CZ*SY*SX+SZ*CX)*X_INIT(2)
     '  +(CZ*SY*CX-SZ*SX)*X_INIT(3)+TRANSLATE(1)
      X_TRANS(2)=SZ*CY*X_INIT(1)-(SZ*SY*SX-CZ*CX)*X_INIT(2)
     '  +(SZ*SY*CX+CZ*SX)*X_INIT(3)+TRANSLATE(2)
      X_TRANS(3)=(-SY)*X_INIT(1)-CY*SX*X_INIT(2)+CY*CX*X_INIT(3)
     '  +TRANSLATE(3)


      CALL EXITS('COORD_TRANS')
      RETURN
 9999 CALL ERRORS('COORD_TRANS',ERROR)
      CALL EXITS('COORD_TRANS')
      RETURN 1
      END


