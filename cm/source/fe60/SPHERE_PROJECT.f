      SUBROUTINE SPHERE_PROJECT(RAW_XYZ,XYZ,IERROR,ERROR,*)

C#### Subroutine: SPHERE_PROJECT
C###  Description:
C###    Projects a point (RAW_XYZ) onto the surface of a unit sphere
C###    (XYZ).  IERROR not zero indicates an error.

C***  Created by: David Bullivant, January 2002
C***  Last modified: 14 January 2002

      IMPLICIT NONE
!     Parameter List
      INTEGER IERROR
      REAL*8 RAW_XYZ(3),XYZ(3)
      REAL*8 RADIUS
      CHARACTER ERROR*(*)

      CALL ENTERS('SPHERE_PROJECT',*9999)

      IERROR=0
      RADIUS=DSQRT(RAW_XYZ(1)*RAW_XYZ(1)+RAW_XYZ(2)*RAW_XYZ(2)+
     '  RAW_XYZ(3)*RAW_XYZ(3))
      IF(RADIUS.GT.0.d0) THEN
        XYZ(1)=RAW_XYZ(1)/RADIUS
        XYZ(2)=RAW_XYZ(2)/RADIUS
        XYZ(3)=RAW_XYZ(3)/RADIUS
      ELSE
        XYZ(1)=1.d0
        XYZ(2)=0.d0
        XYZ(3)=0.d0
        IERROR=1
      ENDIF

      CALL EXITS('SPHERE_PROJECT')
      RETURN
 9999 CALL ERRORS('SPHERE_PROJECT',ERROR)
      CALL EXITS('SPHERE_PROJECT')
      RETURN 1
      END


