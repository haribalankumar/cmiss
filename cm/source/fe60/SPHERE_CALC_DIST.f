      SUBROUTINE SPHERE_CALC_DIST(XYZ_1,XYZ_2,DISTANCE,IERROR,ERROR,*)

C#### Subroutine: SPHERE_CALC_DIST
C###  Description:
C###    Calculates the DISTANCE between two points, XYZ_1 and XYZ_2,
C###    on the surface of a unit sphere.  The DISTANCE between them is
C###    the shorter distance along the great circle through them.
C###    IERROR not zero indicates an error.

C***  Created by: David Bullivant, January 2002
C***  Last modified: 14 February 2002

      IMPLICIT NONE
!     Parameter List
      INTEGER IERROR
      REAL*8 DISTANCE,XYZ_1(3),XYZ_2(3)
      CHARACTER ERROR*(*)

      CALL ENTERS('SPHERE_CALC_DIST',*9999)

      IERROR=0
C     dot product
      DISTANCE=XYZ_1(1)*XYZ_2(1)+XYZ_1(2)*XYZ_2(2)+XYZ_1(3)*XYZ_2(3)
      DISTANCE=ACOS(DISTANCE)
      IF (DISTANCE.LT.0.d0) THEN
        DISTANCE= -DISTANCE
      ENDIF

      CALL EXITS('SPHERE_CALC_DIST')
      RETURN
 9999 CALL ERRORS('SPHERE_CALC_DIST',ERROR)
      CALL EXITS('SPHERE_CALC_DIST')
      RETURN 1
      END


