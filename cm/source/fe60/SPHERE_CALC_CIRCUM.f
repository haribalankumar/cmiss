      SUBROUTINE SPHERE_CALC_CIRCUM(XYZ_1,XYZ_2,XYZ_3,XYZ_CENTRE,
     '  IERROR,ERROR,*)

C#### Subroutine: SPHERE_CALC_CIRCUM
C###  Description:
C###  Calculates the circumcentre (XYZ_CENTRE) of three points
C###   (XYZ_1, XYZ_2 and XYZ_3) on the surface of a sphere.  There are
C###   two choices for the centre.  The straight line in 3-space
C###   between the two choices goes through the centre of the sphere.
C###   One is chosen based on the order of the vertices.  IERROR not
C###   zero indicates an error.

C***  Created by: David Bullivant, January 2002
C***  Last modified: 3 February 2002

      IMPLICIT NONE
!     Parameter List
      INTEGER IERROR
      REAL*8 XYZ_CENTRE(3),XYZ_1(3),XYZ_2(3),XYZ_3(3)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 XYZ(3),X12,X13,Y12,Y13,Z12,Z13

      CALL ENTERS('SPHERE_CALC_CIRCUM',*9999)

      IERROR=0
      X12=XYZ_1(1)-XYZ_2(1)
      Y12=XYZ_1(2)-XYZ_2(2)
      Z12=XYZ_1(3)-XYZ_2(3)
      X13=XYZ_1(1)-XYZ_3(1)
      Y13=XYZ_1(2)-XYZ_3(2)
      Z13=XYZ_1(3)-XYZ_3(3)
      XYZ(1)=Y12*Z13-Z12*Y13
      XYZ(2)=Z12*X13-X12*Z13
      XYZ(3)=X12*Y13-Y12*X13
      CALL SPHERE_PROJECT(XYZ,XYZ_CENTRE,IERROR,ERROR,*9999)

      CALL EXITS('SPHERE_CALC_CIRCUM')
      RETURN
 9999 CALL ERRORS('SPHERE_CALC_CIRCUM',ERROR)
      CALL EXITS('SPHERE_CALC_CIRCUM')
      RETURN 1
      END


