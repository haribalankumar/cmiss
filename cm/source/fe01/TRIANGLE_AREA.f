      REAL*8 FUNCTION TRIANGLE_AREA(A,B,C)

C#### Function: TRIANGLE_AREA
C###  Type: REAL*8
C###  Description:
C###    Returns the area of the triangle specified in 3D by points A,B,C
C**** Created by Peter Bier, April 2003

      IMPLICIT NONE
!     Parameter List
      REAL*8 A(3),B(3),C(3)
!     Local variables
      INTEGER nj
      REAL*8 AB(3),AC(3),AREA,CR(3)
!     Functions
      DO nj = 1,3
        AB(nj) = B(nj) - A(nj)
        AC(nj) = C(nj) - A(nj)
      ENDDO

      CALL CROSS(AB,AC,CR)
      AREA = DSQRT(CR(1)**2 + CR(2)**2 + CR(3)**2)
      TRIANGLE_AREA = 0.5d0 * AREA
      RETURN
      END


