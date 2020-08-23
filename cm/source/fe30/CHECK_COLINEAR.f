      SUBROUTINE CHECK_COLINEAR(POINT1,POINT2,POINT3,COLINEAR,ERROR,*)

C#### Subroutine: CHECK_COLINEAR
C###  Description:
C###    CHECK_COLINEAR checks whether two vectors are colinear.
C***  Adapted from PLANE_FROM_3_PTS, Mar 2001

      IMPLICIT NONE

      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter list
      REAL*8 POINT1(3),POINT2(3),POINT3(3)
      LOGICAL COLINEAR
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nj
      REAL*8 ERR1(3),ERR2(3),LU,LV,U(3),V(3)

      CALL ENTERS('CHECK_COLINEAR',*9999)

      COLINEAR=.FALSE.
      LU=0.d0
      LV=0.d0
      DO nj=1,NJT
        U(nj)=POINT2(nj)-POINT1(nj)
        V(nj)=POINT3(nj)-POINT1(nj)
        LU=LU+U(nj)**2.d0
        LV=LV+V(nj)**2.d0
      ENDDO !nj
      LU=DSQRT(LU)
      LV=DSQRT(LV)
      ! GDR 11/2/05 if 2 of the points are the same then LU and LV
      ! can be zero causing div by zero below and resulting in
      ! the wrong answer (on Linux) 
      IF((LU.NE.0).AND.(LV.NE.0)) THEN
        DO nj=1,NJT
          ERR1(nj)=DABS(U(nj)/LU-V(nj)/LV)
          ERR2(nj)=DABS(U(nj)/LU+V(nj)/LV)
        ENDDO !nj
        IF((ERR1(1).LE.ZERO_TOL.AND.ERR1(2).LE.ZERO_TOL.AND.ERR1(3).LE.
     '   ZERO_TOL).OR.(ERR2(1).LE.ZERO_TOL.AND.ERR2(2).LE.ZERO_TOL.AND.
     '   ERR2(3).LE.ZERO_TOL)) COLINEAR=.TRUE.
      ELSE
        COLINEAR=.TRUE.
      ENDIF 
      
      CALL EXITS('CHECK_COLINEAR')
      RETURN
 9999 CALL ERRORS('CHECK_COLINEAR',ERROR)
      CALL EXITS('CHECK_COLINEAR')
      RETURN 1
      END

