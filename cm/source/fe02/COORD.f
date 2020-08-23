      SUBROUTINE COORD(ICOOR0,ICOOR1,X0,X1,ERROR,*)

C#### Subroutine: COORD
C###  Description:
C###    <HTML>
C###    COORD transforms the point coordinates X0 in the coordinate
C###    system ICOOR0 to X1 in the coordinate system ICOOR1, where
C###    ICOOR0 and ICOOR1 are:
C###    <PRE>
C###      1 rectangular cartesian coords (x,y,z)
C###      2 cylindrical polar coords  (r,theta,z)
C###      3 spherical polar coords    (r,theta,phi)
C###      4 prolate-spheroidal coords (lamda,mu,theta)
C###      5 oblate -spheroidal coords (lamda,mu,theta)
C###    </PRE> </HTML>

      IMPLICIT NONE
!     Parameter List
      INTEGER ICOOR0,ICOOR1
      REAL*8 X0(3),X1(3)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 Z(3)

      CALL ENTERS('COORD',*9999)
      IF(ICOOR0.EQ.ICOOR1) THEN
        X1(1)=X0(1)
        X1(2)=X0(2)
        X1(3)=X0(3)
      ELSE IF(ICOOR0.EQ.1) THEN
        CALL ZX(ICOOR1,X0,X1)
      ELSE IF(ICOOR1.EQ.1) THEN
        CALL XZ(ICOOR0,X0,X1)
      ELSE
        CALL XZ(ICOOR0,X0,Z)
        CALL ZX(ICOOR1,Z,X1)
      ENDIF

      CALL EXITS('COORD')
      RETURN
 9999 CALL ERRORS('COORD',ERROR)
      CALL EXITS('COORD')
      RETURN 1
      END


