      SUBROUTINE XZ(ICOORD,Y,Z)

C#### Subroutine: XZ
C###  Description:
C###    XZ performs a coordinate transformation from the coordinate
C###    system identified by ICOORD at the point with coordinates X to
C###    rectangular cartesian coordinates Z.
C**** ICOORD=1 for rectangular cartesian coordinates
C****        2 for cylindrical polar coordinates
C****        3 for spherical polar coordinates
C****        4 for prolate spheriodal coordinates
C****        5 for oblate spheroidal coordinates

      IMPLICIT NONE
      INCLUDE 'b14.cmn'
!     Parameter List
      INTEGER ICOORD
      REAL*8 Y(3),Z(3)
!     Local Variables
      REAL*8 X(3)

C     CALL ENTERS('XZ',*9999)
      X(1)=Y(1)
      X(2)=Y(2)
      X(3)=Y(3)
      GO TO (1,2,3,4,5),ICOORD
    1   Z(1)=X(1)
        Z(2)=X(2)
        Z(3)=X(3)
        GO TO 99
    2   Z(1)=X(1)*DCOS(X(2))
        Z(2)=X(1)*DSIN(X(2))
        Z(3)=X(3)
        GO TO 99
    3   Z(1)=X(1)*DCOS(X(2))*DCOS(X(3))
        Z(2)=X(1)*DSIN(X(2))*DCOS(X(3))
        Z(3)=X(1)*DSIN(X(3))
        GO TO 99
    4   Z(1)=FOCUS*DCOSH(X(1))*DCOS(X(2))
        Z(2)=FOCUS*DSINH(X(1))*DSIN(X(2))*DCOS(X(3))
        Z(3)=FOCUS*DSINH(X(1))*DSIN(X(2))*DSIN(X(3))
        GO TO 99
    5   Z(1)=FOCUS*DCOSH(X(1))*DCOS(X(2))*DCOS(X(3))
        Z(2)=FOCUS*DSINH(X(1))*DSIN(X(2))
        Z(3)=FOCUS*DCOSH(X(1))*DCOS(X(2))*DSIN(X(3))

C99   CALL EXITS('XZ')
 99   RETURN
      END


