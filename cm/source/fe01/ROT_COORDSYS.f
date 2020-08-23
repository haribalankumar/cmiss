      SUBROUTINE ROT_COORDSYS(id,ANGLE,COORD_MATRIX,ERROR,*)

C#### Subroutine: ROT_COORDSYS
C###  Description:
C###    ROT_COORDSYS returns 3D orthogonal rotation of input
C###    COORD_MATRIX about ID-axis. NOTE this corresponds to a
C###    post-multiplication of COORD_MATRIX by the transpose of
C###    the orthogonal rotation matrix for a rotation of the coordinate
C###    system.

      IMPLICIT NONE
!     Parameter List
      INTEGER id
      REAL*8 ANGLE,COORD_MATRIX(3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER id1,id2
      REAL*8 a(3),b(3),COS_ANGLE,SIN_ANGLE

      CALL ENTERS('ROT_COORDSYS',*9999)
      IF(id.EQ.1) THEN !rotate about 1st axis
        id1=2
        id2=3
      ELSE IF(id.EQ.2) THEN !rotate about 2nd axis
C KAT 30Nov00: Changed to rotate about 2nd axis instead of -2nd axis.
        id1=3
        id2=1
C        id1=1
C        id2=3
      ELSE IF(id.EQ.3) THEN !rotate about 3rd axis
        id1=1
        id2=2
      ENDIF
      COS_ANGLE=DCOS(ANGLE)
      SIN_ANGLE=DSIN(ANGLE)
      a(1)= COS_ANGLE*COORD_MATRIX(1,id1)+SIN_ANGLE*COORD_MATRIX(1,id2)
      a(2)= COS_ANGLE*COORD_MATRIX(2,id1)+SIN_ANGLE*COORD_MATRIX(2,id2)
      a(3)= COS_ANGLE*COORD_MATRIX(3,id1)+SIN_ANGLE*COORD_MATRIX(3,id2)
      b(1)=-SIN_ANGLE*COORD_MATRIX(1,id1)+COS_ANGLE*COORD_MATRIX(1,id2)
      b(2)=-SIN_ANGLE*COORD_MATRIX(2,id1)+COS_ANGLE*COORD_MATRIX(2,id2)
      b(3)=-SIN_ANGLE*COORD_MATRIX(3,id1)+COS_ANGLE*COORD_MATRIX(3,id2)
      COORD_MATRIX(1,id1)=a(1)
      COORD_MATRIX(2,id1)=a(2)
      COORD_MATRIX(3,id1)=a(3)
      COORD_MATRIX(1,id2)=b(1)
      COORD_MATRIX(2,id2)=b(2)
      COORD_MATRIX(3,id2)=b(3)

      CALL EXITS('ROT_COORDSYS')
      RETURN
 9999 CALL ERRORS('ROT_COORDSYS',ERROR)
      CALL EXITS('ROT_COORDSYS')
      RETURN 1
      END


