      SUBROUTINE D_ROT_COORDSYS(id,NITB,ANGLE,COORD_MATRIX,
     '  D_ANGLE,D_COORD_MATRIX,ERROR,*)

C#### Subroutine: D_ROT_COORDSYS
C###  Description:
C###    Given an input COORD_MATRIX and its NITB derivatives
C###    D_COORD_MATRIX, D_ROT_COORDSYS performs 3D orthogonal rotation on
C###    COORD_MATRIX about the id-axis and calculates the derivatives of
C###    the rotated COORD_MATRIX. NOTE this corresponds to a
C###    post-multiplication of COORD_MATRIX by the transpose of the
C###    orthogonal rotation matrix and differentiation of the product.

      IMPLICIT NONE
!     Parameter List
      INTEGER id,NITB
      REAL*8 ANGLE,COORD_MATRIX(3,3),D_ANGLE(3),D_COORD_MATRIX(3,3,3)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER id1,id2,ni,njj
      REAL*8 A(3),B(3),COS_ANGLE,SIN_ANGLE,D_A,D_B

      CALL ENTERS('D_ROT_COORDSYS',*9999)

      IF(id.EQ.1) THEN !rotate about 1st axis
        id1=2
        id2=3
      ELSE IF(id.EQ.2) THEN !rotate about 2nd axis
        id1=1
        id2=3
      ELSE IF(id.EQ.3) THEN !rotate about 3rd axis
        id1=1
        id2=2
      ENDIF

      IF(ANGLE.NE.0.0d0) THEN
        COS_ANGLE=DCOS(ANGLE)
        SIN_ANGLE=DSIN(ANGLE)
C       Rotate COORD_MATRIX
        DO njj=1,3
          A(njj)= COS_ANGLE*COORD_MATRIX(njj,id1)
     '      +SIN_ANGLE*COORD_MATRIX(njj,id2)
          B(njj)=-SIN_ANGLE*COORD_MATRIX(njj,id1)
     '      +COS_ANGLE*COORD_MATRIX(njj,id2)
          COORD_MATRIX(njj,id1)=A(njj)
          COORD_MATRIX(njj,id2)=B(njj)
        ENDDO !njj
C       Rotate D_COORD_MATRIX for each derivative
        DO ni=1,NITB
          DO njj=1,3
            D_A= COS_ANGLE*D_COORD_MATRIX(njj,id1,ni)
     '        +SIN_ANGLE*D_COORD_MATRIX(njj,id2,ni)
            D_B=-SIN_ANGLE*D_COORD_MATRIX(njj,id1,ni)
     '        +COS_ANGLE*D_COORD_MATRIX(njj,id2,ni)
            D_COORD_MATRIX(njj,id1,ni)=D_A
            D_COORD_MATRIX(njj,id2,ni)=D_B
          ENDDO !njj
        ENDDO !ni
      ELSE !zero angle
        DO njj=1,3
          A(njj)=COORD_MATRIX(njj,id1)
          B(njj)=COORD_MATRIX(njj,id2)
        ENDDO !njj
      ENDIF !angle
C     Apply the derivatives of rotation angle
      DO ni=1,NITB
        IF(D_ANGLE(ni).NE.0.0d0) THEN
          DO njj=1,3
            D_COORD_MATRIX(njj,id1,ni)=
     '        D_COORD_MATRIX(njj,id1,ni)+B(njj)*D_ANGLE(ni)
            D_COORD_MATRIX(njj,id2,ni)=
     '        D_COORD_MATRIX(njj,id2,ni)-A(njj)*D_ANGLE(ni)
          ENDDO !njj
        ENDIF !non-zero angle deriv
      ENDDO !ni

      CALL EXITS('D_ROT_COORDSYS')
      RETURN
 9999 CALL ERRORS('D_ROT_COORDSYS',ERROR)
      CALL EXITS('D_ROT_COORDSYS')
      RETURN 1
      END


