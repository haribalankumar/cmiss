      SUBROUTINE MAKE_TETRAHEDRAL(ELEM,FACE,IREGION,TETRA,CIRC,
     '   NODE,RSQ,PLANE,ERROR,*)

C#### Subroutine: MAKE_TETRAHEDRAL
C###  Description:
C###    Creates an initial tetrahedral enclosing the boundary object.
C###    Genemesh routine.
CC JMB 15-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4), FACE(LDFACE,NFACES), IREGION(0:N_GM),
     '  TETRA(0:LDTETRA,4)
      REAL*8 CIRC(LDCIRC,3), NODE(LDNODE,3), RSQ(N_GM)
      LOGICAL PLANE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J
      REAL*8 MAXVAL(3),MINVAL(3)

      CALL ENTERS('MAKE_TETRAHEDRAL',*9999)

      ! Determine boundary limits
      DO J=1,NJT
        MAXVAL(J)=NODE(1,J)
        MINVAL(J)=MAXVAL(J)

        DO I=2,IREGION(0) ! npts
          MAXVAL(J)=MAX(MAXVAL(J),NODE(I,J))
          MINVAL(J)=MIN(MINVAL(J),NODE(I,J))
        ENDDO !i
      ENDDO !j
      ! Construct the initial tetrahedral
      I=IREGION(0)+1
      DO J=1,NJT
        NODE(I,J)=MINVAL(J)-3.d0*(MAXVAL(J)-MINVAL(J))
      ENDDO !j
      I=I+1
      NODE(I,1)=MAXVAL(1)+3.0d0*(MAXVAL(1)-MINVAL(1))
      NODE(I,2)=MINVAL(2)-3.0d0*(MAXVAL(2)-MINVAL(2))
      NODE(I,3)=MINVAL(3)-3.0d0*(MAXVAL(3)-MINVAL(3))
      I=I+1
      NODE(I,1)=0.5d0*(MAXVAL(1)+MINVAL(1))
      NODE(I,2)=MAXVAL(2)+3.0d0*(MAXVAL(2)-MINVAL(2))
      NODE(I,3)=MINVAL(3)-3.0d0*(MAXVAL(3)-MINVAL(3))
      I=I+1
      NODE(I,1)=0.5d0*(MAXVAL(1)+MINVAL(1))
      NODE(I,2)=0.5d0*(MAXVAL(2)+MINVAL(2))
      NODE(I,3)=MAXVAL(3)+4.0d0*(MAXVAL(3)-MINVAL(3))
      ! Initialize the tetrahedral
      IELEM=0
      CALL ADD_TETRAHEDRAL(TETRA,ERROR,*9999)
      ! Assign the global nodes to the created tetrahedral
      DO J=1,NPTS
        ELEM(IELEM,J)=(I-NPTS)+J
        FACE(IELEM,J)=0 ! Initialize face
      ENDDO !j
      CALL GET_CIRCUMDATA(ELEM(IELEM,1),CIRC(IELEM,1),NODE,
     '  RSQ(IELEM),PLANE,ERROR,*9999)! Determine the circumdata
      DO J=IREGION(0)+1,I! Update
        IREGION(J)=0 ! Temporary
      ENDDO !j
      IREGION(0)=I

      CALL EXITS('MAKE_TETRAHEDRAL')
      RETURN
 9999 CALL ERRORS('MAKE_TETRAHEDRAL',ERROR)
      CALL EXITS('MAKE_TETRAHEDRAL')
      RETURN 1
      END


