      SUBROUTINE FLIP_EDGE_TO_FACE(AXIS,ELEM,FACE,ITETRA,TETRA,
     '  CIRC,NODE,RSQ,PLANE,ERROR,*)

C#### Subroutine: FLIP_EDGE_TO_FACE
C###  Description:
C###    Determines the 3 periphery nodes, stores the face
C###    neighbours data structure, removes the three tetrahedrals from
C###    the triangulation, adds two new tetrahedral amd uses the face
C###    neighbours to restore neighbour integrity.
CC JMB 18-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER AXIS(2),ELEM(LDELEM,4),FACE(LDFACE,NFACES),
     '  ITETRA(3),TETRA(0:LDTETRA,4)
      REAL*8 CIRC(LDCIRC,3),NODE(LDNODE,3),RSQ(N_GM)
      LOGICAL PLANE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,K,NEIB,NEIGH(2,3),NPER,OPP(4),PER(4),
     '  PTR(2,3),IPTR(2,3),REF,WORK(2)
      LOGICAL BOOL

      CALL ENTERS('FLIP_EDGE_TO_FACE',*9999)

      NPER=1
      DO I=1,2
        DO J=1,NPTS
          REF=ELEM(ITETRA(I),J)
          IF((REF.NE.AXIS(1)).AND.(REF.NE.AXIS(2)))THEN
            PER(NPER)=REF
            NPER=NPER+1 ! Increment counter
          ENDIF !ref
        ENDDO !j
      ENDDO !i
      ! If the second node of the tetrahedral B was in tetrahedral A,
      !  use the third node instead
      IF((PER(1).EQ.PER(3)).OR.(PER(2).EQ.PER(3))) PER(3)=PER(4)
      ! We know that PER(3) is opposite A, therefore we need to
      !  establish which of PER(1) and PER(2) are opposite B and C
      OPP(1)=PER(3)
      REF=PER(1)
      BOOL=.TRUE.
      DO I=1,NPTS
        BOOL=BOOL.AND.(ELEM(ITETRA(2),I).NE.REF)
      ENDDO !i
      IF(BOOL)THEN
        OPP(2)=REF
        OPP(3)=PER(2)
      ELSE
        OPP(2)=PER(2)
        OPP(3)=REF
      ENDIF !bool
      ! Store the face neighbour data structure
      DO I=1,NJT
        DO J=1,2
          DO K=1,NPTS
            IF(ELEM(ITETRA(I),K).EQ.AXIS(J))THEN
              NEIB=FACE(ITETRA(I),K)
            ENDIF !elem
          ENDDO !k
          NEIGH(NJT-J,I)=NEIB
          IF(NEIB.NE.0)THEN
            DO K=1,NFACES
              IF(FACE(NEIB,K).EQ.ITETRA(I))THEN
                PTR(J,I)=NEIB
                IPTR(I,J)=K
              ENDIF !face
            ENDDO !k
          ENDIF !neib
        ENDDO !j
      ENDDO !i
      ! Delete the tetrahedral
      DO I=1,NJT
        CALL DELETE_TETRAHEDRAL(ITETRA(I),TETRA,ERROR,*9999)
      ENDDO !i
      ! Create the two new tetrahedra
      DO I=1,2
       CALL ADD_TETRAHEDRAL(TETRA,ERROR,*9999)
       ELEM(IELEM,1)=AXIS(I)! Element definition
       DO J=2,NPTS
         ELEM(IELEM,J)=OPP(J-1)
       ENDDO !j
       CALL GET_CIRCUMDATA(ELEM(IELEM,1),CIRC(IELEM,1),NODE,
     '   RSQ(IELEM),PLANE,ERROR,*9999) ! Calculate circumdata
       WORK(I)=IELEM
      ENDDO !i
      ! Update the other neighbours
      DO I=1,2
        ! Face definition
        FACE(WORK(I),1)=WORK(NJT-I)
        FACE(WORK(I),2)=NEIGH(I,1)
        DO J=2,NFACES
          FACE(WORK(I),J)=NEIGH(I,J-1)
          IF(PTR(I,J).NE.0) FACE(PTR(I,J),IPTR(I,J))=WORK(I)
        ENDDO !j
      ENDDO !i

      CALL EXITS('FLIP_EDGE_TO_FACE')
      RETURN
 9999 CALL ERRORS('FLIP_EDGE_TO_FACE',ERROR)
      CALL EXITS('FLIP_EDGE_TO_FACE')
      RETURN 1
      END


