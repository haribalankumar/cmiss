      SUBROUTINE INTERNAL_NODES(ELEM,FACE,FOUND,IREGION,TETRA,
     '  CIRC,GAPSQ,NODE,RSQ,ERROR,*)

C#### Subroutine: INTERNAL_NODES
C###  Description:
C###    Automatically generates internal nodes. Genmesh.
CC JMB 20-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),FACE(LDFACE,NFACES),FOUND(0:N_GM),
     '  IREGION(0:N_GM),TETRA(0:LDTETRA,4)
      REAL*8 CIRC(LDCIRC,3),GAPSQ,NODE(LDNODE,3),RSQ(N_GM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER COUNT,I,J,INODE,REF
      REAL*8 FSPACE(0:N_GM)
      LOGICAL PLANE,REFINE
!     Local Functions
      LOGICAL DEGREE,ENCASPULATED,INTERNAL,SPACE

      ! Initialization
      INODE=IREGION(0)
      REFINE=.TRUE.
      COUNT=0
      CALL ENTERS('INTERNAL_NODES',*9999)
      CALL CALCULATE_FSPACE(ELEM,IREGION,TETRA,FSPACE,NODE,ERROR,*9999)
      DO WHILE(REFINE)
        REFINE=.FALSE.
        COUNT=COUNT+1
        TETRA(0,NEXT)=TETRA(0,HEAD)
        DO WHILE(TETRA(0,NEXT).NE.0)! Increment
          REF=TETRA(0,NEXT)
          TETRA(0,NEXT)=TETRA(REF,NEXT)
          IF(INTERNAL(ELEM(REF,1),IREGION).AND.
     '      ENCASPULATED(ELEM(REF,1),CIRC(REF,1),NODE).AND.
     '      SPACE(ELEM(REF,1),CIRC(REF,1),FSPACE,GAPSQ,NODE,
     '      RSQ(REF)).AND.DEGREE(FACE,REF,CIRC,RSQ))THEN
            ! Add the internal node - centred at the circumcentre
            INODE=INODE+1
            CALL DCOPY(NJT,CIRC(REF,1),LDCIRC,NODE(INODE,1),LDNODE)
            FSPACE(INODE)=FSPACE(0)
            ! Generate insertion polygon
            FOUND(0)=0
            CALL CHECK_TETRAHEDRAL(FACE,FOUND,REF,TETRA,CIRC,
     '        NODE(INODE,1),RSQ,ERROR,*9999)
            ! Check and discard insertion polygon for boundary violation
            DO I=1,FOUND(0)
              REF=FOUND(I)
              DO J=1,NPTS
                IF(IREGION(ELEM(REF,J)).EQ.1)THEN
                  GOTO 10
                ENDIF !iregion
              ENDDO !j
            ENDDO !i
            ! All conditions for refinement have been meet
            REFINE=.TRUE.
            CALL ADD_NODE(ELEM,FACE,FOUND,INODE,TETRA,CIRC,NODE,
     '         RSQ,PLANE,ERROR,*9999)
            IF(PLANE)THEN
              CALL FIX_TETRAHEDRAL(ELEM,FACE,TETRA,CIRC,NODE,RSQ,
     '          PLANE,ERROR,*9999)
              CALL CHECK_INTEGRITY(ELEM,FACE,TETRA,CIRC,NODE,RSQ,
     '          ERROR,*9999)
            ENDIF !plane
          ENDIF !internal
 10     ENDDO !while tetra
      ENDDO !while refine
      ! Update
      IREGION(0)=INODE

      CALL EXITS('INTERNAL_NODES')
      RETURN
 9999 CALL ERRORS('INTERNAL_NODES',ERROR)
      CALL EXITS('INTERNAL_NODES')
      RETURN 1
      END


