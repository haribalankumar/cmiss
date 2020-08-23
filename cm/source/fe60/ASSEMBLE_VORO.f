      SUBROUTINE ASSEMBLE_VORO(ELEM,FACE,FOUND,INODE,TETRA,CIRC,NODE,
     '  RSQ,ERROR,*)

C#### Subroutine: ASSEMBLE_VORO
C###  Description:
C###    Genmesh routine.
CC JMB 15-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),FACE(LDFACE,NFACES),FOUND(0:N_GM),
     '  INODE,TETRA(0:LDTETRA,4)
      REAL*8 CIRC(LDCIRC,3),NODE(LDNODE,3),RSQ(N_GM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER INEXT,REF
!     Local Functions
      INTEGER SEARCH_NEXT

      CALL ENTERS('ASSEMBLE_VORO',*9999)

      ! Determine the tetrahedral containing the INODE node
      INEXT=1
      REF=TETRA(0,HEAD)
      DO WHILE(INEXT.NE.0)
        INEXT=SEARCH_NEXT(ELEM(REF,1),INODE,NODE)
        IF(INEXT.NE.0)THEN
          REF=FACE(REF,INEXT)
        ENDIF !inext
      ENDDO !while inext
      FOUND(0)=0 ! Initialize search
      ! Recursive tree search phase
      CALL CHECK_TETRAHEDRAL(FACE,FOUND,REF,TETRA,CIRC,
     '  NODE(INODE,1),RSQ,ERROR,*9999)

      CALL EXITS('ASSEMBLE_VORO')
      RETURN
 9999 CALL ERRORS('ASSEMBLE_VORO',ERROR)
      CALL EXITS('ASSEMBLE_VORO')
      RETURN 1
      END


