      SUBROUTINE FIND_TETRAHEDRAL_TAIL(TETRA,ERROR,*)

C#### Subroutine: FIND_TETRAHEDRAL_TAIL
C###  Description:
C###    Determines the tail of the tetrahedrals. Genmesh routine.
CC JMB 19-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER TETRA(0:LDTETRA,4)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER REF

      CALL ENTERS('FIND_TETRAHEDRAL_TAIL',*9999)

      IF(TETRA(0,HEAD).NE.0)THEN
        REF=TETRA(0,HEAD)
        DO WHILE(TETRA(REF,NEXT).NE.0)
          REF=TETRA(REF,NEXT)
        ENDDO !while(tetra(ref,next))
        TETRA(0,TAIL)=REF
      ELSE
        TETRA(0,TAIL)=0
      ENDIF !tetra(0,head)

      CALL EXITS('FIND_TETRAHEDRAL_TAIL')
      RETURN
 9999 CALL ERRORS('FIND_TETRAHEDRAL_TAIL',ERROR)
      CALL EXITS('FIND_TETRAHEDRAL_TAIL')
      RETURN 1
      END


