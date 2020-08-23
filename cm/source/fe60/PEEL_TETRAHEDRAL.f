      SUBROUTINE PEEL_TETRAHEDRAL(FACE,INDEX,TETRA,ERROR,*)

C#### Subroutine: PEEL_TETRAHEDRAL
C###  Description:
C###    Removes the peeled tetrahedral from the list. Genmesh routine.
CC JMB 19-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER FACE(LDFACE,NFACES),INDEX,TETRA(0:LDTETRA,4)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,J,NEIB,REF

      CALL ENTERS('PEEL_TETRAHEDRAL',*9999)

      IF(INDEX.EQ.TETRA(0,HEAD))THEN
        TETRA(0,HEAD)=TETRA(INDEX,NEXT)
      ELSE
        REF=TETRA(0,HEAD)
        DO WHILE((REF.NE.0).AND.(TETRA(REF,NEXT).NE.INDEX))
          TETRA(REF,NEXT)=TETRA(TETRA(REF,NEXT),NEXT)
          REF=TETRA(REF,NEXT)
        ENDDO !while
      ENDIF !index
      DO I=1,NFACES
        NEIB=FACE(INDEX,I)
        IF(NEIB.NE.0)THEN
          DO J=1,NFACES
            IF(FACE(NEIB,J).EQ.INDEX)THEN
              FACE(NEIB,J)=0
            ENDIF !face
          ENDDO !j
        ENDIF !neib
      ENDDO !i

      CALL EXITS('PEEL_TETRAHEDRAL')
      RETURN
 9999 CALL ERRORS('PEEL_TETRAHEDRAL',ERROR)
      CALL EXITS('PEEL_TETRAHEDRAL')
      RETURN 1
      END


