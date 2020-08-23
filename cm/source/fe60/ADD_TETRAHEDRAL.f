      SUBROUTINE ADD_TETRAHEDRAL(TETRA,ERROR,*)

C#### Subroutine: ADD_TETRAHEDRAL
C###  Description:
C###    Adds a tethrahedral (element) to the data set. Genmesh routine.
CC JMB 16-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER TETRA(0:LDTETRA,4)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER REF

      CALL ENTERS('ADD_TETRAHEDRAL',*9999)

      ! IELEM is the index to a newly defined tetrahedral
      IELEM=IELEM+1
      CALL ASSERT_VORO(IELEM.LE.LDTETRA,
     &  'Indexed past array bounds - increase LDTETRA in genmesh.cmn',
     &  ERROR,*9999)
      TETRA(IELEM,NEXT)=0
      IF(TETRA(0,HEAD).EQ.0)THEN
        ! No existing defined tetrahedral
        ! Initialization
        TETRA(0,HEAD)=IELEM
        TETRA(0,TAIL)=IELEM
        TETRA(IELEM,PREV)=0
        TETRA(IELEM,NEXT)=0
      ELSE
        REF=TETRA(0,TAIL)
        TETRA(REF,NEXT)=IELEM
        TETRA(0,TAIL)=IELEM
        ! Append the ordering to the stack
        TETRA(IELEM,PREV)=REF
      ENDIF !tetra

      CALL EXITS('ADD_TETRAHEDRAL')
      RETURN
 9999 CALL ERRORS('ADD_TETRAHEDRAL',ERROR)
      CALL EXITS('ADD_TETRAHEDRAL')
      RETURN 1
      END


