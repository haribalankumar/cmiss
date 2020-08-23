      SUBROUTINE PEEL_AWAY(ELEM, FACE, IREGION, TETRA,ERROR,*)

C#### Subroutine: PEEL_AWAY
C###  Description:
C###    Peels away the boundary and initial tetrahedral. Genmesh.
CC JMB 18-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
!     Parameter List
      INTEGER ELEM(LDELEM,4),FACE(LDFACE,NFACES),IREGION(0:N_GM),
     '  TETRA(0:LDTETRA,4)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,REF
      LOGICAL BOOL

      CALL ENTERS('PEEL_AWAY',*9999)

      ! Peel away the boundary and first tetrahedrals
      REF=TETRA(0,HEAD)
      DO WHILE(REF.NE.0)
        ! Peel away tetrahedral containing any node associated with
        !  the first tetrahedral
        BOOL=.FALSE.
        DO I=1,NPTS
          BOOL=BOOL.OR.(IREGION(ELEM(REF,I)).EQ.0)
        ENDDO !i
        IF(BOOL)THEN
          CALL PEEL_TETRAHEDRAL(FACE,REF,TETRA,ERROR,*9999)
        ENDIF !bool
        IF(.NOT.BOOL)THEN
          ! Peel away any tetrahedral containing only boundary nodes
          BOOL=.TRUE.
          DO I=1,NPTS
            BOOL=BOOL.AND.(IREGION(ELEM(REF,I)).EQ.1)
          ENDDO !i
          IF(BOOL)THEN
            CALL PEEL_TETRAHEDRAL(FACE,REF,TETRA,ERROR,*9999)
          ENDIF !bool
        ENDIF !.not.bool
        REF=TETRA(REF,NEXT)
      ENDDO !while(ref.ne.0)
      CALL FIND_TETRAHEDRAL_TAIL(TETRA,ERROR,*9999)
      ! Peel away the first tetrahedral nodal data
      IREGION(0)=IREGION(0)-NPTS

      CALL EXITS('PEEL_AWAY')
      RETURN
 9999 CALL ERRORS('PEEL_AWAY',ERROR)
      CALL EXITS('PEEL_AWAY')
      RETURN 1
      END


