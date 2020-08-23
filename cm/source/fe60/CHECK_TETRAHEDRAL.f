      SUBROUTINE CHECK_TETRAHEDRAL(FACE,FOUND,REF,TETRA,CIRC,NODE,
     '  RSQ,ERROR,*)

C#### Subroutine: CHECK_TETRAHEDRAL
C###  Description:
C###    Gemesh routine.
CC JMB 16-NOV-2001

      IMPLICIT NONE
      INCLUDE 'genmesh.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER FACE(LDFACE,NFACES),FOUND(0:N_GM),REF,TETRA(0:LDTETRA,4)
      REAL*8 CIRC(LDCIRC,3),NODE(LDNODE,3),RSQ(N_GM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I
      REAL*8 SUM
!     Local Functions
      REAL*8 DLAMCH
      LOGICAL FOUND_LIST

      CALL ENTERS('CHECK_TETRAHEDRAL',*9999)

      SUM=ZERO
      DO I=1,NJT
        SUM=SUM+(NODE(1,I)-CIRC(REF,I))**2.d0
      ENDDO !i

      IF(SUM.LT.(RSQ(REF)-DSQRT(DLAMCH('E'))))THEN
        IF(FOUND_LIST(FOUND,REF))THEN
          DO I=1,NFACES
            IF(FACE(REF,I).NE.0)THEN
              CALL CHECK_TETRAHEDRAL2(FACE,FOUND,FACE(REF,I),TETRA,
     '          CIRC,NODE,RSQ,ERROR,*9999)
            ENDIF !face
          ENDDO !i
        ENDIF !found_list
      ENDIF !sum

      CALL EXITS('CHECK_TETRAHEDRAL')
      RETURN
 9999 CALL ERRORS('CHECK_TETRAHEDRAL',ERROR)
      CALL EXITS('CHECK_TETRAHEDRAL')
      RETURN 1
      END


