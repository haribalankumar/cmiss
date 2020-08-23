      SUBROUTINE ASSERT_VORO(LOEXPR,ERR,ERROR,*)

C#### Subroutine: ASSERT_VORO
C###  Description:
C###    Error handling for a Genmesh routine.
CC JMB 18-NOV-2001

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERR*(*)
      LOGICAL LOEXPR
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I
!     Local Functions
      INTEGER LNBLNK

      CALL ENTERS('ASSERT_VORO',*9999)

      IF(.NOT.LOEXPR)THEN
        ! Terminate the program
        I=LNBLNK(ERR)
        WRITE(0,*) '>> ERROR : ', ERR(1:I)
        STOP
      ENDIF !.not.loexpr

      CALL EXITS('ASSERT_VORO')
      RETURN
 9999 CALL ERRORS('ASSERT_VORO',ERROR)
      CALL EXITS('ASSERT_VORO')
      RETURN 1
      END


