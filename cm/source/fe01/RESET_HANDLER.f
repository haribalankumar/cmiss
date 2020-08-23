      SUBROUTINE RESET_HANDLER(ERROR,*)

C#### Subroutine: RESET_HANDLER
C###  Description:
C###    Resets the error handler to the old (default) error handling
C###    routines.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('RESET_HANDLER',*9999)

      CALL RESETFATALHANDLER()

      CALL EXITS('RESET_HANDLER')
      RETURN
 9999 CALL ERRORS('RESET_HANDLER',ERROR)
      CALL EXITS('RESET_HANDLER')
      RETURN 1
      END


