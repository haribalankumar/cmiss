      SUBROUTINE SET_HANDLER(ERROR,*)

C#### Subroutine: SET_HANDLER
C###  Description:
C###    Sets the error handling routines from the old (default) routines

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables

      CALL ENTERS('SET_HANDLER',*9999)

      CALL SETFATALHANDLER()

      CALL EXITS('SET_HANDLER')
      RETURN
 9999 CALL ERRORS('SET_HANDLER',ERROR)
      CALL EXITS('SET_HANDLER')
      RETURN 1
      END


