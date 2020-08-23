      SUBROUTINE ERRORS(ROUTINE_NAME,ERROR)

C#### Subroutine: ERRORS
C###  Description:
C###    ERRORS handles output of the message in ERROR, resets ERROR to
C###    the blank string, and handles storage and output of the call stack.
C###    Either ERRORS or ERRORIN should be called when exiting any
C###    routine after an error.  When called with an empty ROUTINENAME,
C###    the buffer containing the call stack from the error is flushed.

C KAT 2/6/00: Flushing ERROR string when it is full.
C     Also printing error message as soon as there is one.
C     This enables multiple error messages to be output is there is a
C     problem returning from an error.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERROR*(*),ROUTINE_NAME*(*)
!     Local Variables
!     Functions
      INTEGER LEN_TRIM

C KAT 2/6/00: I think ENTERS/EXITS are commented out because ERRORS is
C     called from some routines before ENTERS is suitably initialized.
C     Probably could be solved by having a common variable to say
C     whether ENTERS had been initialized or not.
C      CALL ENTERS('ERRORS',*9999)

C     In case FLAG_ERROR has not been called
      IF(ERROR.NE.' ') THEN
        CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
        ERROR = ' '
      ENDIF

      CALL ERRORIN(ROUTINE_NAME)

C KAT 2/6/00: Old
C      CALL STRING_TRIM(ERROR,IBEG,IEND)
C      ERROR=ERROR(IBEG:IEND)//'>'//ROUTINENAME

C     CALL EXITS('ERRORS')
      RETURN
      END


