      SUBROUTINE FLAG_ERROR(ERR_CODE,MESSAGE)

C#### Subroutine: FLAG_ERROR
C###  Description:
C###    Standard routine to indicate an error, displaying MESSAGE to the
C###    user.  If MESSAGE is blank, then no newline is output so that a
C###    message may be added.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER ERR_CODE
      CHARACTER*(*) MESSAGE
!     Local Variables
      INTEGER ERR

C KAT 2/6/00: I think ENTERS/EXITS are commented out because ERRORS is
C     called from some routines before ENTERS is suitably initialized.
C     Probably could be solved by having a common variable to say
C     whether ENTERS had been initialized or not.

      CALL ERRORIN(' ') ! Flush error call stack

      CALL WRITE_CHAR(IOER,'>>ERROR: ',ERR)
      IF(ERR_CODE.NE.0) THEN
        CALL WRITE_INT(IOER,ERR_CODE,ERR)
        CALL WRITE_CHAR(IOER,': ',ERR)
      ENDIF
      IF(MESSAGE.NE.' ') THEN
        CALL WRITE_STRING(IOER,LEN(MESSAGE),MESSAGE,ERR)
        CALL WRITE_STRING(IOER,1,NEWLINE,ERR)
      ENDIF

      RETURN
      END


