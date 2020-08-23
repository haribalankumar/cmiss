      SUBROUTINE ERRORIN(ROUTINE_NAME)

C#### Subroutine: ERRORIN
C###  Description:
C###    ERRORIN handles storage and output of the call stack after an
C###    error has occurred.  Either ERRORS or ERRORIN should be called
C###    when exiting any routine after an error.  When called with an
C###    empty ROUTINENAME, the call stack is completed with a newline.
C###    The intention is that a routine should call FLAG_ERROR with an
C###    error message and error code, then call this routine on exit
C###    from each routine, then the main loop calls this routine with an
C###    empty ROUTINENAME to complete the call stack.
C###  See-Also: ERRORS

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      CHARACTER ROUTINE_NAME*(*)
!     Local Variables
      CHARACTER*2 PREFIX,SEPARATOR
      PARAMETER (PREFIX='  ',SEPARATOR=' >')
      INTEGER ERR
!     Static Variables
      ! Try to format the call stack so that routine names are not
      ! wrapped, so keep a record of how much is on the current line.
      INTEGER LINE_LEN
      DATA LINE_LEN/0/

C KAT 2/6/00: I think ENTERS/EXITS are commented out because ERRORS is
C     called from some routines before ENTERS is suitably initialized.
C     Probably could be solved by having a common variable to say
C     whether ENTERS had been initialized or not.
C      CALL ENTERS('ERRORIN',*9999)

      IF(ROUTINE_NAME.EQ.' ') THEN
C       Add trailing newline if necessary
        IF(LINE_LEN.NE.0) THEN
          CALL WRITE_CHAR(IOER,NEWLINE,ERR)
          LINE_LEN=0
        ENDIF
      ELSE
        IF(LINE_LEN.NE.0.AND.
     &    LINE_LEN+LEN(SEPARATOR)+LEN(ROUTINE_NAME).GT.79) THEN
C         Start a new line
          CALL WRITE_CHAR(IOER,NEWLINE,ERR)
          LINE_LEN=0
        ENDIF
        IF(LINE_LEN.EQ.0) THEN
          CALL WRITE_CHAR(IOER,PREFIX,ERR)
          LINE_LEN=LEN(PREFIX)
        ENDIF
        CALL WRITE_CHAR(IOER,SEPARATOR,ERR)
        CALL WRITE_CHAR(IOER,ROUTINE_NAME,ERR)
        LINE_LEN=LINE_LEN+LEN(SEPARATOR)+LEN(ROUTINE_NAME)
      ENDIF

C     CALL EXITS('ERRORIN')
      RETURN
      END


