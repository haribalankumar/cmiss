      SUBROUTINE INTERPRET_COMMAND_LINE(QUIT,COMMAND_LINE,ERR_CODE)

C#### Subroutine: INTERPRET_COMMAND_LINE
C###  Description:
C###    Calls the interpreter to interpret and execute a command line.
C     Created: KAT 6/2/00

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='INTERPRET_COMMAND_LINE')

      INCLUDE 'comm00.cmn'
C      INCLUDE 'cmiss$reference:mxch.inc'
C      INCLUDE 'cmiss$reference:cbpr00.cmn'
C      INCLUDE 'cmiss$reference:gtstr00.cmn'

!     Parameter List
      LOGICAL QUIT
      CHARACTER*(*) COMMAND_LINE
      INTEGER ERR_CODE
!     Local Variables
      INTEGER*4 COMMAND_PTR
      INTEGER STATUS
CC DBs 25/3/01.  F77 does not allow functions in parameter statements
C      CHARACTER SEPARATOR
CC DBe  PARAMETER (SEPARATOR=ACHAR(0))
C      INTEGER ERROR_LENGTH,LINE_BEG,LINE_END,POSITION
C      CHARACTER ERROR_STRING*(ERRSTRLEN)

!     Functions
      EXTERNAL EXECUTE_COMMAND
C      INTEGER LEN_TRIM1


CC DBs 25/3/01.  F77 does not allow functions in parameter statements
C      SEPARATOR=ACHAR(0)
CC DBe
      CALL ENTERS(ROUTINENAME,*1)
 1    CONTINUE

      COMMAND_PTR=0
      CALL CREATE_CSTRING(COMMAND_PTR,COMMAND_LINE)
      IF(COMMAND_PTR.EQ.0) GOTO 999

      CALL INTERPRET_COMMAND(%VAL(INTERPRETER_PTR),%VAL(COMMAND_PTR),
     &  %VAL(0),QUIT,EXECUTE_COMMAND,STATUS)

      CALL DESTROY_CSTRING(COMMAND_PTR)

      IF(STATUS.EQ.0) GOTO 999

C        IF(ERROR_STRING.NE.' ') THEN
CC         This is for outputing error messages and call stack from
CC         CPB interpreter.
C          ERROR_LENGTH=LEN_TRIM1(ERROR_STRING)
C          POSITION = INDEX(ERROR_STRING(:ERROR_LENGTH),SEPARATOR)
C          IF(POSITION.EQ.0) THEN
C            LINE_END = ERROR_LENGTH
C          ELSE
C            LINE_END = POSITION - 1
C          ENDIF
C          IF(POSITION.GT.1) THEN
C            CALL FLAG_ERROR(ERR_CODE,ERROR_STRING(:LINE_END))
C          ENDIF
C          DO WHILE(LINE_END.LT.ERROR_LENGTH)
C            LINE_BEG = LINE_END + 2
C            POSITION=INDEX(ERROR_STRING(LINE_BEG:ERROR_LENGTH),CHAR(0))
C            IF(POSITION.EQ.0) THEN
C              LINE_END = ERROR_LENGTH
C            ELSE
C              LINE_END=LINE_BEG+POSITION-2
C            ENDIF
C            CALL ERRORIN(ERROR_STRING(LINE_BEG:LINE_END))
C          ENDDO
C        ENDIF

      ERR_CODE=0
      CALL EXITS(ROUTINENAME)
      RETURN

 999  CALL ERRORIN(ROUTINENAME)
      ERR_CODE=1
      CALL EXITS(ROUTINENAME)
      END


