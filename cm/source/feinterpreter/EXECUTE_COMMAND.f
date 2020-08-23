      SUBROUTINE EXECUTE_COMMAND(INTERPRETED_LINE,USER_DATA,QUIT,STATUS)

C#### Subroutine: EXECUTE_COMMAND
C###  Description:
C###    Called from the interpreter to execute a command after
C###    control statements, variables, and expressions have been
C###    interpreted.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='EXECUTE_COMMAND')
      INCLUDE 'b00.cmn'

      !Argument Variables
      CHARACTER INTERPRETED_LINE(*) !actually C string
      INTEGER USER_DATA !unused
      LOGICAL QUIT
      INTEGER STATUS
!     Functions
      INTEGER C_STRLEN

      CALL ENTERS(ROUTINENAME,*999)

C     Subscript on INTERPRETED_LINE is to work around an internal
C     compiler error with xlf on AIX when bounds checking and
C     -qsmp=omp:noopt are specified.
      CALL PARSE(QUIT,C_STRLEN(%REF(INTERPRETED_LINE(1))),
     '  INTERPRETED_LINE,*999)

      STATUS=1
      CALL EXITS(ROUTINENAME)
      RETURN
 999  CALL ERRORIN(ROUTINENAME)
      CALL ERRORIN(' ') !flush call stack
      STATUS=0
      CALL EXITS(ROUTINENAME)
      RETURN
      END


