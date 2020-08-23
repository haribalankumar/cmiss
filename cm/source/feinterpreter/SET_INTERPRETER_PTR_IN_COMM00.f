      SUBROUTINE SET_INTERPRETER_PTR_IN_COMM00(VALUE)

C#### Subroutine: SET_INTERPRETER_PTR_IN_COMM00
C###  Description:
C###    Stores the interpreter pointer in the common block so it can
C###    be passed to all the interpreter commands.
C     Created: SAB 2004-01-26

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='SET_INTERPRETER_PTR_IN_COMM00')

      INCLUDE 'comm00.cmn'

      !Argument Variables
      POINTER VALUE

      CALL ENTERS(ROUTINENAME,*9999)

      INTERPRETER_PTR = VALUE

 9999 CALL EXITS(ROUTINENAME)
      END ! SUBROUTINE SET_INTERPRETER_PTR_IN_COMM00


