      SUBROUTINE SET_USER_INTEGER(NAME,VALUE,ERR_CODE)

C#### Subroutine: SET_USER_INTEGER
C###  Description:
C###    Sets an interpreter variable NAME to integer VALUE.
C     Created: KAT 2001-09-19

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='SET_USER_INTEGER')

      INCLUDE 'comm00.cmn'

      !Argument Variables
      CHARACTER*(*) NAME
      INTEGER VALUE
      INTEGER ERR_CODE
      !Local Variables
      INTEGER*4 NAME_PTR
      INTEGER STATUS

      CALL ENTERS(ROUTINENAME,*999)

      NAME_PTR=0
      CALL CREATE_CSTRING(NAME_PTR,NAME)
      IF(NAME_PTR.EQ.0) GOTO 999

      CALL INTERPRETER_SET_INTEGER(%VAL(INTERPRETER_PTR),
     &  %VAL(NAME_PTR),VALUE,STATUS)

      CALL DESTROY_CSTRING(NAME_PTR)

      IF(STATUS.EQ.0) GOTO 999

      ERR_CODE=0
      CALL EXITS(ROUTINENAME)
      RETURN

 999  CALL ERRORIN(ROUTINENAME)
      ERR_CODE=1
      CALL EXITS(ROUTINENAME)
      END ! SUBROUTINE SET_USER_INTEGER


C KAT 2001-09-18: Code only required with CPB interpreter

C###  Routine: INPUT_INTERPRETER_INPUT       prompted input from user
C###  Routine: OUTPUT_INTERPRETER_ERROR      error output
C###  Routine: OUTPUT_INTERPRETER_OUTPUT     general output
C###  Routine: OUTPUT_INTERPRETER_WARNING    warning output
C###  Routine: PREPROCESS_COMMAND_LINE       output / store command

C      SUBROUTINE INPUT_INTERPRETER_INPUT(PROMPT,RETURNVALUE,ERR,ERROR)

CC#### Subroutine: INPUT_INTERPRETER_INPUT
CC###  Description:
CC###    Inputs a string for the interpreter.

C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:cbdi02.cmn'

C      !Argument Variables
C      CHARACTER*(*) PROMPT,RETURNVALUE
C      INTEGER ERR
C      CHARACTER ERROR*(*)

C      CALL ENTERS('INPUT_INTERPRETER_INPUT',*998)

C      IF(PROMPT.NE.' ') THEN
C        OP_STRING(1) = PROMPT//' ?'
C        CALL WRITES(IOOP,OP_STRING,ERROR,*998)
C      ENDIF
C      READ(*,'(A)',IOSTAT=ERR) RETURNVALUE
C      IF(ERR.NE.0) GOTO 999

C      ERR=0
C      CALL EXITS('INPUT_INTERPRETER_INPUT')
C      RETURN
C 998  ERR=1
C 999  CALL ERRORS('INPUT_INTERPRETER_INPUT',ERROR)
C      CALL EXITS('INPUT_INTERPRETER_INPUT')
C      RETURN
C      END


C      SUBROUTINE OUTPUT_INTERPRETER_ERROR(POSITION,COMMAND_LINE,ERR,
C     '  ERROR)

CC#### Subroutine: OUTPUT_INTERPRETER_ERROR
CC###  Description:
CC###    Indicates the position of an interpreter error.

CC***  Extended KAT 20/3/00 from CPB F90 version.

C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:gtstr00.cmn'

C      INTEGER OUTPUTWIDTH !should be in a global include file
C      PARAMETER(OUTPUTWIDTH=80) !LEN(BLANK)
C      INTEGER ERRCOLGOAL
C      PARAMETER(ERRCOLGOAL=OUTPUTWIDTH/2+2) !fairly arbitrary but near middle
C      !Argument Variables
C      INTEGER POSITION
C      CHARACTER*(*) COMMAND_LINE
C      INTEGER ERR
C      CHARACTER ERROR*(*)
C      !Local Variables
C      INTEGER ECHOBEG,ECHOEND,IEND,NSPACES

C      CALL ENTERS('OUTPUT_INTERPRETER_ERROR',*999)

CC      NSPACES=POSITION+1 !position + len('> ') - 1
CC      IF(NEST.GT.0) NSPACES=NSPACES+1 !extra space before prompt in readcom
CC      IF(POSITION.LE.LEN(BLANK)) THEN
CC        OP_STRING(1)=BLANK(:NSPACES)//'^'
CC        CALL WRITES(IOER,OP_STRING,ERROR,*999)
CC      ENDIF

CC     Try to output the line aligned with its original state
C      NSPACES=2 ! len('> ')
C      IF(NEST.GT.0) NSPACES=NSPACES+1 !extra space before prompt in readcom

CC     Ouput the whole line if possible
C      ECHOBEG=1
C      ECHOEND=LEN(COMMAND_LINE) !could trim comments if keen

C      IF(ECHOEND+NSPACES.LE.OUTPUTWIDTH) THEN
C      ELSE
C        NSPACES=1
C        IF(ECHOEND+NSPACES.LE.OUTPUTWIDTH) THEN
C        ELSE IF(POSITION+NSPACES.LE.ERRCOLGOAL) THEN
CC         Ouput the beginning of the line if the hat will end up to the left
CC         of the goal column
C          ECHOEND=OUTPUTWIDTH-NSPACES
C        ELSE IF(ECHOEND-POSITION.LE.OUTPUTWIDTH-ERRCOLGOAL) THEN
CC         Ouput the end of the line if the hat will end up to the right
CC         of the goal column
C          ECHOBEG=ECHOEND-OUTPUTWIDTH+NSPACES+1
C        ELSE
CC         Put the hat in the goal column
C          ECHOBEG=POSITION-ERRCOLGOAL+NSPACES+1
C          ECHOEND=OUTPUTWIDTH-NSPACES+ECHOBEG-1
C        ENDIF
C      ENDIF
C      OP_STRING(1)=BLANK(:NSPACES)//COMMAND_LINE(ECHOBEG:POSITION)
C      IEND=0
C      CALL APPENDC(IEND,
C     '  BLANK(:POSITION-ECHOBEG+NSPACES)//'^',OP_STRING(2))
C      IF(ECHOEND.GT.POSITION) CALL APPENDC(IEND,
C     '  COMMAND_LINE(POSITION+1:ECHOEND),OP_STRING(2))
C      CALL WRITES(IOER,OP_STRING,ERROR,*999)

C      ERR=0
C      CALL EXITS('OUTPUT_INTERPRETER_ERROR')
C      RETURN
C 999  ERR=1
C      CALL ERRORS('OUTPUT_INTERPRETER_ERROR',ERROR)
C      CALL EXITS('OUTPUT_INTERPRETER_ERROR')
C      RETURN
C      END


C      SUBROUTINE OUTPUT_INTERPRETER_OUTPUT(STRING,ERR,ERROR)

CC#### Subroutine: OUTPUT_INTERPRETER_OUTPUT
CC###  Description:
CC###    Outputs a string of general output from the interpreter.

C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:cbdi02.cmn'

C      !Argument Variables
C      CHARACTER*(*) STRING
C      INTEGER ERR
C      CHARACTER ERROR*(*)

C      CALL ENTERS('OUTPUT_INTERPRETER_OUTPUT',*999)

C      OP_STRING(1) = STRING
C      CALL WRITES(IOOP,OP_STRING,ERROR,*999)

C      ERR=0
C      CALL EXITS('OUTPUT_INTERPRETER_OUTPUT')
C      RETURN
C 999  ERR=1
C      CALL ERRORS('OUTPUT_INTERPRETER_OUTPUT',ERROR)
C      CALL EXITS('OUTPUT_INTERPRETER_OUTPUT')
C      RETURN
C      END


C      SUBROUTINE OUTPUT_INTERPRETER_WARNING(WARNING,ERR,ERROR)

CC#### Subroutine: OUTPUT_INTERPRETER_WARNING
CC###  Description:
CC###    Outputs an interpreter warning.

C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:cbdi02.cmn'

C      !Argument Variables
C      CHARACTER*(*) WARNING
C      INTEGER ERR
C      CHARACTER ERROR*(*)

C      CALL ENTERS('OUTPUT_INTERPRETER_WARNING',*999)

C      OP_STRING(1)='>>WARNING: '//WARNING
C      CALL WRITES(IOER,OP_STRING,ERROR,*999)

C      ERR=0
C      CALL EXITS('OUTPUT_INTERPRETER_WARNING')
C      RETURN
C 999  ERR=1
C      CALL ERRORS('OUTPUT_INTERPRETER_WARNING',ERROR)
C      CALL EXITS('OUTPUT_INTERPRETER_WARNING')
C      RETURN
C      END

C      SUBROUTINE PREPROCESS_COMMAND_LINE(INFUNCTION,INLOOP,INREAD,
C     '  COMMAND_LINE,ERR,ERROR)

CC#### Subroutine: PREPROCESS_COMMAND_LINE
CC###  Description:
CC###    Implementation dependent command preprocessing e.g. adds the
CC###    given command to any implementation dependent lists (e.g.
CC###    learn command files etc.). INLOOP is true if the command
CC###    originates inside a "loop" (do or while) and INREAD is true
CC###    if the command originates inside a comfile read.

CC***  Adapted KAT 2Dec99 from CPB F90 version.

C      IMPLICIT NONE

C!     Parameter List
C      LOGICAL INFUNCTION,INLOOP,INREAD
C      CHARACTER*(*) COMMAND_LINE
C      INTEGER ERR
C      CHARACTER ERROR*(*)
CC!     Local Variables
CC      INTEGER COMLINE_END,POSITION
CC      CHARACTER COMMENT*1
CC      PARAMETER (COMMENT='#')

CC      CALL ENTER(ROUTINENAME)
C      CALL ENTERS('PREPROCESS_COMMAND_LINE',*1)
C 1    CONTINUE

CC     Just let the command happen if it is from a loop
CC      IF(.NOT.(INLOOP.OR.INFUNCTION)) THEN
CC        IF(INREAD) THEN
CC KAT 28/3/00: Done in READCOM
CC          !Echo command if being read from a file
CC          COMLINE_END=LEN(COMMAND_LINE)
CC          POSITION=INDEX(COMMAND_LINE,COMMENT)
CC          IF(POSITION.NE.0.AND.POSITION.LT.COMLINE_END) THEN
CC            IF(COMMAND_LINE(POSITION+1:POSITION+1).EQ.COMMENT)
CC     '        COMLINE_END=POSITION-1
CC          ENDIF
CC          OP_STRING(1)=' > '//COMMAND_LINE(:COMLINE_END)
CC          CALL WRITES(IOOP,OP_STRING,ERROR,*999)
CC        ELSE
CC          !Add command to command buffer
CC          CALL ADDSTRTOBUFF(COMMAND_LINE,ERROR,*999)
CC          !Send command to learn file if appropriate
CC          IF(LEARN) THEN
CC            WRITE(UNIT=COM_UNIT,FMT='(A)',IOSTAT=ERR) COMMAND_LINE
CC            IF(ERR.NE.0) THEN
CC              ERROR='Error writing to learn file'
CC              GOTO 999
CCC              CALL SET_IOSTAT_ERROR(ERR,ERRORSTRING)
CCC              CALL FLAG_ERROR("Error writing learn file, "//&
CCC              & ERRORSTRING(1:LEN_TRIM(ERRORSTRING)),ERR,ERROR,*999)
CC            ENDIF !ERR
CC          ENDIF !LEARN
CC        ENDIF !INREAD
CC      ENDIF !not INLOOP

C      ERR=0
C      CALL EXITS('PREPROCESS_COMMAND_LINE')
C      RETURN
CC 999  ERR=1
CC      CALL ERRORS('PREPROCESS_COMMAND_LINE',ERROR)
CC      CALL EXITS('PREPROCESS_COMMAND_LINE')
CC      RETURN
C      END


