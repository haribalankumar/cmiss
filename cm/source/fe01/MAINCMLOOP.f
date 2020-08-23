      SUBROUTINE MAINCMLOOP(BATCH_MODE,EXNAME,COMFNAME,
     '  PARAMETERSFILENAME,FIRST,EXECUTESTRING,ERR)

C#### Subroutine: MAINCMLOOP
C###  Description:
C###    Main CMISS Loop

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER (ROUTINENAME='MAINCMLOOP')
      INCLUDE 'fsklib.inc'
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbpr00.cmn'
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'gtstr00.cmn'
C       INCLUDE 'host00.cmn'
C       INCLUDE 'host00.inc'
      INCLUDE 'ktyp00.cmn'
!     Parameter List
      INTEGER BATCH_MODE,EXNAME(*),COMFNAME(*),
     '  FIRST,EXECUTESTRING(*),ERR
      CHARACTER PARAMETERSFILENAME(*) !actually C string
!     Local Variables
      INTEGER CODE,COMLEN,CURRENT_TIME,EXLEN,
     '  i,IEND,LENGTH,
     '  PRIMARY,SECONDARY,SLEEP_DELAY
      CHARACTER COMFILENAME*255,EXAMPLENUM*255,
     '  ERROR*(ERRSTRLEN),STRG*(MXCH),COMMAND_LINE*(MXCH),
     '  INTERPRETED_LINE(MXCH)
      LOGICAL CONTINUE,QUIT
!     Functions
      INTEGER C_STRLEN,LEN_TRIM

      CONTINUE=.TRUE.
      QUIT=.FALSE.
      NEST=0 !COMFILE_READ_LEVEL
      ERR=0

      IF(FIRST.NE.0) THEN

C       Initialise the cm-cmgui link
        IF(CMGUI_LINK) CALL CMGUI_LINK_INITIALISE(ERROR,*9999)

C KAT 28Jan00: Now done in PARSE
C        DO NOSG=1,MXSG
C          ISEG(NOSG)=0
C          CSEG(NOSG)=' '
C        ENDDO
CC!!! KATs 28Jan00: This should be removed when the new interpreter is used.
C        WRITE(STRG,'(A,1X,E11.5)') 'assign PI',PI
C        CALL PARSE(QUIT,STRG,ERROR,*9999)
CC         COD(2)='PI'
CC         WRITE(COD(3),'(E11.5)') PI
CC         CALL ASSIGN(COMMAND_LINE,1,COD,COQUD,ERROR,*9999)
CC!!! KATe 28Jan00

C     Subscript on INTERPRETED_LINE is to work around an internal
C     compiler error with xlf on AIX when bounds checking and
C     -qsmp=omp are specified.
        LENGTH=C_STRLEN(%REF(PARAMETERSFILENAME(1)))
        USEPARAMFILE=LENGTH.NE.0
        IF(USEPARAMFILE) THEN
          IEND=0
          CALL APPENDCA(IEND,MXCH,'fem define parameters;r;',
     '      INTERPRETED_LINE,ERR)
          IF(ERR.NE.0) THEN
            ERR=-1
            CALL FLAG_ERROR(ERR,'Buffer not large enough for command')
            GOTO 9998
          ENDIF
          IF(IEND+LENGTH.GT.MXCH) THEN
            LENGTH=MXCH-IEND
            IEND=0
            CALL APPENDC(IEND,
     '        'Parameters filename cannot be longer than ',ERROR)
            CALL APPENDI(IEND,LENGTH,ERROR)
            CALL APPENDC(IEND,' characters',ERROR)
            CALL FLAG_ERROR(1,ERROR(:IEND))
            GOTO 9998
          ENDIF
          DO i=1,LENGTH
            IEND=IEND+1
            INTERPRETED_LINE(IEND)=PARAMETERSFILENAME(i)
          ENDDO          
C         temporarily set FIRST params to .FALSE. so that arrays
C         won't be allocated in FEM until after params read in.
          FIRST_SYNTAX=.FALSE.
          CALL PARSE(QUIT,IEND,INTERPRETED_LINE,*9998)
          FIRST_SYNTAX=.TRUE. !reset FIRST params
        ENDIF
      ENDIF !FIRST

C KAT 28Jan00: Socket connection to cmgui no longer used
C      IF(.NOT.USE_SOCKET) THEN
      IF(CMGUI_LINK) THEN

        CONTINUE=.TRUE.
        DO WHILE (CONTINUE)
C         Perform any periodic tasks
          CALL PERIODICTASK
C         Check for available messages (rest in SYNTAX)
          CALL WH_INPUT_F_UPDATE(CMGUI_PROMPT_I,CODE)
          IF(CODE.EQ.0) THEN
            ERROR='>>Could not update prompt input'
            GOTO 9999
          ENDIF
          CALL WH_INPUT_F_UPDATE(CMGUI_COMMAND_I,CODE)
          IF(CODE.EQ.0) THEN
            ERROR='>>Could not update command input'
            GOTO 9999
          ENDIF
          CALL WH_OUTPUT_F_UPDATE(CMGUI_COMMAND_O,CODE)
          IF(CODE.EQ.0) THEN
            ERROR='>>Could not update command output'
            GOTO 9999
          ENDIF
          CALL WH_OUTPUT_F_CAN_OPEN(CMGUI_COMMAND_O,CODE)
          IF(CODE.NE.0) THEN
            CMGUI_MESSAGE_TIME=0
C           We have a command to process
            CALL WH_OUTPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_O,
     '        PRIMARY,SECONDARY,CODE)
            IF(CODE.EQ.0) THEN
              ERROR='>>Could not open command message'
              GOTO 9999
            ENDIF
            IF(SECONDARY.EQ.1) THEN !idle message - do nothing
            ELSEIF(SECONDARY.EQ.2) THEN !command
              CALL WH_OUTPUT_F_NUM_ITEMS(CMGUI_COMMAND_O,LENGTH)
              IF(LENGTH.LE.MXCH) THEN
CC               Zero the string
C                COMMAND_LINE=' '
                CALL WH_OUTPUT_F_GET_CHAR(CMGUI_COMMAND_O,
     '            LENGTH,INTERPRETED_LINE,CODE)
                IF(CODE.EQ.0) THEN
                  ERROR='>>Could not get command'
                  GOTO 9999
                ENDIF
C KAT 16Nov99:  This was already commented out.
C               Send acknowledgement of command receipt
C               CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_I,
C               '                0,4,CODE)
C               IF(CODE.EQ.0) THEN
C               ERROR='>>Could not open command receipt message'
C               GOTO 9999
C               ENDIF
C               CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_I,
C               '                CODE)
C               IF(CODE.EQ.0) THEN
C               ERROR='>>Could not close command receipt message'
C               GOTO 9999
C               ENDIF
                CALL PARSE(QUIT,LENGTH,INTERPRETED_LINE,*170)
                GOTO 220
C               Handle error condition
 170            CALL ERRORIN(' ') ! Flush error call stack
C               CALL OUTPUT_ERROR_STRING(1,ERROR)
C               CALL STRING_TRIM(ERROR,IBEG,IEND)
C               WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//
C               '              '>MAIN_LOOP'
C               CALL WRITES(IOER,OP_STRING,ERROR,*9998)
C               ERROR(1:)=' '
 220            CONTINUE
C KAT 16Nov99:  This was already commented out.
C               Check for quitting
                IF (QUIT) THEN !time to quit
                  CONTINUE=.FALSE.
                ELSE
C                 Send command completion
                  CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_I,0,5,CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR='>>Could not open command completion message'
                    GOTO 9999
                  ENDIF
                  CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_I,CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR='>>Could not close command completion message'
                    GOTO 9999
                  ENDIF
                  CALL WH_INPUT_F_UPDATE(CMGUI_COMMAND_I,CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR=
     '                '>>Could not update command input for completion'
                    GOTO 9999
                  ENDIF
                ENDIF
C               Tell CMGUI to pop down the prompt window if
C               necessary
                CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_PROMPT_I,0,3,CODE)
                IF(CODE.EQ.0) THEN
                  ERROR='>>Could not open prompt message'
                  GOTO 9999
                ENDIF
                CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_PROMPT_I,CODE)
                IF(CODE.EQ.0) THEN
                  ERROR='>>Could not close prompt message'
                  GOTO 9999
                ENDIF
              ELSE
                WRITE(OP_STRING,'('' Command is longer '
     '            //'than '',I5)')
     '            MXCH
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                CALL WH_OUTPUT_F_GET_REMAINDER(CMGUI_COMMAND_O,CODE)
                IF(CODE.EQ.0) THEN
                  ERROR='>>Could not use up command'
                  GOTO 9999
                ENDIF
              ENDIF
            ELSE
              ERROR='>>Invalid primary type'
              GOTO 9999
            ENDIF !primary
            CALL WH_OUTPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_O,CODE)
            IF(CODE.EQ.0) THEN
              ERROR='>>Could not close command message'
              GOTO 9999
            ENDIF
          ELSE !no command
C           CMGUI_MESSAGE_TIME is when we entered the idle loop
            IF(CMGUI_MESSAGE_TIME.EQ.0) THEN
              CALL GET_SECONDS(CMGUI_MESSAGE_TIME)
              IDLE_SENT=.FALSE.
            ELSE
C GMH 29/1/97 Ensure that we dont have cmiss running with no cmgui
              CALL GET_SECONDS(CURRENT_TIME)
              IF(CURRENT_TIME-CMGUI_MESSAGE_TIME.GT.
     '          CMGUI_IDLE_TIME) THEN
                IF(IDLE_SENT) THEN
                  IF(CURRENT_TIME-CMGUI_MESSAGE_TIME
     '              .GT.2*CMGUI_IDLE_TIME) THEN !time to quit
                    CONTINUE=.FALSE.
                  ENDIF
                ELSE
                  CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_I,
     '              0,1,CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR='>>Could not open idle message'
                    GOTO 9999
                  ENDIF
                  CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_I,
     '              CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR='>>Could not close idle message'
                    GOTO 9999
                  ENDIF
                  IDLE_SENT=.TRUE.
                ENDIF
              ENDIF !past idle
C             Try to send any information we may have in idle time
              CALL WH_INPUT_F_UPDATE(CMGUI_DATA_I,CODE)
              IF(CODE.EQ.0) THEN
                ERROR='>>Could not update data input'
                GOTO 9999
              ENDIF
              CALL WH_OUTPUT_F_UPDATE(CMGUI_DATA_O,CODE)
              IF(CODE.EQ.0) THEN
                ERROR='>>Could not update data output'
                GOTO 9999
              ENDIF
              SLEEP_DELAY=1 !seconds
              CALL SLEEPER(SLEEP_DELAY)
            ENDIF
          ENDIF

        ENDDO !main loop
C       Notify the front end by sending an empty message Secondary ID=3
        CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_I,0,3,CODE)
        IF(CODE.EQ.0) THEN
          ERROR='>>Could not open quit message'
          GOTO 9999
        ENDIF
        CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_I,CODE)
        IF(CODE.EQ.0) THEN
          ERROR='>>Could not close quit message'
          GOTO 9999
        ENDIF
        CALL WH_INPUT_F_UPDATE(CMGUI_COMMAND_I,CODE)
        IF(CODE.EQ.0) THEN
          ERROR='>>Could not update command input for quit message'
          GOTO 9999
        ENDIF

      ELSE !NOT cmgui_link

        IF(FIRST.NE.0) THEN
C         The fatal signal handler is not activated by default for
C         multithreaded versions as we don't want a non-master thread
C         trying to longjmp to a setjmp in the master thread.
C$        IF(.FALSE.) THEN
          IF(BATCH_MODE.EQ.0) THEN
C           Interactive: activate fatal signal handler to revert to user
C           input. Note that this is only activated after the
C           interpreter is initialized as we need an interpreter for
C           user input.
            CALL SET_HANDLER(ERROR,*9999)
          ENDIF
C$        ENDIF

C         SAB 10 Oct 2000
C         Execute the command line string
          LENGTH=C_STRLEN(EXECUTESTRING)
          IF(LENGTH.GT.0) THEN
            IF(LENGTH.GT.LEN(COMMAND_LINE)) THEN
              IEND=0
              CALL APPENDC(IEND,
     '          'execute string cannot be longer than ',ERROR)
              CALL APPENDI(IEND,LEN(COMMAND_LINE),ERROR)
              CALL APPENDC(IEND,' characters',ERROR)
              CALL FLAG_ERROR(0,ERROR(:IEND))
              GOTO 9998
            ENDIF
            CALL C2FSTRING(EXECUTESTRING,LENGTH,COMMAND_LINE)
            CALL STORE_COMMAND(COMMAND_LINE(:LENGTH),ERR,ERROR,*9999)
            CALL INTERPRET_COMMAND_LINE(QUIT,COMMAND_LINE(:LENGTH),ERR)
          ENDIF

          COMLEN=C_STRLEN(COMFNAME)
          EXLEN=C_STRLEN(EXNAME)

          IF(ERR.EQ.0.AND.(EXLEN.NE.0.OR.COMLEN.NE.0)) THEN
C           Get command file name
            IF(COMLEN.NE.0) CALL C2FSTRING(COMFNAME,COMLEN,COMFILENAME)
            IF(EXLEN.NE.0) THEN
C             Get example number
              CALL C2FSTRING(EXNAME,EXLEN,EXAMPLENUM)
C             set example directory
              IF(COMLEN.EQ.0) THEN
C               No comfilename supplied; get it from examplenum
                CALL SETEXAMPLEDIR(EXAMPLENUM(:EXLEN),COMFILENAME,
     '            ERROR,*9999)
                COMLEN=LEN_TRIM(COMFILENAME)
                IF(COMLEN.EQ.0) THEN
                  COMFILENAME='example_'//EXAMPLENUM(:EXLEN)
                  COMLEN=LEN_TRIM(COMFILENAME)
                ENDIF
              ELSE
C               STRG will not be used; use the comfilename supplied
                CALL SETEXAMPLEDIR(EXAMPLENUM(:EXLEN),STRG,
     '            ERROR,*9999)
              ENDIF
              COMMAND_LINE=
     '          'read command;'//COMFILENAME(:COMLEN)//';example'
            ELSE
              COMMAND_LINE='read command;'//COMFILENAME(:COMLEN)
            ENDIF
            IEND=LEN_TRIM(COMMAND_LINE)
            CALL STORE_COMMAND(COMMAND_LINE(:IEND),ERR,ERROR,*9999)
            CALL INTERPRET_COMMAND_LINE(QUIT,COMMAND_LINE(:IEND),ERR)

            IF(QUIT) THEN !quit issued from command file
              CONTINUE=.FALSE.
            ELSE
C             Initialize prompt
              PR(1:1)=CHAR(0)
            ENDIF
          ENDIF
        ENDIF

        IF(ERR.NE.0) THEN
          IF(BATCH_MODE.NE.0) THEN
            GO TO 9998
          ELSE
            CALL ERRORIN(' ') ! Flush error call stack
          ENDIF !batch mode
        ELSEIF(BATCH_MODE.NE.0) THEN
          CONTINUE=.FALSE.
        ENDIF

        DO WHILE(CONTINUE) !is main program loop

C GMH 15/11/95 making call conditional
          IF(USE_GRAPHICS.EQ.1) THEN
            CALL GXWAIT(0.0,ERR) !Update graphics
C KAT 24Feb99: reseting ERR.  Should something be done on error?
          ENDIF

          CALL GETSTR1(PR,COMMAND_LINE,ERROR,*9999)
          IEND=LEN_TRIM(COMMAND_LINE)
          IF(IEND.EQ.0) THEN !blank line
C           Reset prompt
            PR(1:1)=CHAR(0)
C           Fortran doesn't like empty character variables.
            IEND=1
          ENDIF
          CALL STORE_COMMAND(COMMAND_LINE(:IEND),ERR,ERROR,*150)
          CALL INTERPRET_COMMAND_LINE(QUIT,COMMAND_LINE(:IEND),ERR)
          IF(ERR.NE.0) GOTO 150
C         CALL PARSE(0,ISEG,CSEG,COMMAND_LINE,QUIT,ERROR,*150)
C         COMMAND_LINE=' '

          IF(QUIT) THEN !The QUIT command has been used to quit
            CONTINUE=.FALSE.
          ENDIF

          GOTO 200
C***      Handle error condition
 150      CALL ERRORIN(' ') ! Flush error call stack
C          CALL OUTPUT_ERROR_STRING(ERR,ERROR)
 200      CONTINUE
        ENDDO !end of main program loop
      ENDIF !cmgui_link

C KAT 28Jan00: Socket connection to cmgui no longer used
C      ELSE !USE_SOCKET = .TRUE.

C        CONNID1=0
C        CONNID2=1

C        IF (FSKLISTEN(CONNID1,PORT1) .EQ. -1) GOTO 9998
C        IF (FSKLISTEN(CONNID2,PORT2) .EQ. -1) GOTO 9998

C        FIRST=1

C        USEPARAMFILE=PARAMETERS.EQ.1
C        IF(USEPARAMFILE) THEN
C          CALL CSTRINGLEN(CSTRLEN,PARAMETERSFILENAME)
C          CALL C2FSTRING(PARAMETERSFILENAME,CSTRLEN,PARAMETERSFILE)
C          PARAMETERSFILE(CSTRLEN+1:)=' '
C          CALL STRING_TRIM(PARAMETERSFILE,IBEG,IEND)
C          STRG='fem define parameters;r;'//PARAMETERSFILE(IBEG:IEND)
CC         temporarily set FIRST params to .FALSE. so that arrays won't
CC         be allocated in FEM until after params read in.
C          FIRST_SYNTAX=.FALSE.
C          CALL PARSE(0,ISEG,CSEG,STRG,END,ERROR,*9999)
C          FIRST_SYNTAX=.TRUE. !reset FIRST params
C        ENDIF

C        CONTINUE=.TRUE.
C        DO WHILE (CONTINUE)
CC CPB 6/4/95 Adding timeout for sockets to enable periodic tasks
C          RETURNVAL=0
C          DO WHILE(RETURNVAL.EQ.0)
C            RETURNVAL=FSKSELECT(CONNID1,200) ! timeout after 200ms
C            IF(RETURNVAL.EQ.-1) GOTO 9998
C            CALL PERIODICTASK
C    ENDDO
C          IF (FSKREAD(LEN,SK_LONG_INT,1,CONNID1) .EQ. -1) GOTO 9998
C          IF (FSKREAD(INTSTR,SK_CHAR,LEN+1,CONNID1) .EQ. -1) GOTO 9998
C          CALL FSKC2F(STRING,LEN,INTSTR)

C          CALL STRING_TRIM(STRING,IBEG,IEND)
C          WRITE(OP_STRING,'(A)') STRING(IBEG:IEND)
C          CALL WRITES(IOOP,OP_STRING,ERROR,*9998)
C          IF(STRING(IBEG:IBEG+3).NE.'QUIT') THEN
C            CALL PARSE(0,ISEG,CSEG,STRING,END,ERROR,*160)
C            WRITE(OP_STRING,'(/1X,A)') '> '
C            CALL WRITES(IOOP,OP_STRING,ERROR,*9998)
C          ENDIF

C          GOTO 210

CC***      Handle error condition
C 160      CALL STRING_TRIM(ERROR,IBEG,IEND)
C          WRITE(OP_STRING,'(A)') ' '//ERROR(IBEG:IEND)//'>MAIN_LOOP'
C          CALL WRITES(IOER,OP_STRING,ERROR,*9998)
C          ERROR(1:)=' '

CC***      Tell front end that CMISS is back in the main loop (used for
CC***      popping down the prompt window)
C 210      IF(FSKWRITE(-1,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9998

C          IF (STRING(1:4) .EQ. 'QUIT') THEN
CC***        First tell front-end to close down sockets, then close down
CC***        own sockets.
C            IF (FSKWRITE(-2,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9998
CC LKC 1-JUL-1999 unused?
CC            IRET=FSKCLOSE(CONNID1)
CC            IRET=FSKCLOSE(CONNID2)
C            CONTINUE=.FALSE.
C          ENDIF

C        ENDDO !main loop

C      ENDIF

      ERR=0

      GOTO 10000
 9998 ERROR=' '
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL ERRORIN(' ') ! Flush error call stack
      IF(BATCH_MODE.EQ.0) THEN
        CALL WRITE_CHAR(IOER,
     &    '>>Fatal Error in Main Loop'//NEWLINE,ERR)
      ENDIF
      IF(ERR.EQ.0) ERR=1

10000 CONTINUE

C***  Do finishing things.

C C ??? MPN 27-oct-94 There must be a better place for this?
C C MPN 13-Sep-94: Parallel element stiffness matrix computations
C       IF(KTYP1A.EQ.2) THEN  !parallel element stiffness matrix calcs
C         DO nhost=1,NUMHOSTS_USED
C           IF(SOCKET_OPEN(nhost)) THEN
C C           Signal to slave processes to stop
C             IF(FSKWRITE(QUIT_PROCESS,SK_LONG_INT,1,ICONNID(nhost))
C      '        .EQ.-1) CALL FLAG_ERROR(0,'Failed to signal slave')
C C LKC 1-JUL-1999 unused?
C C            IRET=FSKCLOSE(ICONNID(nhost)) !close socket
C           ENDIF
C         ENDDO
C       ENDIF

C     On AIX open files are not closed (flushed) automatically on exit
C     as we don't have a fortran main program so we must explicitly
C     close open files.
      IF(ECHO_OUTPUT) THEN
        ECHO_OUTPUT=.FALSE.
        CALL CLOSEF(IOOUT,ERROR,*10099)
        GOTO 10100
10099     CALL ERRORS(ROUTINENAME,ERROR)
          CALL ERRORIN(' ') ! Flush error call stack
          IF(ERR.EQ.0) ERR=1
10100   CONTINUE
      ENDIF
      RETURN
      END


