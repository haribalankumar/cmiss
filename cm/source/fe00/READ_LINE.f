      SUBROUTINE READ_LINE(IUNIT,IOSTAT,NCHAR,CLINE,ERROR,*)

C#### Subroutine: READ_LINE
C###  Description:
C###    READ_LINE reads CLINE from IUNIT, and returns number of
C###    characters NCHAR.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'cmis00.cmn'
!old  INCLUDE '[-.sklib]fsklib.inc'
      INCLUDE 'fsklib.inc'
!     Parameter List
      INTEGER IOSTAT,IUNIT,NCHAR
      CHARACTER CLINE*(*),ERROR*(*)
!     Local Variables
      INTEGER CLEN,CODE,CURRENT_TIME,IBEG,IEND,INTSTR(1024),
     '  PRIMARY,SECONDARY,SLEEP_DELAY
      LOGICAL UNFOUND

      CALL ENTERS('READ_LINE',*9999)
      IF(.NOT.USE_SOCKET) THEN
        IF(CMGUI_LINK) THEN
          IOSTAT=0 !No errors
          UNFOUND=.TRUE.
          DO WHILE (UNFOUND)
C         First update the prompt wormholes
            CALL WH_INPUT_F_UPDATE(CMGUI_PROMPT_I,CODE)
            IF(CODE.EQ.0) THEN
              ERROR='Could not update prompt input'
              GOTO 9999
            ENDIF
            CALL WH_OUTPUT_F_UPDATE(CMGUI_PROMPT_O,CODE)
            IF(CODE.EQ.0) THEN
              ERROR='Could not update prompt output'
              GOTO 9999
            ENDIF
            CALL WH_OUTPUT_F_CAN_OPEN(CMGUI_PROMPT_O,
     '        CODE)
            IF(CODE.NE.0) THEN
              CMGUI_MESSAGE_TIME=0
C             We have an input to process
              CALL WH_OUTPUT_F_OPEN_MESSAGE(CMGUI_PROMPT_O,
     '          PRIMARY,SECONDARY,CODE)
              IF(CODE.EQ.0) THEN
                ERROR='Could not open prompt message'
                GOTO 9999
              ENDIF
              IF(SECONDARY.EQ.1) THEN !idle_message  - do nothing
              ELSEIF(SECONDARY.EQ.2) THEN !response
                CALL WH_OUTPUT_F_NUM_ITEMS(CMGUI_PROMPT_O,
     '            NCHAR)
                IF(NCHAR.LE.LEN(CLINE)) THEN
C                 Zero the string
                  CLINE=' '
                  CALL WH_OUTPUT_F_GET_CHAR(CMGUI_PROMPT_O,
     '              NCHAR,CLINE,CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR='Could not get prompt'
                    GOTO 9999
                  ENDIF
C                 Check for a default response
                  IF(CLINE(1:9).EQ.'<default>') THEN
                    NCHAR=0
                    CLINE=' '
                  ENDIF
                ELSE
                  NCHAR=0
                  CLINE=' '
                  WRITE(OP_STRING,'(''Command is longer than'',I5)')
     '              LEN(CLINE)
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  CALL WH_OUTPUT_F_GET_REMAINDER(CMGUI_PROMPT_O,
     '              CODE)
                  IF(CODE.EQ.0) THEN
                    ERROR='Could not use up prompt'
                    GOTO 9999
                  ENDIF
                ENDIF
                UNFOUND=.FALSE.
              ELSE
                ERROR='Invalid secondary type'
                GOTO 9999
              ENDIF !primary
              CALL WH_OUTPUT_F_CLOSE_MESSAGE(CMGUI_PROMPT_O,
     '          CODE)
              IF(CODE.EQ.0) THEN
                ERROR='Could not close command message'
                GOTO 9999
              ENDIF
            ELSE !no input
C             CMGUI_MESSAGE_TIME is when we entered the idle loop
              IF(CMGUI_MESSAGE_TIME.EQ.0) THEN
                CALL GET_SECONDS(CMGUI_MESSAGE_TIME)
                IDLE_SENT=.FALSE.
              ELSE
C GMH 29/1/97 Ensure that we dont have cmiss running with no cmgui
                CALL GET_SECONDS(CURRENT_TIME)
                IF(CURRENT_TIME-CMGUI_MESSAGE_TIME.GT.
     '            CMGUI_IDLE_TIME) THEN
                  IF(IDLE_SENT) THEN
                    IF(CURRENT_TIME-CMGUI_MESSAGE_TIME
     '                .GT.2*CMGUI_IDLE_TIME) THEN !fake a 'r' response
                      CLINE='r'
                      NCHAR=1
                      UNFOUND=.FALSE.
                    ENDIF
                  ELSE
                    CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_I,
     '                0,1,CODE)
                    IF(CODE.EQ.0) THEN
                      ERROR='Could not open idle message'
                      GOTO 9999
                    ENDIF
                    CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_I,
     '                CODE)
                    IF(CODE.EQ.0) THEN
                      ERROR='Could not close idle message'
                      GOTO 9999
                    ENDIF
                    IDLE_SENT=.TRUE.
                  ENDIF
                ENDIF !past idle
                SLEEP_DELAY=1 !seconds
                CALL SLEEPER(SLEEP_DELAY)
              ENDIF
            ENDIF
          ENDDO !not found
        ELSE
C         READ(UNIT=IUNIT,FMT='(Q,A)',IOSTAT=IOSTAT) NCHAR,
C    '      CLINE(1:1+NCHAR)
          CLINE=' '
CGMH 15-May-95 READ(UNIT=IUNIT,FMT='(A)',IOSTAT=IOSTAT) CLINE(1:)
          READ(UNIT=IUNIT,FMT='(A)',IOSTAT=IOSTAT) CLINE
          CALL STRING_TRIM(CLINE,IBEG,IEND)
C MPN 30Apr97: bug
C old     IF(ICHAR(CLINE(1:1)).EQ.32) THEN !Return key pressed only
          IF(ICHAR(CLINE(1:1)).EQ.32.OR.
     '      IBEG.EQ.1.AND.IEND.EQ.1.AND.CLINE(1:1).EQ.' ') THEN !Return key pressed only
            NCHAR=0
          ELSE
            NCHAR=1+IEND-IBEG
          ENDIF
        ENDIF
      ELSE IF(USE_SOCKET) THEN
        IF(FSKREAD(CLEN,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
        IF(FSKREAD(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1) GOTO 9999
        CALL FSKC2F(CLINE,CLEN,INTSTR)
        IOSTAT=0
        NCHAR=CLEN
      ENDIF

      CALL EXITS('READ_LINE')
      RETURN
 9999 CALL ERRORS('READ_LINE',ERROR)
      CALL EXITS('READ_LINE')
      RETURN 1
      END


