      SUBROUTINE WRITE_STRING(IUNIT,STRLEN,STRING,ERR)

C#### Subroutine: WRITE_STRING
C###  Description:
C###    <HTML>
C###    Writes the first STRLEN characters from STRING to output
C###    identified by IUNIT as follows:
C###    <PRE>
C###      IOOP listing
C###      IODI diagnostics
C###      IOTR trace
C###      IOER errors
C###      IOH1 ? help
C###      IOH2 ?? help
C###      IOH3 ??? help
C###
C###    IUNIT>9 is used for file output
C###    </PRE>
C###    ERR is 0 on success and 1 on failure
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'error0.inc'
      INCLUDE 'fsklib.inc'
!     Parameter List
      INTEGER IUNIT,STRLEN
      CHARACTER STRING(STRLEN)
      INTEGER ERR
!     Local Variables
      INTEGER CODE,ERR2,i,IBEG,IEND,j
C     Socket stuff
      INTEGER DATA_TYPE
      PARAMETER(DATA_TYPE = 0)  !Cmiss Output
      INTEGER CLEN
      CHARACTER C_STRING(1024)

C**** Note: cannot call Enters or Exits (since would get into
C**** an infinite loop).

      IF(IUNIT.GE.IOOUT.OR..NOT.(USE_SOCKET.OR.CMGUI_LINK)) THEN
        ! file output or similar

C       xlf on IRIX 5.1 has a limit (about 32768) on the number of
C       characters in a record when writing to a pipe (??&!!).
C       A new record apparently is only started when the output format
C       contains no $ descriptor.
C       Therefore start a new record when it seems appropriate.
C       Starting a new record on each newline should also convert the
C       newline character to the platform's line separator.
        i=1
        IBEG=1
        IEND=MIN(9999,STRLEN)
        DO WHILE(i.LE.STRLEN)
          IF(STRING(i).EQ.NEWLINE) THEN
            IEND=i-1
          ENDIF
          IF(i.GE.IEND) THEN
            WRITE(IUNIT,'($,9999A)') (STRING(j),j=IBEG,IEND)
            IF(i.GT.IEND) WRITE(IUNIT,'()')
            IBEG=i+1
            IEND=MIN(i+9999,STRLEN)
          ENDIF
          i=i+1
        ENDDO !i
        ERR=0
C       This simpler code works but the above seems a little faster on AIX
C         i=1
C         DO WHILE(i.LE.STRLEN)
C           IF(STRING(i).EQ.NEWLINE) THEN
C             WRITE(IUNIT,'()')
C           ELSE
C             WRITE(IUNIT,'($,A)') STRING(i)
C           ENDIF
C           i=i+1
C         ENDDO !i

      ELSE !socket or cmgui link

        ERR=1 ! until succeeded
        IF(USE_SOCKET) THEN !transfer string to frontend thru socket
          IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1)
     '      GOTO 999
C!!! This seems inconsistent with DISPLAY_PROMPT
          IF(FSKWRITE(IUNIT,SK_LONG_INT,1,CONNID2).EQ.-1)
     '      GOTO 999
          CLEN=MIN(STRLEN,1023)
          DO i=1,CLEN
            C_STRING(i)=STRING(i)
          ENDDO
          CLEN=CLEN+1
          C_STRING(CLEN)=CHAR(0)
          IF(FSKWRITE(CLEN,SK_LONG_INT,1,CONNID2).EQ.-1)
     '      GOTO 999
          IF(FSKWRITE(C_STRING,SK_CHAR,CLEN,CONNID2).EQ.-1)
     '      GOTO 999

        ELSE ! IF(CMGUI_LINK) THEN !write out directly
          CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_COMMAND_I,0,2,CODE)
          IF(CODE.EQ.0) GOTO 999
C??? What does this int indicate?
C??? Should it be different for error/informational messages?
          CALL WH_INPUT_F_ADD_INT(CMGUI_COMMAND_I,1,1,CODE)
          IF(CODE.EQ.0) GOTO 999
C         Check for a blank string
          IF(STRLEN.GT.0) THEN !STRLEN may be 0 in f90
            CALL WH_INPUT_F_ADD_CHAR(CMGUI_COMMAND_I,STRLEN,STRING,
     '        CODE)
          ELSE !send a space as a separator
            CALL WH_INPUT_F_ADD_CHAR(CMGUI_COMMAND_I,1,' ',CODE)
          ENDIF !zero length
          IF(CODE.EQ.0) GOTO 999
          CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_COMMAND_I,CODE)
          IF(CODE.EQ.0) GOTO 999
C cpb 4/3/97 Update command window with output.
C Only update the command window if the text seems complete.
C i.e. ends in a newline.
          IF(STRING(STRLEN).EQ.NEWLINE) THEN
            CALL WH_INPUT_F_UPDATE(CMGUI_COMMAND_I,CODE)
            IF(CODE.EQ.0) GOTO 999
          ENDIF

        ENDIF !socket/cmgui_link
        ERR=0

 999    CONTINUE ! ERR=1
        IF(ERR.NE.0) THEN
          IF(IUNIT.EQ.IOER) THEN
C           Try to ensure error messages are sent somewhere
            WRITE(*,*) STRING
          ELSE
            CALL ERRORIN('WRITE_STRING')
          ENDIF
        ENDIF

      ENDIF

C     Output strings to echo output file if required
      IF(ECHO_OUTPUT.AND.IUNIT.LT.IOOUT) THEN
        CALL WRITE_STRING_WRAPPER(IOOUT,STRLEN,STRING,ERR2)
        ! Currently this can't return nonzero ERR2,
        ! so not bothering to check error code.
      ENDIF
        
      END


