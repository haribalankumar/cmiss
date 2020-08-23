      SUBROUTINE DISPLAY_PROMPT(FORMAT,ERROR,*)

C#### Subroutine: DISPLAY_PROMPT
C###  Description:
C###    DISPLAY_PROMPT displays a prompt with format FORMAT.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='DISPLAY_PROMPT')
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cmis00.cmn'
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'fsklib.inc'
!     Parameter List
      CHARACTER FORMAT*(*),ERROR*(*)
!     Local Variables
      INTEGER DATA_TYPE
      PARAMETER(DATA_TYPE = 1)  !Prompt Data
      INTEGER CLEN,CODE,IBEG,IEND,INTSTR(1024)

      CALL ENTERS(ROUTINENAME,*9999)

      IF(.NOT.USE_SOCKET) THEN
        IF(CMGUI_LINK) THEN
          CALL WH_INPUT_F_OPEN_MESSAGE(CMGUI_PROMPT_I,
     '      0,2,CODE)
          IF(CODE.EQ.0) THEN
            CALL FLAG_ERROR(0,'Could not open prompt message')
            GOTO 9999
          ENDIF
          CALL STRING_TRIM(FORMAT,IBEG,IEND)
C         Check for a blank FORMAT
          IF(IEND.GT.0) THEN
            CALL WH_INPUT_F_ADD_CHAR(CMGUI_PROMPT_I,
     '        IEND,FORMAT,CODE)
          ELSE
C           send a space as a separator
            CALL WH_INPUT_F_ADD_CHAR(CMGUI_PROMPT_I,
     '        1,FORMAT,CODE)
          ENDIF !zero length
          IF(CODE.EQ.0) THEN
            CALL FLAG_ERROR(0,'Could not add string to prompt message')
            GOTO 9999
          ENDIF
          CALL WH_INPUT_F_CLOSE_MESSAGE(CMGUI_PROMPT_I,
     '      CODE)
          IF(CODE.EQ.0) THEN
            CALL FLAG_ERROR(0,'Could not close prompt message')
            GOTO 9999
          ENDIF
        ELSE
          WRITE(UNIT=IOOP,FMT=FORMAT)
        ENDIF
      ELSE IF(USE_SOCKET) THEN
C  Tell the front end it's getting prompt data
        IF(FSKWRITE(DATA_TYPE,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
        CLEN=FSKLEN(FORMAT)
        CALL FSKF2C(FORMAT,CLEN,INTSTR)
        IF(FSKWRITE(CLEN+1,SK_LONG_INT,1,CONNID2).EQ.-1) GOTO 9999
        IF(FSKWRITE(INTSTR,SK_CHAR,CLEN+1,CONNID2).EQ.-1) GOTO 9999
      ENDIF

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORIN(ROUTINENAME)
      ERROR=' '
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


