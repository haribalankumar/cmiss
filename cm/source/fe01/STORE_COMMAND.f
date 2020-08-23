      SUBROUTINE STORE_COMMAND(COMMAND_LINE,ERR,ERROR,*)

C#### Subroutine: STORE_COMMAND
C###  Description:
C###    Adds the given COMMAND_LINE to the command buffer and
C###    learn command files.
C***  KAT 28/3/00

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'learn00.cmn'

!     Parameter List
      CHARACTER*(*) COMMAND_LINE
      INTEGER ERR
      CHARACTER*(*) ERROR

      CALL ENTERS('STORE_COMMAND',*998)

      !Add command to command buffer
      CALL ADDSTRTOBUFF(COMMAND_LINE,ERROR,*998)
      !Send command to learn file if appropriate
      IF(LEARN) THEN
        WRITE(UNIT=COM_UNIT,FMT='(A)',IOSTAT=ERR) COMMAND_LINE
        IF(ERR.NE.0) THEN
          ERROR='Error writing to learn file'
          GOTO 999
        ENDIF !ERR
      ENDIF !LEARN

      CALL EXITS('STORE_COMMAND')
      RETURN
 998  ERR=-1
 999  CALL ERRORS('STORE_COMMAND',ERROR)
      CALL EXITS('STORE_COMMAND')
      RETURN
      END


