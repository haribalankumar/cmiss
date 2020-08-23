      SUBROUTINE UNKNOW(noco,CO,ERROR,*)

C#### Subroutine: UNKNOW
C###  Description:
C###    UNKNOW gives error message if command is ambiguous or unknown.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER noco
      CHARACTER CO(*)*(*),ERROR*(*)
!     Local Variables
      INTEGER ERR,IBEG,IEND

      CALL ENTERS('UNKNOW',*9998)
      IF(CO(noco).NE.' ') THEN
        CALL STRING_TRIM(CO(noco),IBEG,IEND)
        CALL FLAG_ERROR(0,' ')
        CALL WRITE_CHAR(IOER,'Command "',ERR)
        CALL WRITE_CHAR(IOER,CO(noco)(IBEG:IEND),ERR)
        CALL WRITE_CHAR(IOER,'" is ambiguous or unknown.',ERR)
        CALL WRITE_CHAR(IOER,NEWLINE,ERR)
        GOTO 9998
      ENDIF

      CALL EXITS('UNKNOW')
      RETURN
 9998 ERROR=' '
      CALL ERRORS('UNKNOW',ERROR)
      CALL EXITS('UNKNOW')
      RETURN 1
      END


