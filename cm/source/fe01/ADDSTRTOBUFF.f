      SUBROUTINE ADDSTRTOBUFF(STRING,ERROR,*)

C#### Subroutine: ADDSTRTOBUFF
C###  Description:
C###    ADDSTRTOBUFF adds the command string STRING to the buffer of
C###    command strings

      IMPLICIT NONE
c      INCLUDE 'cmiss$reference:gtstr00.cmn'
!     Parameter List
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER CSTRING(100)

      CALL ENTERS('ADDSTRTOBUFF',*9999)

      CALL F2CSTRING(CSTRING,STRING)
      CALL ADDSTRINGTOBUFFER(CSTRING)

      CALL EXITS('ADDSTRTOBUFF')
      RETURN
 9999 CALL ERRORS('ADDSTRTOBUFF',ERROR)
      CALL EXITS('ADDSTRTOBUFF')
      RETURN 1
      END


