      SUBROUTINE WRITE_LINE(IUNIT,LINE,ERROR,*)

C#### Subroutine: WRITE_LINE
C###  Description:
C###    Simple interface to WRITE_STRING that allows writing of a fortran
C###    character variable without supplying a length and adds a newline.
C###  See-Also: WRITE_STRING, WRITE_CHAR

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IUNIT
      CHARACTER LINE*(*),ERROR*(*)
C     Local Variables
      INTEGER ERR

C**** Note: WRITE_LINE cannot call Enters or Exits (since would get into
C**** an infinite loop).

      CALL WRITE_STRING(IUNIT,LEN(LINE),LINE,ERR)
      IF(ERR.EQ.0.OR.IUNIT.EQ.IOER) THEN
        CALL WRITE_STRING(IUNIT,1,NEWLINE,ERR)
      ENDIF
      IF(ERR.NE.0) THEN
        CALL ERRORIN('WRITE_LINE')
        ERROR=' '
        RETURN 1
      ENDIF
      END


