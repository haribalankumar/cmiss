      LOGICAL FUNCTION ISENDBINFILE(FILEID)

C#### Function: ISENDBINFILE
C###  Type: LOGICAL
C###  Description:
C###    ISENDBINFILE returns .TRUE. if the binary file specified by
C###    fileid is at end of file, .FALSE. if not.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER FILEID
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR,RETURNCODE
      CHARACTER ERROR*100

      CALL ISENDBINARYFILE(FILEID,RETURNCODE,ERR,ERROR)
      IF(RETURNCODE.EQ.1) THEN
        ISENDBINFILE=.TRUE.
      ELSE
        ISENDBINFILE=.FALSE.
      ENDIF
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        WRITE(OP_STRING,'(1X,A)') ERROR
      ENDIF

      RETURN
      END


