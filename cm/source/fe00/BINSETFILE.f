      SUBROUTINE BINSETFILE(FILEID,SET_CODE,ERROR,*)

C#### Subroutine: BINSETFILE
C###  Description:
C###    BINSETFILE sets the position of a binary file specified by
C###    fileid. If SET_CODE = 0 the file is positioned at the
C###    beginning of the file, if SET_CODE = 1 the file is positioned
C###    at the current file position and if SET_CODE = 2 then the
C###    file is positioned at the end of the file.

      IMPLICIT NONE
!     Parameter List
      INTEGER FILEID,SET_CODE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER CERROR(50),CERRLEN,ERR

      CALL ENTERS('BINSETFILE',*9999)

      CALL BINARYSETFILE(FILEID,SET_CODE,ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        GOTO 9999
      ENDIF

      CALL EXITS('BINSETFILE')
      RETURN
 9999 CALL ERRORS('BINSETFILE',ERROR)
      CALL EXITS('BINSETFILE')
      CALL BINARYCLOSEFILE(FILEID,ERR,CERROR)
      RETURN 1
      END


