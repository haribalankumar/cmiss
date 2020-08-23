      SUBROUTINE BINOPENFILE(FILEID,COMMAND,FILENAME,ERROR,*)

C#### Subroutine: BINOPENFILE
C###  Description:
C###    BINOPENFILE opens a binary file specified by FILEID (> 0
C###    and < 99) of the given filename for either reading
C###    or writing (specified by command being "READ" or "WRITE").

      IMPLICIT NONE
!     Parameter List
      INTEGER FILEID
      CHARACTER COMMAND*(*),FILENAME*(*),ERROR*(*)
!     Local Variables
      INTEGER CACCESSCODE(2),CERROR(50),CERRLEN,CFNAME(50),ERR
      CHARACTER FACCESSCODE*6

      CALL ENTERS('BINOPENFILE',*9999)

      CALL F2CSTRING(CFNAME,FILENAME)
      IF(COMMAND(1:4).EQ.'READ') THEN
C Why was this set to read/write? Will not work for user
C running example, only CMISS will have permissions to write
C to file.
C         FACCESSCODE='rb+'
        FACCESSCODE='rb'
      ELSE IF(COMMAND(1:5).EQ.'WRITE') THEN
        FACCESSCODE='wb+'
      ELSE
        ERROR='>>Invalid command'
        GOTO 9999
      ENDIF
      CALL F2CSTRING(CACCESSCODE,FACCESSCODE)
      CALL BINARYOPENFILE(FILEID,CFNAME,CACCESSCODE,ERR,CERROR)
      IF(ERR.NE.0) THEN
        CALL CSTRINGLEN(CERRLEN,CERROR)
        CALL C2FSTRING(CERROR,CERRLEN,ERROR)
        GOTO 9999
      ENDIF

      CALL EXITS('BINOPENFILE')
      RETURN
 9999 CALL ERRORS('BINOPENFILE',ERROR)
      CALL EXITS('BINOPENFILE')
      CALL BINARYCLOSEFILE(FILEID,ERR,CERROR)
      RETURN 1
      END


