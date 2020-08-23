      SUBROUTINE BINSKIPCMHEADER(FILEID,IP,ERROR,*)

C#### Subroutine: BINSKIPCMHEADER
C###  Description:
C###    BINSKIPCMHEADER skips CMISS binary files header information.
C###    If IP is 1 then only the identity section is skiped, if IP
C###    is 2 then the identity and machine header sections are skiped,
C###    and if IP is 3 then all header sections (identity, machine,
C###    and file) are skipped.

      IMPLICIT NONE
      INCLUDE 'binf00.cmn'
      INCLUDE 'mach00.cmn'
      INCLUDE 'mach00.inc'
!     Parameter List
      INTEGER FILEID,IP
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER NUMSKIPBYTES,BINVERSION

      CALL ENTERS('BINSKIPCMHEADER',*9999)

      IF(IP.GE.1) THEN
C***    Skip the identity header section
        CALL BINREADFILE(FILEID,CHARTYPE,2,INTDATA,REAL4DATA,
     '    REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
        CALL ASSERT(CHARDATA(1:1).EQ.CHAR(7),
     '    '>>Not a CMISS binary file',ERROR,*9999)
      ENDIF
      IF(IP.GE.2) THEN
C***    Skip the machine header section
        BINVERSION=ICHAR(CHARDATA(2:2))
        IF(BINVERSION.EQ.0) THEN
          NUMSKIPBYTES=3
        ELSE IF(BINVERSION.EQ.1.OR.
     '      BINVERSION.EQ.2) THEN
          CALL BINREADFILE(FILEID,CHARTYPE,1,INTDATA,REAL4DATA,
     '      REAL8DATA,CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
          NUMSKIPBYTES=ICHAR(CHARDATA(1:1))
        ELSE
          ERROR='>>Unknown binary file header identity format'
          GOTO 9999
        ENDIF
        CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
      ENDIF
      IF(IP.GE.3) THEN
C***    Skip the file header section
        NUMSKIPBYTES=INTSIZE
        IF(BINVERSION.EQ.2) THEN
          NUMSKIPBYTES=NUMSKIPBYTES+3*INTSIZE
        ELSE
          NUMSKIPBYTES=NUMSKIPBYTES+SPSIZE
        ENDIF
        CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
        CALL BINREADFILE(FILEID,INTTYPE,1,INTDATA,REAL4DATA,REAL8DATA,
     '    CHARDATA,LOGDATA,SINTDATA,ERROR,*9999)
        NUMSKIPBYTES=INTDATA(1)*CHARSIZE+INTSIZE
        CALL BINSKIPFILE(FILEID,NUMSKIPBYTES,ERROR,*9999)
      ENDIF

      CALL EXITS('BINSKIPCMHEADER')
      RETURN
 9999 CALL ERRORS('BINSKIPCMHEADER',ERROR)
      CALL EXITS('BINSKIPCMHEADER')
      RETURN 1
      END


