      SUBROUTINE GETNEXTCOMFILELINE(CSTRMAXLEN,CSTRING,ERR)

C#### Subroutine: GETNEXTCOMFILELINE
C###  Description:
C###    Returns the next record from the current command file as a
C###    C string.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'gtstr00.cmn'
!     Parameter List
      INTEGER CSTRMAXLEN,CSTRING(*),ERR
!     Local Variables
      INTEGER IOSTAT,MINLEN
      CHARACTER COMMAND*255

C      DATA COMMAND / ' ' /

      IF(IREC_COMFILE(0).LT.0) THEN !no comfile opened
        ERR=1
        RETURN
      ELSE
        IREC_COMFILE(0)=IREC_COMFILE(0)+1
C KAT 1Feb00: END= not available with REC=
C        READ(OPEN_COM_UNIT,FMT='(A)',REC=IREC_COMFILE(0),
C     '    IOSTAT=IOSTAT,END=9998) COMMAND
        READ(OPEN_COM_UNIT,FMT='(A)',REC=IREC_COMFILE(0),
     '    IOSTAT=IOSTAT) COMMAND
        IF(IOSTAT.EQ.0) THEN
          ERR=0
          MINLEN=MIN(LEN(COMMAND),CSTRMAXLEN)
          CALL F2CSTRING(CSTRING,COMMAND(1:MINLEN))
        ELSE IF(IOSTAT.LE.0) THEN !end of file
          IREC_COMFILE(0)=IREC_COMFILE(0)-1
C         If previous line is readable, leave IREC_COMFILE(0) after that line
          READ(OPEN_COM_UNIT,FMT='(A)',REC=IREC_COMFILE(0),
     '      IOSTAT=IOSTAT) COMMAND
          IF(IOSTAT.EQ.0) IREC_COMFILE(0)=IREC_COMFILE(0)+1
          ERR=1
        ELSE
          ERR=-1
        ENDIF
      ENDIF
      RETURN
C KAT 1Feb00: END= not available with REC=
C 9998 IREC_COMFILE(0)=IREC_COMFILE(0)-1
C      ERR=1
C      RETURN
      END


