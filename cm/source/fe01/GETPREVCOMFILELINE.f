      SUBROUTINE GETPREVCOMFILELINE(CSTRMAXLEN,CSTRING,ERR)

C#### Subroutine: GETPREVCOMFILELINE
C###  Description:
C###    Returns the previous record from the current command file as a
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

      IF(IREC_COMFILE(0).LE.0) THEN
        ERR=1
        RETURN
      ELSE
        IREC_COMFILE(0)=IREC_COMFILE(0)-1
        READ(OPEN_COM_UNIT,FMT='(A)',REC=IREC_COMFILE(0),
     '    IOSTAT=IOSTAT) COMMAND
        IF(IOSTAT.EQ.0) THEN
          ERR=0
          MINLEN=MIN(LEN(COMMAND),CSTRMAXLEN)
          CALL F2CSTRING(CSTRING,COMMAND(1:MINLEN))
        ELSE
          ERR=-1
        ENDIF
      ENDIF
      RETURN
      END


