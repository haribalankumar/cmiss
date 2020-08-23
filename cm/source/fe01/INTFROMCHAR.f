      SUBROUTINE INTFROMCHAR(IRESULT,CDATA,ERR)

C#### Subroutine: INTFROMCHAR
C###  Description:
C###    Converts character string CDATA to an integer.
C###    IDATA is calculated as a reasonable guess of an integer
C###    from the leading integer characters in CDATA.  If there are no
C###    leading integer characters IRESULT is set to 0.  If the whole of
C###    CDATA does not look like an integer then an error message is
C###    reported and ERR is set to a non-zero value.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER (ROUTINENAME='INTFROMCHAR')
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IRESULT !out
      CHARACTER CDATA*(*) !in
      INTEGER ERR !out
!     Local Variables
      INTEGER i,NUM_SIGNS,SIGN
      LOGICAL CONTINUE

C     CALL ENTERS(ROUTINENAME,*9999)
      ERR=0
      NUM_SIGNS=0
      SIGN=1
C     Get sign and first digit
      i=0
      CONTINUE=.TRUE.
      DO WHILE(i.LT.LEN(CDATA).AND.CONTINUE)
        i=i+1
        IF(CDATA(i:i).GE.'0'.AND.CDATA(i:i).LE.'9') THEN
          IRESULT=ICHAR(CDATA(i:i))-ICHAR('0')
          CONTINUE=.FALSE.
        ELSEIF(CDATA(i:i).EQ.'-') THEN
          NUM_SIGNS=NUM_SIGNS+1
          SIGN=-SIGN
        ELSEIF(CDATA(i:i).EQ.'+') THEN
          NUM_SIGNS=NUM_SIGNS+1
        ELSEIF(CDATA(i:i).EQ.' ') THEN
          ! do nothing (this allows blanks between the sign and the digits)
        ELSE !unexpected character
          IRESULT=0
          ERR=1
          CONTINUE=.FALSE.
        ENDIF
      ENDDO !i

C     Get digits
      CONTINUE=.TRUE.
      DO WHILE(i.LT.LEN(CDATA).AND.CONTINUE)
        i=i+1
        IF(CDATA(i:i).GE.'0'.AND.CDATA(i:i).LE.'9') THEN
          IRESULT=IRESULT*10+ICHAR(CDATA(i:i))-ICHAR('0')
        ELSEIF(CDATA(i:i).EQ.' ') THEN
          CONTINUE=.FALSE. !no blanks allowed between digits
        ELSE
          ERR=1
          CONTINUE=.FALSE.
        ENDIF
      ENDDO !i
          
C     There should be only one sign.
      IF(NUM_SIGNS.GT.1) ERR=1
      IRESULT=IRESULT*SIGN

C     The remainder should be blanks
      DO WHILE(i.LT.LEN(CDATA).AND.ERR.EQ.0)
        i=i+1
        IF(CDATA(i:i).NE.' ') ERR=1
      ENDDO !i

      IF(ERR.NE.0) THEN
        CALL FLAG_ERROR(0,' ')
        CALL WRITE_CHAR(IOER,'cannot convert ''',ERR)
        CALL WRITE_CHAR(IOER,CDATA,ERR)
        CALL WRITE_CHAR(IOER,''' to an integer'//NEWLINE,ERR)
        CALL ERRORIN(ROUTINENAME)
      ENDIF

C     CALL EXITS(ROUTINENAME)
      RETURN
      END
