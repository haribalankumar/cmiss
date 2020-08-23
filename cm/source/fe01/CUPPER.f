      SUBROUTINE CUPPER(INSTRING,OUTSTRING)

C#### Subroutine: CUPPER
C###  Description:
C###    CUPPER converts a string to upper case.
C###    It has two arguments (+ error args), the first is the
C###    input string and the second is the output string.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      CHARACTER INSTRING*(*),OUTSTRING*(*)
!     Local Variables
      INTEGER i,MINLEN,OFFSET
      CHARACTER ERROR*100
!     Functions
      LOGICAL LOWCAS

C!!! MHT 30/10/96 DON'T call ENTERS from CUPPER as CUPPER is used within
C!!! ENTERS and will hence creates an infinite loop.
C     CALL ENTERS('CUPPER',*9999)

      OFFSET=ICHAR('A')-ICHAR('a')
C     CALL ASSERT(LEN(INSTRING).LE.LEN(OUTSTRING),'>>Input string'
C    '  //' longer than ouput string',ERROR,*9999)
      IF(LEN(INSTRING).LE.LEN(OUTSTRING))THEN
        MINLEN=LEN(INSTRING)
      ELSE
        WRITE(OP_STRING,'(''>>Error in CUPPER: LEN(OUTSTRING)'
     '    //' < LEN(INSTRING): '',A)') INSTRING
        CALL WRITES(IOER,OP_STRING,ERROR,*9999)
        MINLEN=LEN(OUTSTRING)
      ENDIF
      OUTSTRING=INSTRING
      DO i=1,MINLEN
        IF(LOWCAS(INSTRING(i:i))) THEN
          OUTSTRING(i:i)=CHAR(ICHAR(INSTRING(i:i))+OFFSET)
        ENDIF
      ENDDO

C     CALL EXITS('CUPPER')
C     RETURN
C9999 CALL ERRORS('CUPPER',ERROR)
C     CALL EXITS('CUPPER')

      RETURN
 9999 WRITE(*,'('' Writes Error: '',A)') ERROR
      RETURN
      END


