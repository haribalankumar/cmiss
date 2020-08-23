      SUBROUTINE STRING_TRIM(CDATA,IBEGIN,IEND)

C#### Subroutine: STRING_TRIM
C###  Description:
C###   STRING_TRIM finds the first and last non-blank characters in a
C###    string CDATA.
C**** On exit IBEGIN is the position of the first non-blank character
C**** and IEND is the position of the last non-(blank/tab/null)
C**** character.
C**** If CDATA consists of only blank characters IEND will be 1 and
C**** IBEG will be 1.

      IMPLICIT NONE
!     Parameter List
      INTEGER IBEGIN,IEND
      CHARACTER CDATA*(*)
!     Local Variables
      INTEGER l,LENGTH
      CHARACTER NUL*1,TAB*1

      LENGTH=LEN(CDATA)
      NUL=CHAR(0) !ASCII NULL
      TAB=CHAR(9) !ASCII HT (horizontal TAB)
      DO 1 l=1,LENGTH
C MPN 11Mar96 Also ignore tabs and null chars
        IF(CDATA(l:l).NE.' '.AND.CDATA(l:l).NE.TAB.AND.
     '    CDATA(l:l).NE.NUL) THEN
C old        IF(CDATA(l:l).NE.' ') THEN
          IBEGIN=l
          GOTO 2
        ENDIF
 1    CONTINUE
      IBEGIN=LENGTH
 2    DO 3 l=LENGTH,1,-1
C MPN 11Mar96 Also ignore tabs and null chars
        IF(CDATA(l:l).NE.' '.AND.CDATA(l:l).NE.TAB.AND.
     '    CDATA(l:l).NE.NUL) THEN
C old        IF(CDATA(l:l).NE.' ') THEN
          IEND=l
          GOTO 4
        ENDIF
 3    CONTINUE
      IEND=1
 4    IF(IBEGIN.GT.IEND) THEN
        IBEGIN=1
        IEND=1
      ENDIF
      RETURN
      END


