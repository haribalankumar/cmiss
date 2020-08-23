      SUBROUTINE SUBSTR(STRLEN,STRING,DE,ERROR,*)

C#### Subroutine: SUBSTR
C###  Description:
C###    SUBSTR breaks the string STRING into substrings SB(nosb)
C###    separated by single character delimiters found in the string DE.
C**** Portions of STRING surrounded by " are not broken into substrings.
C**** Note: the first character of the resulting substrings are delimiters.
C**** STRING is intent(in) and should not be modified.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'trac00.cmn'
!     Parameter List
      INTEGER STRLEN
      CHARACTER STRING(STRLEN),DE*(*),ERROR*(*)
!     Local Variables
      INTEGER ERR,i,IBEG,IEND,n1,NODE
      CHARACTER CHAR1,DELIMITER
      LOGICAL CELEM,ST

      CALL ENTERS('SUBSTR',*9999)
      NTSB=0
      IBEG=1
      IEND=STRLEN
      DO WHILE(STRING(IEND).EQ.' '.AND.IBEG.LT.IEND)
        IEND=IEND-1
      ENDDO !i
      DO WHILE(STRING(IBEG).EQ.' '.AND.IBEG.LT.IEND)
        IBEG=IBEG+1
      ENDDO !i
      IF(STRING(IBEG).NE.' ') THEN
        NODE=IBEG
        ST=.FALSE. !not in quoted string
        DELIMITER=' '
        DO n1=IBEG,IEND+1
          IF(n1.EQ.IEND+1) THEN
            CHAR1=' '
          ELSE
            CHAR1=STRING(n1)
          ENDIF
          IF(n1-NODE.GT.LEN(SB(1))) THEN
            CALL FLAG_ERROR(-1,
     '        'Word in interpreted command is too long:')
            CALL WRITE_STRING(IOER,n1-NODE,STRING(NODE),ERR)
            CALL WRITE_STRING(IOER,1,NEWLINE,ERR)
            GOTO 9999
          ENDIF
          IF(ST) THEN
            IF(CHAR1.EQ.'"') THEN
              NTSB=NTSB+1
              SB(NTSB)=DELIMITER
              DO i=0,n1-1-NODE
                SB(NTSB)(2+i:2+i)=STRING(NODE+i)
              ENDDO !i
              DELIMITER=' '
              NODE=n1+1
              ST=.FALSE.
              IF(NTSB.EQ.NXSB) THEN
                GOTO 2
              ENDIF
            ENDIF
          ELSE
            IF(CHAR1.EQ.CHAR(9)) CHAR1=' ' ! tabs -> spaces
            IF(CELEM(CHAR1,DE)) THEN
              IF(n1.EQ.NODE) THEN
C               spaces may surround a delimiter
                IF(CHAR1.EQ.' ') THEN
                ELSEIF(DELIMITER.EQ.' ') THEN
                  DELIMITER=CHAR1 !use non space delimiter
                ELSE
                  NTSB=NTSB+1
                  SB(NTSB)=DELIMITER
                  DELIMITER=CHAR1
                  IF(NTSB.EQ.NXSB) GOTO 2
                ENDIF
              ELSE
                NTSB=NTSB+1
                SB(NTSB)=DELIMITER
                DO i=0,n1-1-NODE
                  SB(NTSB)(2+i:2+i)=STRING(NODE+i)
                ENDDO !i
                DELIMITER=CHAR1
                IF(NTSB.EQ.NXSB) THEN
                  GOTO 2
                ENDIF
              ENDIF
              NODE=n1+1
            ELSEIF(n1.EQ.NODE.AND.CHAR1.EQ.'"') THEN
C             Start of quote
              NODE=n1+1
              ST=.TRUE.
            ENDIF
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('SUBSTR')
      RETURN
 2    CALL FLAG_ERROR(-1,'Too many words in command')
 9999 CALL ERRORIN('SUBSTR')
      ERROR=' '
      CALL EXITS('SUBSTR')
      RETURN 1
      END


