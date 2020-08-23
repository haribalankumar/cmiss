      SUBROUTINE APPENDCA(IEND,LENSTR2,CHAR1,STR2,ERR)

C#### Subroutine: APPENDC
C###  Description:
C###    APPENDC writes fortran character CHAR1 into array of characters
C###    STR2 after position IEND and sets IEND to the last copied
C###    character in STRING2.

      IMPLICIT NONE
!     Parameter List
      INTEGER IEND,LENSTR2
      CHARACTER CHAR1*(*),STR2(LENSTR2)
      INTEGER ERR
!     Local Variables
      INTEGER i,LENCHAR1

C      CALL ENTERS('APPENDC',*9999)

      LENCHAR1=LEN(CHAR1)
      IF(IEND+LENCHAR1.GT.LENSTR2) THEN
        ERR=1
      ELSE
        DO i=1,LENCHAR1
          IEND=IEND+1
          STR2(i)=CHAR1(i:i)
        ENDDO
        ERR=0
      ENDIF

C      CALL EXITS('APPENDC')
      RETURN
      END


