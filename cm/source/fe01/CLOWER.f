      CHARACTER*(*) FUNCTION CLOWER(STRING)

C#### Function: CLOWER
C###  Type: CHARACTER
C###  Description:
C###    CLOWER converts the character string STRING to lower case.

      IMPLICIT NONE
!     Parameter List
      CHARACTER STRING*(*)
!     Local Variables
      INTEGER i,OFFSET
      LOGICAL UPCAS

      OFFSET=ICHAR('a')-ICHAR('A')
      CLOWER=STRING
      DO i=1,LEN(STRING)
        IF(UPCAS(STRING(i:i))) THEN
          CLOWER(i:i)=CHAR(ICHAR(STRING(i:i))+OFFSET)
        ENDIF
      ENDDO

      RETURN
      END

      
!newe


C old MLB 19-9-96 changed to subroutine
C      CHARACTER*(*) FUNCTION CUPPER(STRING)
C
CC#### Function: CUPPER
CC###  Type: CHARACTER
CC###  Description:
CC###    CUPPER converts the character string STRING to upper case.
C
C      IMPLICIT NONE
C!     Parameter List
C      CHARACTER STRING*(*)
C!     Local Variables
C      INTEGER i,OFFSET
C      LOGICAL LOWCAS
C
C      OFFSET=ICHAR('A')-ICHAR('a')
C      CUPPER=STRING
C      DO i=1,LEN(STRING)
C        IF(LOWCAS(STRING(i:i))) THEN
C          CUPPER(i:i)=CHAR(ICHAR(STRING(i:i))+OFFSET)
C        ENDIF
C      ENDDO
C      RETURN
C      END


