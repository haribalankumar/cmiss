      INTEGER FUNCTION CLOCAT(SUBSTR,STRING)

C#### Function: CLOCAT
C###  Type: INTEGER
C###  Description:
C###    CLOCAT returns the position of the first occurrence of SUBSTR
C###    within STRING.  If SUBSTR is not found within STRING, CLOCAT
C###    returns zero.

      IMPLICIT NONE
!     Parameter List
      CHARACTER SUBSTR*(*),STRING*(*)
!     Local Variables
      INTEGER l0,l1,LEN,LENSUB,LENSTR

      LENSUB=LEN(SUBSTR)
      LENSTR=LEN(STRING)
      l0=LENSUB-1
      CLOCAT=0
C KAT 13Oct97: changed to check last character in STRING
C      DO l1=1,(LENSTR-LENSUB)
      DO l1=1,(LENSTR-l0)
        IF(SUBSTR.EQ.STRING(l1:l1+l0)) THEN
          CLOCAT=l1
          RETURN
        ENDIF
      ENDDO

      RETURN
      END


