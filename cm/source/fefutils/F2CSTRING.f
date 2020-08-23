      SUBROUTINE F2CSTRING(CSTRING,FSTRING)

C#### Subroutine: F2CSTRING
C###  Description: Converts a fortran string to a c (integer) string.

      IMPLICIT NONE
!     Parameter List
      INTEGER CSTRING(*)
      CHARACTER FSTRING*(*)
!     Local Variables
      INTEGER FSTRINGLEN,I,INTCHAR,LENGTH

      LENGTH=FSTRINGLEN(FSTRING)
      DO I=1,LENGTH
        INTCHAR=ICHAR(FSTRING(I:I))
        CALL PACKCHARACTERS(INTCHAR,I-1,CSTRING)
      ENDDO !I

      INTCHAR=0
      CALL PACKCHARACTERS(INTCHAR,I-1,CSTRING)

      RETURN
      END
