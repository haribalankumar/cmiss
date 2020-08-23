      SUBROUTINE C2FSTRING(CSTRING,CLENGTH,FSTRING)

C#### Subroutine: C2FSTRING
C###  Description: Converts a c (integer) string to a fortran string.

      IMPLICIT NONE
!     Parameter List
      INTEGER CSTRING(*),CLENGTH
      CHARACTER FSTRING*(*)
!     Local Variables
      INTEGER I,INTCHAR

      FSTRING=' '
      DO I=1,MIN(CLENGTH,LEN(FSTRING))
        CALL UNPACKCHARACTERS(INTCHAR,I-1,CSTRING)
        FSTRING(I:I)=CHAR(INTCHAR)
      ENDDO !I

      RETURN
      END


