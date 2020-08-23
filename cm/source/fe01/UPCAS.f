      LOGICAL FUNCTION UPCAS(CHAR)

C#### Function: UPCAS
C###  Type: LOGICAL
C###  Description:
C###    UPCAS returns .TRUE. if the character CHAR is upper case
C###    alphabetic.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CHAR

      IF(LGE(CHAR,'A').AND.LLE(CHAR,'Z')) THEN
        UPCAS=.TRUE.
      ELSE
        UPCAS=.FALSE.
      ENDIF
      RETURN
      END


