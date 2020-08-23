      LOGICAL FUNCTION LOWCAS(CHAR)

C#### Function: LOWCAS
C###  Type: LOGICAL
C###  Description:
C###    LOWCAS returns .TRUE. if the character CHAR is lower
C###    case alphabetic.

      IMPLICIT NONE
!     Parameter List
      CHARACTER CHAR

      IF(LGE(CHAR,'a').AND.LLE(CHAR,'z')) THEN
        LOWCAS=.TRUE.
      ELSE
        LOWCAS=.FALSE.
      ENDIF
      RETURN
      END


