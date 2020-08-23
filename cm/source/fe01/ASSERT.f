      SUBROUTINE ASSERT(LOEXPR,ERRMSG,ERROR,*)

C#### Subroutine: ASSERT
C###  Description:
C###    ASSERT assigns the error message ERRMSG to ERROR if the
C###    logical expression LOEXPR is false, and control is returned
C###    via the alternative return argument.

      IMPLICIT NONE
!     Parameter List
      CHARACTER ERRMSG*(*),ERROR*(*)
      LOGICAL LOEXPR

      CALL ENTERS('ASSERT',*9999)
      IF(LOEXPR) THEN
        ERROR=' '
      ELSE
        ERROR=ERRMSG
        GO TO 9999
      ENDIF

      CALL EXITS('ASSERT')
      RETURN
 9999 CALL EXITS('ASSERT')
      RETURN 1
      END


