      INTEGER FUNCTION CM_FLOOR(X)

C#### Function: CM_FLOOR
C###  Type: INTEGER
C###  Description:
C###    Finds the CM_FLOOR of a REAL*8 variable
CC 14-MAR-2000

      IMPLICIT NONE
!     Parameter List
      REAL*8 X

      IF(X.GT.0.0d0) THEN
        CM_FLOOR=INT(X)
      ELSE
        CM_FLOOR=INT(X)-1
      ENDIF

      RETURN
      END

