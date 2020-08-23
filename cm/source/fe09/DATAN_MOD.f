      REAL*8 FUNCTION DATAN_MOD(X,Y)

C#### Function: DATAN_MOD
C###  Type: REAL*8
C###  Description:
C###    DATAN_MOD returns a modified value of atan which takes
C###    into account if X or Y is zero.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
!     Parameter List
      REAL*8 X,Y
!     Local Variables
      REAL*8 RESULT

      IF(X.EQ.0.0d0) THEN
        IF(Y.EQ.0.0d0) THEN
          RESULT=0.0d0
        ELSE IF(Y.LT.0.0d0) THEN
          RESULT=-PI/2.0d0
        ELSE
          RESULT=PI/2.0d0
        ENDIF
      ELSE IF(Y.EQ.0.0d0) THEN
        IF(X.LE.0.0d0) THEN
          RESULT=PI
        ELSE
          RESULT=0.0d0
        ENDIF
      ELSE
        RESULT=DATAN2(Y,X)
      ENDIF

      DATAN_MOD=RESULT

      RETURN
      END


