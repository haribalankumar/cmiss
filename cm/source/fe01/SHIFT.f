      SUBROUTINE SHIFT(AXIS,AMOUNT,TRANS,ERROR,*)

C#### Subroutine: SHIFT
C###  Description:
C###    SHIFT updates the transformation matrix TRANS with a shifting
C###    by AMOUNT along the axis defined by AXIS.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,AXIS(3),TRANS(3,4)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 AMOUN0

      CALL ENTERS('SHIFT',*9999)
      IF(AXIS(1).NE.0.0D0) THEN
        IF(AXIS(2).NE.0.0D0) THEN
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN0=AMOUNT/DSQRT(AXIS(1)**2+AXIS(2)**2+AXIS(3)**2)
            CALL SHIFT1(AMOUN0*AXIS(1),TRANS)
            CALL SHIFT2(AMOUN0*AXIS(2),TRANS)
            CALL SHIFT3(AMOUN0*AXIS(3),TRANS)
          ELSE
            AMOUN0=AMOUNT/DSQRT(AXIS(1)**2+AXIS(2)**2)
            CALL SHIFT1(AMOUN0*AXIS(1),TRANS)
            CALL SHIFT2(AMOUN0*AXIS(2),TRANS)
          ENDIF
        ELSE
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN0=AMOUNT/DSQRT(AXIS(1)**2+AXIS(3)**2)
            CALL SHIFT1(AMOUN0*AXIS(1),TRANS)
            CALL SHIFT3(AMOUN0*AXIS(3),TRANS)
          ELSE
            CALL SHIFT1(SIGN(AMOUNT,AMOUNT*AXIS(1)),TRANS)
          ENDIF
        ENDIF
      ELSE
        IF(AXIS(2).NE.0.0D0) THEN
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN0=AMOUNT/DSQRT(AXIS(2)**2+AXIS(3)**2)
            CALL SHIFT2(AMOUN0*AXIS(2),TRANS)
            CALL SHIFT3(AMOUN0*AXIS(3),TRANS)
          ELSE
            CALL SHIFT2(SIGN(AMOUNT,AMOUNT*AXIS(2)),TRANS)
          ENDIF
        ELSE
          IF(AXIS(3).NE.0.0D0) THEN
            CALL SHIFT3(SIGN(AMOUNT,AMOUNT*AXIS(3)),TRANS)
          ELSE
C PJH 9/7/98  No error needed - should be able to shift about 0,0,0 point
C             ERROR='The axis of shift is undefined'
C             GO TO 9999
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('SHIFT')
      RETURN
 9999 CALL ERRORS('SHIFT',ERROR)
      CALL EXITS('SHIFT')
      RETURN 1
      END


