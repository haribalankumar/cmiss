      SUBROUTINE TWIST(AXIS,AMOUNT,TRANS,ERROR,*)

C#### Subroutine: TWIST
C###  Description:
C###    TWIST updates the transformation matrix TRANS with a rotation
C###    of AMOUNT radians about the axis defined by AXIS.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,AXIS(3),TRANS(3,4)
      CHARACTER ERROR*(*)
!     Local Variables
      REAL*8 AMOUN1,AMOUN2,AMOUN3

      CALL ENTERS('TWIST',*9999)
      IF(AXIS(1).NE.0.0D0) THEN
        IF(AXIS(2).NE.0.0D0) THEN
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN1=DATAN2(AXIS(2),AXIS(3))
            AMOUN2=DATAN2(DSQRT(AXIS(2)**2+AXIS(3)**2),AXIS(1))
            CALL TWIST1( AMOUN1,TRANS)
            CALL TWIST2(-AMOUN2,TRANS)
            CALL TWIST3( AMOUNT,TRANS)
            CALL TWIST2( AMOUN2,TRANS)
            CALL TWIST1(-AMOUN1,TRANS)
          ELSE
            AMOUN3=DATAN2(AXIS(2),AXIS(1))
            CALL TWIST3(-AMOUN3,TRANS)
            CALL TWIST1( AMOUNT,TRANS)
            CALL TWIST3( AMOUN3,TRANS)
          ENDIF
        ELSE
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN2=DATAN2(AXIS(1),AXIS(3))
            CALL TWIST2(-AMOUN2,TRANS)
            CALL TWIST3( AMOUNT,TRANS)
            CALL TWIST2( AMOUN2,TRANS)
          ELSE
            CALL TWIST1( SIGN( AMOUNT, AXIS(1)*AMOUNT ),TRANS)
          ENDIF
        ENDIF
      ELSE
        IF(AXIS(2).NE.0.0D0) THEN
          IF(AXIS(3).NE.0.0D0) THEN
            AMOUN1=DATAN2(AXIS(3),AXIS(2))
            CALL TWIST1( AMOUN1,TRANS)
            CALL TWIST2( AMOUNT,TRANS)
            CALL TWIST1(-AMOUN1,TRANS)
          ELSE
            CALL TWIST2( SIGN( AMOUNT, AXIS(2)*AMOUNT ),TRANS)
          ENDIF
        ELSE
          IF(AXIS(3).NE.0.0D0) THEN
            CALL TWIST3( SIGN( AMOUNT, AXIS(3)*AMOUNT ),TRANS)
          ELSE
            ERROR='The axis of rotation is undefined'
            GO TO 9999
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('TWIST')
      RETURN
 9999 CALL ERRORS('TWIST',ERROR)
      CALL EXITS('TWIST')
      RETURN 1
      END


