      SUBROUTINE SHIFT1(AMOUNT,TRANS)

C#### Subroutine: SHIFT1
C###  Description:
C###    SHIFT1 updates the transformation matrix TRANS with a shifting
C###    by AMOUNT along the Z(1) axis.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,TRANS(3,4)

C     CALL ENTERS('SHIFT1',*9999)
      IF(amount.NE.0.0D0) THEN
        TRANS(1,4)=TRANS(1,4)+AMOUNT
      ENDIF

C     CALL EXITS('SHIFT1')
      RETURN
      END


