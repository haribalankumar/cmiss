      SUBROUTINE SHIFT3(AMOUNT,TRANS)

C#### Subroutine: SHIFT3
C###  Description:
C###    SHIFT3 updates the transformation matrix TRANS with a shifting
C###    by AMOUNT along the Z(3) axis.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,TRANS(3,4)

C     CALL ENTERS('SHIFT3',*9999)
      IF(amount.NE.0.0D0) THEN
        TRANS(3,4)=TRANS(3,4)+AMOUNT
      ENDIF

C     CALL EXITS('SHIFT3')
      RETURN
      END


