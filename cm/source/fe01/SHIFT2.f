      SUBROUTINE SHIFT2(AMOUNT,TRANS)

C#### Subroutine: SHIFT2
C###  Description:
C###    SHIFT2 updates the transformation matrix TRANS with a shifting
C###    by AMOUNT along the Z(2) axis.

      IMPLICIT NONE
!     Parameter List
      REAL*8 AMOUNT,TRANS(3,4)

C     CALL ENTERS('SHIFT2',*9999)
      IF(amount.NE.0.0D0) THEN
        TRANS(2,4)=TRANS(2,4)+AMOUNT
      ENDIF

C     CALL EXITS('SHIFT2')
      RETURN
      END


